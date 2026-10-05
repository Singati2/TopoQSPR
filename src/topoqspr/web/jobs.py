"""Persistent job store, background worker and job runner of the web interface.

* :class:`JobStore` - SQLite tables ``cohorts`` and ``jobs`` in ``<data_dir>/web.sqlite``
  (WAL mode).  Every job row survives a restart; :meth:`JobStore.recover` turns rows left
  ``running`` by a previous process into ``failed`` with the reason, and ``queued`` rows
  are picked up again.  Nothing is dropped silently.
* :class:`Worker` - one background thread.  It claims the oldest queued job (at most one
  job runs at a time: the claim is a single SQL statement that fails while another row is
  ``running``) and runs it in a **child process** (``python -m topoqspr.web.jobs JOB_DIR``),
  so a hard timeout can kill it and a crash cannot take the server down.
* :func:`run_job` - the child: runs :func:`topoqspr.pipeline.run_study` on the job's
  ``config.json``, reports the stage and progress through ``progress.json`` (from the
  pipeline's log records), then writes ``ladder.json`` (the result table served by the API)
  and ``results.zip`` (every result file, ``LICENCE-DATA.txt`` and ``README.txt``).
  A failure writes ``error.json`` with the message.

Job directory (``<data_dir>/jobs/<uuid4>/``): ``request.json``, ``config.json``,
``cohort/manifest.csv`` + ``manifest.json``, ``progress.json``, ``worker.log``,
``results/`` (the pipeline's results directory), ``ladder.json``, ``results.zip``,
``error.json``.
"""

from __future__ import annotations

import datetime as _dt
import json
import logging
import math
import os
import shutil
import sqlite3
import subprocess
import sys
import threading
import time
import traceback
import uuid
import zipfile
from pathlib import Path
from typing import Any, Mapping

__all__ = [
    "LADDER_SETS",
    "STAGES",
    "JobStore",
    "Worker",
    "build_ladder",
    "licence_note",
    "new_id",
    "now_utc",
    "run_job",
    "valid_id",
]

#: Feature sets of the web ladder, in display order (the mandatory baselines null, size and
#: composition, the verified indices, and RDKit 2D descriptors; EXPERIMENT_PLAN §5).
LADDER_SETS: tuple[str, ...] = ("null", "size", "composition", "TI", "rdkit2d")

#: Pipeline stages in order, with their share of the progress bar.
STAGES: tuple[tuple[str, float], ...] = (
    ("cohort", 0.02), ("standardize", 0.03), ("features", 0.05), ("targets", 0.10),
    ("evaluate", 0.75), ("protocol", 0.02), ("write", 0.03),
)

STATUSES = ("queued", "running", "done", "failed")

#: Fixed time stamp of every member of ``results.zip`` (byte-stable archives).
_ZIP_DATE = (2000, 1, 1, 0, 0, 0)


def now_utc() -> str:
    """Current UTC time, ISO 8601 with a ``Z`` suffix, whole seconds."""
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z")


def new_id() -> str:
    """Unguessable identifier (uuid4, 122 random bits)."""
    return str(uuid.uuid4())


def valid_id(value: str) -> bool:
    """True if ``value`` is a canonical lower-case uuid4 string (the only ids ever issued).

    Every path under the data directory is built from an id that passed this check, so no
    user-supplied string can reach the file system.
    """
    if not isinstance(value, str) or len(value) != 36:
        return False
    try:
        u = uuid.UUID(value)
    except ValueError:
        return False
    return u.version == 4 and str(u) == value


def _write_json_atomic(path: Path, obj: Any) -> None:
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(json.dumps(obj, sort_keys=True, indent=1, allow_nan=False) + "\n",
                   encoding="utf-8")
    os.replace(tmp, path)


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


# --------------------------------------------------------------------------- store
_SCHEMA = """
CREATE TABLE IF NOT EXISTS cohorts (
    token TEXT PRIMARY KEY,
    created_utc TEXT NOT NULL,
    client_ip TEXT NOT NULL,
    info TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS jobs (
    id TEXT PRIMARY KEY,
    created_utc TEXT NOT NULL,
    updated_utc TEXT NOT NULL,
    started_utc TEXT,
    finished_utc TEXT,
    status TEXT NOT NULL CHECK (status IN ('queued', 'running', 'done', 'failed')),
    stage TEXT,
    progress REAL NOT NULL DEFAULT 0,
    message TEXT,
    error TEXT,
    client_ip TEXT NOT NULL,
    cohort_token TEXT NOT NULL,
    request TEXT NOT NULL,
    seq INTEGER NOT NULL
);
CREATE INDEX IF NOT EXISTS jobs_status ON jobs (status, seq);
"""


class JobStore:
    """SQLite-backed job and cohort table (thread-safe; one connection per call)."""

    def __init__(self, path: str | os.PathLike) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        with self._connect() as con:
            con.execute("PRAGMA journal_mode=WAL")
            con.executescript(_SCHEMA)

    def _connect(self) -> sqlite3.Connection:
        con = sqlite3.connect(self.path, timeout=30, isolation_level=None)
        con.row_factory = sqlite3.Row
        return con

    def _exec(self, sql: str, args: tuple = ()) -> list[sqlite3.Row]:
        with self._lock:
            con = self._connect()
            try:
                return con.execute(sql, args).fetchall()
            finally:
                con.close()

    # ---- cohorts
    def add_cohort(self, token: str, info: Mapping, client_ip: str) -> None:
        self._exec("INSERT INTO cohorts (token, created_utc, client_ip, info) VALUES (?, ?, ?, ?)",
                   (token, now_utc(), client_ip, json.dumps(info, sort_keys=True)))

    def get_cohort(self, token: str) -> dict | None:
        rows = self._exec("SELECT * FROM cohorts WHERE token = ?", (token,))
        if not rows:
            return None
        out = dict(rows[0])
        out["info"] = json.loads(out["info"])
        return out

    # ---- jobs
    def create_job(self, job_id: str, cohort_token: str, request: Mapping, client_ip: str,
                   created_utc: str) -> None:
        self._exec(
            "INSERT INTO jobs (id, created_utc, updated_utc, status, stage, progress, message, "
            "client_ip, cohort_token, request, seq) VALUES (?, ?, ?, 'queued', NULL, 0, ?, ?, ?, ?, "
            "(SELECT COALESCE(MAX(seq), 0) + 1 FROM jobs))",
            (job_id, created_utc, created_utc, "waiting in the queue", client_ip, cohort_token,
             json.dumps(request, sort_keys=True)))

    def get_job(self, job_id: str) -> dict | None:
        rows = self._exec("SELECT * FROM jobs WHERE id = ?", (job_id,))
        if not rows:
            return None
        out = dict(rows[0])
        out["request"] = json.loads(out["request"])
        if out["status"] == "queued":
            ahead = self._exec("SELECT COUNT(*) FROM jobs WHERE status = 'queued' AND seq < ?",
                               (out["seq"],))[0][0]
            out["queue_position"] = int(ahead) + 1
        else:
            out["queue_position"] = None
        return out

    def count(self, *statuses: str) -> int:
        marks = ",".join("?" * len(statuses))
        return int(self._exec(f"SELECT COUNT(*) FROM jobs WHERE status IN ({marks})",
                              tuple(statuses))[0][0])

    def active_for_ip(self, client_ip: str) -> int:
        return int(self._exec("SELECT COUNT(*) FROM jobs WHERE client_ip = ? AND status IN "
                              "('queued', 'running')", (client_ip,))[0][0])

    def claim_next(self) -> dict | None:
        """Mark the oldest queued job ``running`` and return it, unless a job is running.

        One UPDATE statement, so two workers (or two processes) can never both run a job.
        """
        t = now_utc()
        rows = self._exec(
            "UPDATE jobs SET status = 'running', started_utc = ?, updated_utc = ?, "
            "stage = 'starting', message = 'starting' "
            "WHERE id = (SELECT id FROM jobs WHERE status = 'queued' ORDER BY seq LIMIT 1) "
            "AND NOT EXISTS (SELECT 1 FROM jobs WHERE status = 'running') RETURNING id",
            (t, t))
        return self.get_job(rows[0]["id"]) if rows else None

    def update_progress(self, job_id: str, stage: str | None, progress: float,
                        message: str | None) -> None:
        self._exec("UPDATE jobs SET stage = ?, progress = ?, message = ?, updated_utc = ? "
                   "WHERE id = ? AND status = 'running'",
                   (stage, float(progress), message, now_utc(), job_id))

    def finish(self, job_id: str, status: str, *, error: str | None = None,
               message: str | None = None) -> None:
        if status not in ("done", "failed"):
            raise ValueError(f"finish status must be done or failed, got {status!r}")
        t = now_utc()
        self._exec("UPDATE jobs SET status = ?, error = ?, message = ?, finished_utc = ?, "
                   "updated_utc = ?, progress = CASE WHEN ? = 'done' THEN 1.0 ELSE progress END "
                   "WHERE id = ?", (status, error, message, t, t, status, job_id))

    def recover(self, reason: str) -> list[str]:
        """Fail every job left ``running`` by a previous process (queued jobs stay queued)."""
        rows = self._exec("SELECT id FROM jobs WHERE status = 'running'")
        ids = [r["id"] for r in rows]
        for job_id in ids:
            self.finish(job_id, "failed", error=reason, message="failed")
        return ids

    def expired(self, before_utc: str) -> list[str]:
        """Finished jobs whose ``finished_utc`` is before ``before_utc``."""
        rows = self._exec("SELECT id FROM jobs WHERE status IN ('done', 'failed') AND "
                          "finished_utc < ?", (before_utc,))
        return [r["id"] for r in rows]

    def delete_job(self, job_id: str) -> None:
        self._exec("DELETE FROM jobs WHERE id = ?", (job_id,))

    def expired_cohorts(self, before_utc: str) -> list[str]:
        rows = self._exec("SELECT token FROM cohorts WHERE created_utc < ? AND token NOT IN "
                          "(SELECT cohort_token FROM jobs)", (before_utc,))
        return [r["token"] for r in rows]

    def delete_cohort(self, token: str) -> None:
        self._exec("DELETE FROM cohorts WHERE token = ?", (token,))


# --------------------------------------------------------------------------- worker
class Worker(threading.Thread):
    """Single background worker: runs queued jobs one at a time in a child process.

    ``jobs_dir`` holds one directory per job; ``timeout_s`` is the hard wall-time limit of
    a job (the child is killed and the job failed with the reason).
    """

    def __init__(self, store: JobStore, jobs_dir: str | os.PathLike, *, timeout_s: float,
                 poll_s: float = 1.0, env: Mapping[str, str] | None = None) -> None:
        super().__init__(name="topoqspr-web-worker", daemon=True)
        self.store = store
        self.jobs_dir = Path(jobs_dir)
        self.timeout_s = float(timeout_s)
        self.poll_s = float(poll_s)
        self.extra_env = dict(env or {})
        self._stop_event = threading.Event()
        self.current: str | None = None

    def stop(self, join_s: float = 10.0) -> None:
        self._stop_event.set()
        if self.is_alive():
            self.join(join_s)

    def run(self) -> None:  # pragma: no cover - exercised through the API tests
        log = logging.getLogger("topoqspr.web")
        while not self._stop_event.is_set():
            try:
                job = self.store.claim_next()
            except Exception:  # noqa: BLE001 - keep the worker alive, log the cause
                log.exception("job claim failed")
                job = None
            if job is None:
                self._stop_event.wait(self.poll_s)
                continue
            self.current = job["id"]
            try:
                self._run(job)
            except Exception as exc:  # noqa: BLE001 - recorded on the job, never lost
                log.exception("job %s crashed in the worker", job["id"])
                self.store.finish(job["id"], "failed", error=f"internal worker error: {exc!r}"[:500],
                                  message="failed")
            finally:
                self.current = None

    def _child_env(self) -> dict[str, str]:
        env = dict(os.environ)
        pkg_root = str(Path(__file__).resolve().parents[2])  # .../src (the topoqspr parent)
        env["PYTHONPATH"] = os.pathsep.join([pkg_root] + ([env["PYTHONPATH"]]
                                                          if env.get("PYTHONPATH") else []))
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            env.setdefault(k, "1")
        env.setdefault("PYTHONHASHSEED", "0")
        env.update(self.extra_env)
        return env

    def _run(self, job: Mapping) -> None:
        job_id = job["id"]
        job_dir = self.jobs_dir / job_id
        if not valid_id(job_id) or not (job_dir / "config.json").exists():
            self.store.finish(job_id, "failed", error="job files missing on disk", message="failed")
            return
        t0 = time.monotonic()
        with open(job_dir / "worker.log", "ab") as logf:
            proc = subprocess.Popen([sys.executable, "-m", "topoqspr.web.jobs", str(job_dir)],
                                    cwd=str(job_dir), stdout=logf, stderr=subprocess.STDOUT,
                                    stdin=subprocess.DEVNULL, env=self._child_env())
            reason = None
            while True:
                try:
                    rc = proc.wait(timeout=self.poll_s)
                    break
                except subprocess.TimeoutExpired:
                    pass
                prog = _read_json(job_dir / "progress.json")
                if isinstance(prog, dict):
                    self.store.update_progress(job_id, prog.get("stage"),
                                               float(prog.get("progress", 0.0)),
                                               prog.get("message"))
                if time.monotonic() - t0 > self.timeout_s:
                    reason = (f"timed out: the job ran longer than the web limit of "
                              f"{self.timeout_s:.0f} s and was stopped")
                elif self._stop_event.is_set():
                    reason = "the server stopped while the job was running; please submit it again"
                if reason:
                    proc.kill()
                    proc.wait()
                    rc = None
                    break
        if reason is None and rc == 0 and (job_dir / "ladder.json").exists():
            self.store.finish(job_id, "done", message="finished")
            return
        if reason is None:
            err = _read_json(job_dir / "error.json")
            reason = (err.get("error") if isinstance(err, dict) and err.get("error")
                      else f"the job process exited with code {rc} (see the server log)")
        shutil.rmtree(job_dir / "results", ignore_errors=True)
        self.store.finish(job_id, "failed", error=str(reason)[:2000], message="failed")


# --------------------------------------------------------------------------- child
class _ProgressHandler(logging.Handler):
    """Turns the pipeline's log records into ``progress.json`` updates."""

    def __init__(self, path: Path, n_targets: int, n_learners: int) -> None:
        super().__init__(logging.INFO)
        self.path = path
        self.n_learners = max(int(n_learners), 1)
        self.n_eval = max(int(n_targets), 1) * self.n_learners
        self.done: set[str] = set()
        self.n_eval_done = 0
        self.write("cohort", "reading the cohort manifest")

    def progress(self) -> float:
        p = 0.0
        for name, w in STAGES:
            if name in self.done:
                p += w
            elif name == "evaluate":
                p += w * min(self.n_eval_done / self.n_eval, 1.0)
        return min(p, 0.99)

    def current(self) -> str:
        if self.n_eval_done >= self.n_eval:
            self.done.add("evaluate")
        for name, _ in STAGES:
            if name not in self.done:
                return name
        return "write"

    def write(self, stage: str, message: str) -> None:
        _write_json_atomic(self.path, {"stage": stage, "progress": round(self.progress(), 4),
                                       "message": message})

    def emit(self, record: logging.LogRecord) -> None:
        try:
            msg, args = record.msg, record.args or ()
            if msg.startswith("stage %s") and args:
                self.done.add(str(args[0]).split(":")[0])
            elif msg.startswith("evaluate %s / %s") and args:
                self.n_eval_done += 1
            elif msg.startswith("target %s not evaluated") and args:
                self.n_eval_done += self.n_learners
            else:
                return
            stage = self.current()
            text = {"standardize": "standardising structures",
                    "features": "computing 102 topological indices and baselines",
                    "targets": "computing / matching target values",
                    "evaluate": f"nested cross-validation ({self.n_eval_done} of {self.n_eval} "
                                "target x learner evaluations done)",
                    "protocol": "applicability-domain report",
                    "write": "writing result files"}.get(stage, stage)
            self.write(stage, text)
        except Exception:  # noqa: BLE001 - progress reporting must never break a run
            pass


def _num(x: Any) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _flag(x: Any) -> bool | None:
    if isinstance(x, bool):
        return x
    s = str(x).strip().lower()
    if s in ("true", "1"):
        return True
    if s in ("false", "0"):
        return False
    return None


def build_ladder(results_dir: str | os.PathLike, learner: str = "selected") -> dict:
    """The baseline-ladder table served by ``/api/jobs/{id}/results``.

    One entry per target (run-record order) with its spec, n, missing count, skip reason,
    and one row per feature set (``LADDER_SETS`` order): fold-mean and pooled R², MAE,
    paired ΔR² and ΔMAE against ``size`` with 95 % CI and p, valid/failed folds, the
    underpowered flag and reason, and the AD coverage (share of test rows outside the
    applicability domain, pooled).  Undefined numbers are ``null`` with the summary's
    ``reason``.  Nothing is interpreted.
    """
    from ..pipeline import read_results

    res = Path(results_dir)
    record = json.loads((res / "run_record.json").read_text(encoding="utf-8"))
    summary = read_results(res, "summary", dtype=str)
    ad = read_results(res, "ad_summary", dtype=str) if (res / "ad_summary.csv").exists() else None
    skipped = (record.get("counts") or {}).get("skipped_targets") or {}
    n_values = (record.get("counts") or {}).get("targets") or {}
    out_targets = []
    for spec in record.get("targets") or []:
        name = spec["name"]
        sub = summary[(summary["target"] == name) & (summary["learner"] == learner)]
        rows = []
        for set_name in LADDER_SETS:
            r = sub[sub["set"] == set_name]
            if r.empty:
                continue
            r = r.iloc[0]
            cov = None
            ad_reason = None
            if ad is not None:
                a = ad[(ad["target"] == name) & (ad["learner"] == learner) & (ad["set"] == set_name)]
                if not a.empty:
                    cov = _num(a.iloc[0].get("share_outside_pooled"))
                    ad_reason = a.iloc[0].get("reason") or None
            is_ref = set_name == r.get("reference")
            rows.append(dict(
                set=set_name, n=_num(r.get("n")), p=_num(r.get("p")),
                r2_fold_mean=_num(r.get("r2_fold_mean")),
                r2_pooled_mean=_num(r.get("r2_pooled_mean")),
                mae_fold_mean=_num(r.get("mae_fold_mean")),
                d_r2_mean=None if is_ref else _num(r.get("d_r2_mean")),
                d_r2_ci_low=None if is_ref else _num(r.get("d_r2_ci_low")),
                d_r2_ci_high=None if is_ref else _num(r.get("d_r2_ci_high")),
                d_r2_p=None if is_ref else _num(r.get("d_r2_p")),
                d_mae_mean=None if is_ref else _num(r.get("d_mae_mean")),
                d_mae_ci_low=None if is_ref else _num(r.get("d_mae_ci_low")),
                d_mae_ci_high=None if is_ref else _num(r.get("d_mae_ci_high")),
                n_valid_folds=_num(r.get("n_valid_folds")),
                n_failed_folds=_num(r.get("n_failed_folds")),
                underpowered=_flag(r.get("underpowered")),
                underpowered_reason=r.get("underpowered_reason") or None,
                reason=r.get("reason") or None,
                is_reference=bool(is_ref),
                ad_share_outside=cov, ad_reason=ad_reason,
            ))
        out_targets.append(dict(
            name=name, target_class=spec.get("target_class"),
            target_class_label=spec.get("target_class_label"), units=spec.get("units"),
            source=spec.get("source"), description=spec.get("description"),
            n_values=n_values.get(name), skipped_reason=skipped.get(name),
            reference="size", learner=learner, rows=rows))
    ev = record["config"]["evaluation"]
    return dict(
        targets=out_targets,
        protocol=dict(mode="exploratory", learner=learner, families=ev["families"],
                      outer_splits=ev["outer_splits"], repeats=len(ev["seeds"]),
                      inner_splits=ev["inner_splits"], seeds=ev["seeds"], split=ev["split"],
                      ad=ev["ad"], min_n=ev["min_n"], min_n_test=ev["min_n_test"],
                      ci="paired corrected resampled t interval over outer folds "
                         "(EVALUATION_PROTOCOL.md); with one repeat it is wide"),
        chembl=record.get("chembl"),
        code=record.get("code"),
        topoqspr_version=record.get("topoqspr_version"),
        config_sha256=record.get("config_sha256"),
        manifest_sha256=record.get("manifest_sha256"),
        counts={k: (record.get("counts") or {}).get(k)
                for k in ("candidates", "standardised", "rejected", "connectivity_groups")},
        timings_s=record.get("timings_s"),
    )


def licence_note(record: Mapping) -> str:
    """``LICENCE-DATA.txt`` of a results archive: ChEMBL and every measured-target source."""
    from .. import expdata
    from ..chembl import CHEMBL_DATA_LICENCE

    lines = [
        "TopoQSPR web results - data licences",
        "=====================================",
        "",
        "Code: TopoQSPR, MIT licence (https://github.com/Singati2/TopoQSPR).",
        "",
        "ChEMBL data",
        "-----------",
    ]
    chembl = record.get("chembl") or {}
    lines += [
        f"manifest.csv, molecules.csv and every file derived from them contain data from "
        f"ChEMBL release {chembl.get('release', '?')} ({chembl.get('release_date', '?')}).",
        f"Licence: {CHEMBL_DATA_LICENCE}.",
        "Redistributed copies must keep this attribution and the same licence.",
        "Cite: Zdrazil B et al. The ChEMBL Database in 2023. Nucleic Acids Res. 52 (2024) "
        "D1180-D1192, doi:10.1093/nar/gkad1004.",
        "",
    ]
    measured = [t for t in record.get("targets") or [] if str(t.get("source", "")).startswith("expdata:")]
    if measured:
        lines += ["Measured targets", "----------------"]
        for t in measured:
            spec = expdata.TARGETS.get(t["name"])
            src = expdata.SOURCES.get(spec.source) if spec else None
            if src is None:
                lines.append(f"{t['name']}: source not recorded")
                continue
            lines += [f"{t['name']}: {src.title}", f"  citation: {src.citation}",
                      f"  doi: {src.doi}", f"  version: {src.version}",
                      f"  licence: {src.licence}", ""]
    lines += ["Computed targets and descriptors (RDKit, BSD licence) carry no third-party data "
              "licence beyond that of the input structures.", ""]
    return "\n".join(lines)


_README = """TopoQSPR web results (exploratory mode)
=======================================

This archive was produced by the TopoQSPR web interface in EXPLORATORY mode:
learner "selected" over ridge / PLS / kNN, {outer}-fold outer cross-validation with
{repeats} repeat(s), {inner}-fold inner loop, no Y-randomisation and no bootstrap.
Confirmatory analyses use the command-line protocol (topoqspr run with the frozen
study configuration; docs/PIPELINE.md and docs/EVALUATION_PROTOCOL.md).

Files: run_record.json (provenance: code commit, package versions, ChEMBL release,
config hash, file hashes, timings), manifest.csv (the cohort), molecules.csv,
features.csv, targets.csv, failures.csv (every undefined value with its reason),
folds.csv, predictions.csv, paired.csv, inner.csv, summary.csv (the baseline ladder),
ad.csv, ad_summary.csv, study.json, LICENCE-DATA.txt.

Read the CSVs with only empty cells as missing values (pandas: keep_default_na=False,
na_values=[""]); otherwise the feature set named "null" is lost.
"""


def _write_zip(job_dir: Path, results: Path, record: Mapping) -> Path:
    ev = record["config"]["evaluation"]
    members = sorted(p for p in results.iterdir() if p.is_file())
    out = job_dir / "results.zip"
    tmp = job_dir / ".results.zip.tmp"
    with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        def add(name: str, data: bytes) -> None:
            info = zipfile.ZipInfo(f"topoqspr-results/{name}", date_time=_ZIP_DATE)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, data)
        add("README.txt", _README.format(outer=ev["outer_splits"], repeats=len(ev["seeds"]),
                                         inner=ev["inner_splits"]).encode())
        add("LICENCE-DATA.txt", licence_note(record).encode())
        for p in members:
            add(p.name, p.read_bytes())
    os.replace(tmp, out)
    return out


def run_job(job_dir: str | os.PathLike) -> int:
    """Child-process entry point: run the job in ``job_dir``; 0 on success."""
    from ..pipeline import run_study

    job_dir = Path(job_dir)
    try:
        cfg = json.loads((job_dir / "config.json").read_text(encoding="utf-8"))
        meta = json.loads((job_dir / "request.json").read_text(encoding="utf-8"))
        n_targets = sum(len(t["names"]) for t in cfg["targets"])
        handler = _ProgressHandler(job_dir / "progress.json", n_targets,
                                   len(cfg["evaluation"]["learners"]))
        plog = logging.getLogger("topoqspr.pipeline")
        plog.setLevel(logging.INFO)
        plog.addHandler(handler)
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
        results = job_dir / "results"
        if results.exists():
            shutil.rmtree(results)
        run_study(cfg, results, timestamp_utc=meta["created_utc"], base_dir=job_dir)
        handler.done.update(name for name, _ in STAGES)
        record = json.loads((results / "run_record.json").read_text(encoding="utf-8"))
        handler.write("write", "packaging results")
        _write_zip(job_dir, results, record)
        _write_json_atomic(job_dir / "ladder.json", build_ladder(results))
        _write_json_atomic(job_dir / "progress.json",
                           {"stage": "done", "progress": 1.0, "message": "finished"})
        return 0
    except Exception as exc:  # noqa: BLE001 - reported to the user through error.json
        traceback.print_exc()
        msg = str(exc) or type(exc).__name__
        _write_json_atomic(job_dir / "error.json",
                           {"error": f"{type(exc).__name__}: {msg}"[:2000]})
        return 1


if __name__ == "__main__":  # pragma: no cover - run by the worker
    if len(sys.argv) != 2:
        print("usage: python -m topoqspr.web.jobs JOB_DIR", file=sys.stderr)
        sys.exit(2)
    sys.exit(run_job(sys.argv[1]))
