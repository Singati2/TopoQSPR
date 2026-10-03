"""Run records: everything needed to identify, audit and re-run a study.

A run record (``run_record.json``, written last into a results directory by
:func:`topoqspr.pipeline.run_study`) holds

* ``code``: the git commit of the source tree the package was imported from,
  whether tracked files differ from that commit (``dirty``) and the number of
  untracked files, plus ``topoqspr.__version__``. An installed copy outside
  any git tree (the Docker image) falls back to the commit written into
  ``topoqspr/_build_info.py`` at build time (:func:`build_info`); ``source``
  says which (``"git"``, ``"build_info"`` or ``None``);
* ``chembl`` also names the data licence (``data_licence``) of ChEMBL-derived
  results;
* ``environment``: Python version and implementation, platform, the version of
  every installed distribution (``importlib.metadata``) and the environment
  variables that can change numerical results (``PYTHONHASHSEED`` and the
  BLAS / OpenMP thread counts);
* ``chembl``: the release and its date as read from the client (``None`` when
  the cohort did not come from ChEMBL);
* ``config``: the full resolved configuration and ``config_sha256``, the
  SHA-256 of the canonical JSON of its *study* part (see :func:`config_sha256`);
* ``seeds``, ``inputs`` (SHA-256 of every input file), ``manifest_sha256``,
  ``outputs`` (SHA-256 and size of every file written), ``timings_s`` (wall
  time per stage) and ``created_utc``, which is **passed in by the caller**
  so that nothing in the record depends on reading the clock here.

Canonical JSON (:func:`canonical_json`) is ``json.dumps`` with sorted keys,
``(",", ":")`` separators, ``ensure_ascii=True`` and ``allow_nan=False``:
the same object always gives the same bytes, independent of key order and
of ``PYTHONHASHSEED``.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from importlib import metadata
from pathlib import Path
from typing import Any, Mapping

from . import __version__

__all__ = [
    "RUN_RECORD_SCHEMA",
    "RESULTS_ENV_VARS",
    "canonical_json",
    "sha256_bytes",
    "sha256_file",
    "config_sha256",
    "git_info",
    "build_info",
    "package_versions",
    "environment_info",
    "file_records",
    "build_run_record",
    "write_run_record",
    "read_run_record",
]

#: Version of the run-record layout; bump on incompatible changes.
RUN_RECORD_SCHEMA = 1

#: Environment variables that can change numerical results or their order.
RESULTS_ENV_VARS = (
    "PYTHONHASHSEED",
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OMP_PROC_BIND",
)

#: Top-level config keys that do not affect results and are left out of the hash.
RUNTIME_KEY = "runtime"


def canonical_json(obj: Any) -> str:
    """Deterministic JSON text: sorted keys, compact separators, ASCII, no NaN."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                      allow_nan=False)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | os.PathLike) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def config_sha256(config: Mapping) -> str:
    """SHA-256 of the canonical JSON of ``config`` without its ``runtime`` part.

    ``runtime`` holds settings that cannot change any result (``n_jobs``,
    cache locations, offline mode), so two runs of the same study on
    different machines share the hash. Pass a *resolved* config
    (:func:`topoqspr.pipeline.resolve_config`) so that defaults written out
    and defaults left implicit hash the same.
    """
    study = {k: v for k, v in config.items() if k != RUNTIME_KEY}
    return sha256_bytes(canonical_json(study).encode("ascii"))


def _git(args: list[str], cwd: Path) -> str | None:
    try:
        out = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True,
                             timeout=30, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout if out.returncode == 0 else None


def git_info(path: str | os.PathLike | None = None) -> dict:
    """Commit, dirty flag and untracked-file count of the repository at ``path``.

    ``path`` defaults to the directory this package was imported from, so the
    record names the code that actually ran. Outside a git checkout (an
    installed wheel) the commit baked in at build time is used when there is
    one (:func:`build_info`); otherwise every field is ``None`` and ``reason``
    says why.
    """
    here = Path(path) if path is not None else Path(__file__).resolve().parent
    top = _git(["rev-parse", "--show-toplevel"], here)
    if top is None:
        if path is None:
            built = build_info()
            if built is not None:
                return built
        return dict(commit=None, dirty=None, untracked=None, describe=None, root=None,
                    source=None,
                    reason="not inside a git work tree (or git unavailable) and no "
                           "topoqspr/_build_info.py (see the Dockerfile)")
    commit = (_git(["rev-parse", "HEAD"], here) or "").strip() or None
    tracked = _git(["status", "--porcelain", "--untracked-files=no"], here)
    untracked = _git(["ls-files", "--others", "--exclude-standard"], here)
    describe = _git(["describe", "--tags", "--always", "--dirty"], here)
    return dict(
        commit=commit,
        dirty=None if tracked is None else bool(tracked.strip()),
        untracked=None if untracked is None else len([x for x in untracked.splitlines() if x]),
        describe=None if describe is None else describe.strip(),
        root=top.strip(),
        source="git",
        reason="",
    )


def build_info() -> dict | None:
    """The commit baked into an installed copy at build time, or ``None``.

    ``docker build --build-arg GIT_COMMIT=... --build-arg GIT_DIRTY=...``
    writes ``topoqspr/_build_info.py`` before the package is installed, so a
    container (where the package lives in site-packages, outside any git
    tree) still names the code it runs. Same keys as :func:`git_info`, with
    ``source="build_info"``.
    """
    try:
        from . import _build_info  # type: ignore[attr-defined]
    except ImportError:
        return None
    commit = str(getattr(_build_info, "GIT_COMMIT", "") or "").strip().lower()
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        return None
    dirty = {"1": True, "true": True, "0": False, "false": False}.get(
        str(getattr(_build_info, "GIT_DIRTY", "")).strip().lower())
    return dict(commit=commit, dirty=dirty, untracked=None,
                describe=str(getattr(_build_info, "GIT_DESCRIBE", "") or "") or None,
                root=None, source="build_info",
                reason="installed copy; commit recorded at build time (topoqspr/_build_info.py)")


def package_versions() -> dict[str, str]:
    """``{distribution name (lower case): version}`` for every installed distribution."""
    out: dict[str, str] = {}
    for dist in metadata.distributions():
        name = dist.metadata.get("Name")
        if name:
            out.setdefault(name.lower(), dist.version)
    return dict(sorted(out.items()))


def environment_info() -> dict:
    return dict(
        python=sys.version,
        python_version=platform.python_version(),
        implementation=platform.python_implementation(),
        executable=sys.executable,
        platform=platform.platform(),
        machine=platform.machine(),
        system=platform.system(),
        env={k: os.environ.get(k) for k in RESULTS_ENV_VARS},
        packages=package_versions(),
    )


def file_records(root: str | os.PathLike, names: list[str]) -> dict[str, dict]:
    """``{name: {"sha256", "bytes"}}`` for files under ``root`` (sorted by name)."""
    root = Path(root)
    return {n: dict(sha256=sha256_file(root / n), bytes=(root / n).stat().st_size)
            for n in sorted(names)}


def build_run_record(*, config: Mapping, timestamp_utc: str, chembl: Mapping | None,
                     seeds: Mapping, inputs: Mapping, manifest_sha256: str | None,
                     outputs: Mapping, timings_s: Mapping, counts: Mapping,
                     base_dir: str | None, code: Mapping | None = None,
                     extra: Mapping | None = None) -> dict:
    """Assemble the run record (pure function of its arguments + the environment).

    ``code`` is :func:`git_info` taken when the run *started* (the code that was
    imported); it is read now if not given.
    """
    if not isinstance(timestamp_utc, str) or not timestamp_utc:
        raise ValueError("timestamp_utc must be supplied by the caller (ISO 8601 string)")
    return dict(
        schema=RUN_RECORD_SCHEMA,
        created_utc=timestamp_utc,
        topoqspr_version=__version__,
        code=dict(code) if code is not None else git_info(),
        environment=environment_info(),
        chembl=dict(chembl) if chembl is not None else None,
        config=dict(config),
        config_sha256=config_sha256(config),
        base_dir=base_dir,
        seeds=dict(seeds),
        inputs=dict(inputs),
        manifest_sha256=manifest_sha256,
        counts=dict(counts),
        outputs=dict(outputs),
        timings_s={k: round(float(v), 3) for k, v in timings_s.items()},
        **(dict(extra) if extra else {}),
    )


def write_run_record(record: Mapping, path: str | os.PathLike) -> Path:
    """Write the record as indented, key-sorted JSON (UTF-8, trailing newline)."""
    path = Path(path)
    path.write_text(json.dumps(record, sort_keys=True, indent=1, allow_nan=False) + "\n",
                    encoding="utf-8")
    return path


def read_run_record(path: str | os.PathLike) -> dict:
    """Load a run record and check its schema and config hash."""
    rec = json.loads(Path(path).read_text(encoding="utf-8"))
    if rec.get("schema") != RUN_RECORD_SCHEMA:
        raise ValueError(f"unsupported run-record schema {rec.get('schema')!r} "
                         f"(this version reads {RUN_RECORD_SCHEMA})")
    want = rec.get("config_sha256")
    got = config_sha256(rec["config"])
    if want != got:
        raise ValueError(f"run record config hash mismatch: recorded {want}, recomputed {got}")
    return rec
