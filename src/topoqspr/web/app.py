"""TopoQSPR v1.0 web interface (FastAPI; exploratory mode).  See docs/WEB.md.

Flow (mirrors the paper's method):

1. ``GET  /api/headings?q=`` - MeSH heading suggestions from the release's
   ``drug_indication`` heading index (cohort.suggest_headings; cached per release).
2. ``POST /api/cohort`` - builds the controlled cohort for the confirmed MeSH IDs and
   returns the manifest preview (included and excluded parents with reasons, counts,
   ChEMBL release).  No model is run.
3. ``POST /api/run`` - after the user has confirmed the manifest: queue an exploratory
   job (learner ``selected`` over ridge / PLS / kNN, 5-fold outer CV x 1 repeat, no
   permutation or bootstrap).
4. ``GET  /api/jobs/{id}`` - status and progress; ``GET /api/jobs/{id}/results`` - the
   baseline ladder (null, size, composition, TI, rdkit2d) per target.
5. Downloads: ``/api/cohorts/{token}/manifest.csv``, ``/api/jobs/{id}/manifest.csv``,
   ``/api/jobs/{id}/results.zip``, ``/api/jobs/{id}/run_record.json``, and
   ``POST /api/indices`` (indices of a pasted / uploaded SMILES list, as a zip).

``app`` is built from environment variables (:meth:`Settings.from_env`); tests and
embedders call :func:`create_app` with explicit :class:`Settings`.
"""

from __future__ import annotations

import asyncio
import csv
import io
import json
import logging
import math
import os
import re
import shutil
import threading
import time
import zipfile
from collections import deque
from contextlib import asynccontextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, ConfigDict, Field, field_validator

from .. import __version__
from . import jobs as jb

__all__ = ["Settings", "create_app", "app"]

log = logging.getLogger("topoqspr.web")

STATIC_DIR = Path(__file__).resolve().parent / "static"

#: Computed targets offered on the web (class 1 MolWt is a sanity check only).
RDKIT_TARGETS = ("MolLogP", "TPSA", "MolMR", "MolWt")
#: Mesh descriptor IDs: ``D`` + 6-9 digits.
_MESH_RE = re.compile(r"^D\d{6,9}$")
_CONTROL = re.compile(r"[\x00-\x1f\x7f]")


def _env_bool(name: str, default: bool) -> bool:
    v = os.environ.get(name)
    if v is None or v == "":
        return default
    return v.strip().lower() in ("1", "true", "yes", "on")


def _env_list(name: str) -> tuple[str, ...]:
    v = os.environ.get(name, "")
    return tuple(x.strip() for x in v.split(",") if x.strip())


@dataclass(frozen=True)
class Settings:
    """Server settings (environment variables ``TOPOQSPR_WEB_*``; docs/WEB.md §3)."""

    data_dir: Path = Path(os.path.expanduser("~/.cache/topoqspr/web"))
    chembl_cache_dir: Optional[str] = None      # None = topoqspr default (~/.cache/topoqspr/chembl)
    expdata_cache_dir: Optional[str] = None
    offline: bool = False                       # serve ChEMBL from the cache only
    chembl_release: Optional[str] = None        # pin, e.g. "ChEMBL_37"
    allowed_origins: tuple[str, ...] = ()       # CORS; empty = same origin only
    mesh_allowlist: tuple[str, ...] = ()        # empty = every heading; else only these IDs
    max_molecules: int = 300                    # included parents per cohort that may run
    min_molecules: int = 10
    max_indices_molecules: int = 300            # /api/indices
    max_smiles_length: int = 1000
    max_body_bytes: int = 256 * 1024
    max_queue: int = 5                          # queued + running jobs, all users
    max_active_per_ip: int = 1
    job_timeout_s: float = 1800.0
    rate_api_per_min: int = 120                 # every /api request
    rate_cohort_per_min: int = 10
    rate_indices_per_min: int = 10
    rate_runs_per_hour: int = 6
    trust_proxy: bool = False                   # take the client IP from X-Forwarded-For
    retention_days: float = 14.0                # finished jobs / unused cohorts kept this long
    n_jobs: int = 2                             # outer-fold parallelism inside a job
    worker: bool = True                         # start the background worker
    # evaluation of the exploratory mode (overridable for tests only)
    families: tuple[str, ...] = ("ridge", "pls", "knn")
    outer_splits: int = 5
    inner_splits: int = 3
    seeds: tuple[int, ...] = (0,)

    @classmethod
    def from_env(cls) -> "Settings":
        e = os.environ.get
        d = cls()
        return cls(
            data_dir=Path(e("TOPOQSPR_WEB_DATA_DIR") or d.data_dir),
            chembl_cache_dir=e("TOPOQSPR_WEB_CHEMBL_CACHE_DIR") or None,
            expdata_cache_dir=e("TOPOQSPR_WEB_EXPDATA_CACHE_DIR") or None,
            offline=_env_bool("TOPOQSPR_WEB_OFFLINE", False),
            chembl_release=e("TOPOQSPR_WEB_CHEMBL_RELEASE") or None,
            allowed_origins=_env_list("TOPOQSPR_WEB_ORIGINS"),
            mesh_allowlist=tuple(m.upper() for m in _env_list("TOPOQSPR_WEB_MESH_ALLOWLIST")),
            max_molecules=int(e("TOPOQSPR_WEB_MAX_MOLECULES") or d.max_molecules),
            max_indices_molecules=int(e("TOPOQSPR_WEB_MAX_INDICES") or d.max_indices_molecules),
            max_body_bytes=int(e("TOPOQSPR_WEB_MAX_BODY_BYTES") or d.max_body_bytes),
            max_queue=int(e("TOPOQSPR_WEB_MAX_QUEUE") or d.max_queue),
            job_timeout_s=float(e("TOPOQSPR_WEB_JOB_TIMEOUT_S") or d.job_timeout_s),
            rate_api_per_min=int(e("TOPOQSPR_WEB_RATE_API_PER_MIN") or d.rate_api_per_min),
            rate_runs_per_hour=int(e("TOPOQSPR_WEB_RATE_RUNS_PER_HOUR") or d.rate_runs_per_hour),
            trust_proxy=_env_bool("TOPOQSPR_WEB_TRUST_PROXY", False),
            retention_days=float(e("TOPOQSPR_WEB_RETENTION_DAYS") or d.retention_days),
            n_jobs=max(1, min(4, int(e("TOPOQSPR_WEB_N_JOBS") or d.n_jobs))),
            worker=_env_bool("TOPOQSPR_WEB_WORKER", True),
        )

    @property
    def jobs_dir(self) -> Path:
        return Path(self.data_dir) / "jobs"

    @property
    def cohorts_dir(self) -> Path:
        return Path(self.data_dir) / "cohorts"


# --------------------------------------------------------------------------- limits
class RateLimiter:
    """Sliding-window request counter per (bucket, client IP), in memory."""

    def __init__(self) -> None:
        self._hits: dict[tuple[str, str], deque] = {}
        self._lock = threading.Lock()

    def check(self, bucket: str, ip: str, limit: int, window_s: float) -> float | None:
        """Record a hit; return ``None`` if allowed, else seconds until the next slot."""
        if limit <= 0:
            return None
        now = time.monotonic()
        with self._lock:
            q = self._hits.setdefault((bucket, ip), deque())
            while q and now - q[0] >= window_s:
                q.popleft()
            if len(q) >= limit:
                return max(window_s - (now - q[0]), 1.0)
            q.append(now)
            if len(self._hits) > 50_000:  # bound memory under abuse
                for k in [k for k, v in self._hits.items() if not v][:25_000]:
                    del self._hits[k]
            return None


class _BodyLimit:
    """ASGI middleware: reject request bodies above ``max_bytes`` with 413.

    Checks ``Content-Length`` up front and counts the bytes of chunked bodies as they arrive.
    """

    def __init__(self, app, max_bytes: int) -> None:
        self.app = app
        self.max_bytes = int(max_bytes)

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            return await self.app(scope, receive, send)
        headers = dict(scope.get("headers") or [])
        cl = headers.get(b"content-length")
        too_big = JSONResponse({"detail": f"request body larger than {self.max_bytes} bytes"},
                               status_code=413)
        if cl is not None:
            try:
                if int(cl) > self.max_bytes:
                    return await too_big(scope, receive, send)
            except ValueError:
                return await JSONResponse({"detail": "invalid Content-Length"},
                                          status_code=400)(scope, receive, send)
        seen = 0

        class _TooLarge(Exception):
            pass

        async def limited_receive():
            nonlocal seen
            msg = await receive()
            if msg["type"] == "http.request":
                seen += len(msg.get("body", b""))
                if seen > self.max_bytes:
                    raise _TooLarge
            return msg

        started = False

        async def tracking_send(msg):
            nonlocal started
            if msg["type"] == "http.response.start":
                started = True
            await send(msg)

        try:
            await self.app(scope, limited_receive, tracking_send)
        except _TooLarge:
            if not started:
                await too_big(scope, receive, send)


_CSP = ("default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; "
        "connect-src 'self'; font-src 'self'; object-src 'none'; base-uri 'none'; "
        "form-action 'self'; frame-ancestors 'none'")


# --------------------------------------------------------------------------- models
_PHASES = (0.5, 1.0, 2.0, 3.0, 4.0)
_REF_TYPES = Literal["ATC", "ClinicalTrials", "DailyMed", "EMA", "FDA", "INN", "USAN"]


class CohortRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mesh_ids: list[str] = Field(min_length=1, max_length=10)
    min_phase: Optional[float] = 4.0
    evidence_sources: Optional[list[_REF_TYPES]] = Field(default=None, max_length=7)
    exclude_sole_evidence: Optional[list[_REF_TYPES]] = Field(default=None, max_length=7)
    exclude_diagnostics: bool = False
    disease_label: Optional[str] = Field(default=None, max_length=200)

    @field_validator("mesh_ids")
    @classmethod
    def _mesh(cls, v: list[str]) -> list[str]:
        out = sorted({str(m).strip().upper() for m in v})
        bad = [m for m in out if not _MESH_RE.match(m)]
        if bad:
            raise ValueError(f"not MeSH descriptor IDs (D followed by 6-9 digits): {bad[:5]}")
        return out

    @field_validator("min_phase")
    @classmethod
    def _phase(cls, v):
        if v is not None and v not in _PHASES:
            raise ValueError(f"min_phase must be one of {list(_PHASES)} or null (no phase rule)")
        return v

    @field_validator("evidence_sources", "exclude_sole_evidence")
    @classmethod
    def _nonempty(cls, v):
        if v is not None and not v:
            return None
        return sorted(set(v)) if v else v

    @field_validator("disease_label")
    @classmethod
    def _label(cls, v):
        if v is None:
            return None
        v = v.strip()
        if _CONTROL.search(v):
            raise ValueError("disease_label contains control characters")
        return v or None


class Targets(BaseModel):
    model_config = ConfigDict(extra="forbid")

    computed: list[Literal["MolLogP", "TPSA", "MolMR", "MolWt"]] = Field(default_factory=list,
                                                                       max_length=4)
    measured: list[str] = Field(default_factory=list, max_length=8)

    @field_validator("measured")
    @classmethod
    def _measured(cls, v: list[str]) -> list[str]:
        from .. import expdata

        bad = [x for x in v if x not in expdata.TARGETS]
        if bad:
            raise ValueError(f"unknown measured targets {bad[:5]}; available: "
                             f"{sorted(expdata.TARGETS)}")
        return v


class RunRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    cohort_token: str = Field(min_length=36, max_length=36)
    targets: Targets
    mode: Literal["exploratory"] = "exploratory"
    confirmed: bool = Field(description="the user has reviewed the cohort manifest")


class IndicesRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    text: str = Field(min_length=1, max_length=200_000)


# --------------------------------------------------------------------------- helpers
def _client_ip(request: Request, settings: Settings) -> str:
    if settings.trust_proxy:
        xff = request.headers.get("x-forwarded-for", "")
        parts = [p.strip() for p in xff.split(",") if p.strip()]
        if parts:
            return parts[-1][:64]  # the address the trusted proxy saw
    return (request.client.host if request.client else "unknown")[:64]


def _clean(x: Any) -> Any:
    """JSON-safe value: NaN/inf -> None, numpy scalars -> Python."""
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if hasattr(x, "item") and not isinstance(x, (str, bytes)):
        try:
            return _clean(x.item())
        except (ValueError, AttributeError):
            return str(x)
    return x


def _dict_rows(df, cols) -> list[dict]:
    return [{c: _clean(v) for c, v in zip(cols, rec)} for rec in df[list(cols)].itertuples(index=False)]


def parse_smiles_text(text: str, max_rows: int, max_len: int) -> list[tuple[str, str]]:
    """Parse a pasted / uploaded list into ``(id, smiles)`` pairs.

    Accepted: a CSV whose header has a ``smiles`` column (and optionally ``id``), or one
    molecule per line as ``SMILES [name]`` (the .smi convention).  Blank lines and lines
    starting with ``#`` are skipped.  Raises ``ValueError`` with a clear message.
    """
    lines = [ln for ln in text.replace("\r\n", "\n").replace("\r", "\n").split("\n")]
    body = [ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
    if not body:
        raise ValueError("no molecules found")
    header = [h.strip().lower() for h in next(csv.reader([body[0]]))]
    pairs: list[tuple[str, str]] = []
    if "smiles" in header and ("," in body[0] or len(header) == 1):
        si = header.index("smiles")
        ii = header.index("id") if "id" in header else None
        for k, row in enumerate(csv.reader(body[1:]), start=1):
            if not row:
                continue
            if si >= len(row):
                raise ValueError(f"row {k}: no smiles column")
            rid = row[ii].strip() if ii is not None and ii < len(row) else f"mol{k}"
            pairs.append((rid, row[si].strip()))
    else:
        for k, ln in enumerate(body, start=1):
            parts = ln.strip().split(None, 1)
            pairs.append(((parts[1].strip() if len(parts) > 1 else f"mol{k}"), parts[0]))
    if not pairs:
        raise ValueError("no molecules found")
    if len(pairs) > max_rows:
        raise ValueError(f"{len(pairs)} molecules; the web limit is {max_rows} "
                         "(use the command line `topoqspr indices` for larger sets)")
    ids = [p[0] for p in pairs]
    if len(set(ids)) != len(ids):
        dup = sorted({i for i in ids if ids.count(i) > 1})[:5]
        raise ValueError(f"duplicate ids: {dup}")
    for rid, smi in pairs:
        if not smi:
            raise ValueError(f"{rid}: empty SMILES")
        if len(smi) > max_len:
            raise ValueError(f"{rid}: SMILES longer than {max_len} characters")
        if len(rid) > 100 or _CONTROL.search(rid):
            raise ValueError(f"id {rid[:20]!r}...: ids must be at most 100 printable characters")
    return pairs


def _job_view(job: dict) -> dict:
    req = job["request"]
    return dict(job_id=job["id"], status=job["status"], stage=job["stage"],
                progress=round(float(job["progress"] or 0.0), 4), message=job["message"],
                error=job["error"], created_utc=job["created_utc"],
                started_utc=job["started_utc"], finished_utc=job["finished_utc"],
                queue_position=job.get("queue_position"), mode=req.get("mode"),
                targets=req.get("targets"), cohort=req.get("cohort"))


# --------------------------------------------------------------------------- app
def create_app(settings: Settings | None = None) -> FastAPI:
    """Build the ASGI app.  The worker starts and crashed jobs are recovered at startup."""
    settings = settings or Settings.from_env()
    state: dict[str, Any] = {"headings": {}, "worker": None}
    limiter = RateLimiter()
    cohort_sem = threading.BoundedSemaphore(2)
    indices_sem = threading.BoundedSemaphore(1)

    def store() -> jb.JobStore:
        return state["store"]

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        Path(settings.data_dir).mkdir(parents=True, exist_ok=True)
        settings.jobs_dir.mkdir(parents=True, exist_ok=True)
        settings.cohorts_dir.mkdir(parents=True, exist_ok=True)
        st = jb.JobStore(Path(settings.data_dir) / "web.sqlite")
        state["store"] = st
        failed = st.recover("the server restarted while this job was running; "
                            "please submit it again")
        for job_id in failed:
            shutil.rmtree(settings.jobs_dir / job_id / "results", ignore_errors=True)
        if failed:
            log.warning("marked %d interrupted job(s) failed: %s", len(failed), failed)
        _purge(st)
        if settings.worker:
            w = jb.Worker(st, settings.jobs_dir, timeout_s=settings.job_timeout_s)
            w.start()
            state["worker"] = w
        try:
            yield
        finally:
            w = state.get("worker")
            if w is not None:
                await asyncio.to_thread(w.stop)

    def _purge(st: jb.JobStore) -> None:
        if settings.retention_days <= 0:
            return
        import datetime as dt
        cutoff = (dt.datetime.now(dt.timezone.utc)
                  - dt.timedelta(days=settings.retention_days)).replace(microsecond=0)
        before = cutoff.isoformat().replace("+00:00", "Z")
        for job_id in st.expired(before):
            if jb.valid_id(job_id):
                shutil.rmtree(settings.jobs_dir / job_id, ignore_errors=True)
            st.delete_job(job_id)
        for token in st.expired_cohorts(before):
            if jb.valid_id(token):
                shutil.rmtree(settings.cohorts_dir / token, ignore_errors=True)
            st.delete_cohort(token)

    app = FastAPI(title="TopoQSPR web (v1.0, exploratory)", version=__version__,
                  lifespan=lifespan, docs_url=None, redoc_url=None, openapi_url="/api/openapi.json")
    app.state.settings = settings
    app.state.web = state

    if settings.allowed_origins:
        app.add_middleware(CORSMiddleware, allow_origins=list(settings.allowed_origins),
                           allow_methods=["GET", "POST"], allow_headers=["Content-Type"],
                           allow_credentials=False, max_age=600)

    @app.middleware("http")
    async def guard(request: Request, call_next):
        path = request.url.path
        if path.startswith("/api/"):
            ip = _client_ip(request, settings)
            wait = limiter.check("api", ip, settings.rate_api_per_min, 60.0)
            if wait is not None:
                return JSONResponse({"detail": "too many requests; slow down"}, status_code=429,
                                    headers={"Retry-After": str(int(math.ceil(wait)))})
        response = await call_next(request)
        response.headers.setdefault("Content-Security-Policy", _CSP)
        response.headers.setdefault("X-Content-Type-Options", "nosniff")
        response.headers.setdefault("Referrer-Policy", "no-referrer")
        response.headers.setdefault("X-Frame-Options", "DENY")
        if path.startswith("/api/"):
            response.headers.setdefault("Cache-Control", "no-store")
        return response

    app.add_middleware(_BodyLimit, max_bytes=settings.max_body_bytes)

    @app.exception_handler(RequestValidationError)
    async def _validation(request: Request, exc: RequestValidationError):
        msgs = []
        for err in exc.errors():
            loc = ".".join(str(x) for x in err.get("loc", ()) if x != "body")
            msgs.append(f"{loc}: {err.get('msg')}" if loc else str(err.get("msg")))
        return JSONResponse({"detail": "; ".join(msgs)[:2000]}, status_code=422)

    def limit(request: Request, bucket: str, n: int, window: float, what: str) -> None:
        wait = limiter.check(bucket, _client_ip(request, settings), n, window)
        if wait is not None:
            raise HTTPException(429, f"rate limit: at most {n} {what}; try again in "
                                     f"{int(math.ceil(wait))} s",
                                headers={"Retry-After": str(int(math.ceil(wait)))})

    def chembl_client():
        from ..chembl import ChemblClient

        return ChemblClient(settings.chembl_cache_dir, offline=settings.offline,
                            release=settings.chembl_release, timeout=30.0, max_retries=2)

    def chembl_error(exc: Exception) -> HTTPException:
        from ..chembl import CacheMiss
        import requests

        if isinstance(exc, CacheMiss):
            return HTTPException(503, "this server runs offline and the ChEMBL data for the "
                                      "request is not in its cache")
        if isinstance(exc, requests.RequestException):
            return HTTPException(502, "the ChEMBL web service could not be reached; try again "
                                      "later")
        if isinstance(exc, (ValueError, KeyError)):
            return HTTPException(400, str(exc)[:1000])
        log.exception("unexpected ChEMBL error")
        return HTTPException(500, f"internal error ({type(exc).__name__}); please report it")

    def heading_index() -> tuple[str, list[dict]]:
        """(release, heading rows); memoised per release, rechecked hourly when online."""
        cached = state["headings"].get("current")
        if cached and (settings.offline or time.monotonic() - cached[2] < 3600):
            return cached[0], cached[1]
        with chembl_client() as client:
            rows = client.mesh_headings()
            release = client.release
        state["headings"]["current"] = (release, rows, time.monotonic())
        return release, rows

    class _HeadingSource:
        def __init__(self, rows):
            self._rows = rows

        def mesh_headings(self):
            return self._rows

    # ---------------------------------------------------------------- routes
    @app.get("/health")
    def health():
        w = state.get("worker")
        try:
            queued = store().count("queued")
            running = store().count("running")
            db = "ok"
        except Exception as exc:  # noqa: BLE001
            queued = running = None
            db = f"error: {type(exc).__name__}"
        ok = db == "ok" and (not settings.worker or (w is not None and w.is_alive()))
        return JSONResponse(dict(status="ok" if ok else "degraded", version=__version__,
                                 worker_alive=bool(w is not None and w.is_alive()),
                                 database=db, queued=queued, running=running),
                            status_code=200 if ok else 503)

    @app.get("/api/meta")
    def meta():
        from .. import expdata
        from ..chembl import CHEMBL_DATA_LICENCE
        from ..target_providers import RDKitTargetProvider

        rd = {s.name: s for s in RDKitTargetProvider(list(RDKIT_TARGETS)).specs()}
        measured = []
        for name, t in expdata.TARGETS.items():
            src = expdata.SOURCES[t.source]
            measured.append(dict(name=name, endpoint=t.endpoint, units=t.unit,
                                 source=src.title, citation=src.citation, doi=src.doi,
                                 licence=src.licence))
        return dict(
            version=__version__,
            mode=dict(name="exploratory", learner="selected", families=list(settings.families),
                      outer_splits=settings.outer_splits, repeats=len(settings.seeds),
                      inner_splits=settings.inner_splits, permutation=False, bootstrap=False,
                      feature_sets=list(jb.LADDER_SETS), reference="size",
                      note="Fast exploratory run. Confirmatory analyses use the command-line "
                           "protocol (topoqspr run, repeated nested CV, Y-randomisation, "
                           "bootstrap, multiplicity control)."),
            limits=dict(max_molecules=settings.max_molecules, min_molecules=settings.min_molecules,
                        max_indices_molecules=settings.max_indices_molecules,
                        max_body_bytes=settings.max_body_bytes, max_queue=settings.max_queue,
                        job_timeout_s=settings.job_timeout_s,
                        runs_per_hour=settings.rate_runs_per_hour),
            chembl=dict(release_pin=settings.chembl_release, offline=settings.offline,
                        licence=CHEMBL_DATA_LICENCE),
            mesh_allowlist=list(settings.mesh_allowlist) or None,
            computed_targets=[dict(name=n, target_class=rd[n].target_class, units=rd[n].units,
                                   description=rd[n].description) for n in RDKIT_TARGETS],
            measured_targets=measured,
            evidence_types=["ATC", "ClinicalTrials", "DailyMed", "EMA", "FDA", "INN", "USAN"],
        )

    @app.get("/api/headings")
    def headings(q: str = Query(min_length=2, max_length=100),
                 limit_n: int = Query(10, ge=1, le=25, alias="limit")):
        from ..cohort import suggest_headings

        if _CONTROL.search(q):
            raise HTTPException(422, "q contains control characters")
        try:
            release, rows = heading_index()
            df = suggest_headings(q, _HeadingSource(rows), limit=limit_n)
        except Exception as exc:  # noqa: BLE001 - mapped to a clear status
            raise chembl_error(exc) from None
        out = _dict_rows(df, ["mesh_id", "mesh_heading", "n_rows", "score", "matched_via"])
        allow = set(settings.mesh_allowlist)
        for r in out:
            r["allowed"] = not allow or r["mesh_id"] in allow
        return dict(release=release, query=q, headings=out)

    @app.post("/api/cohort")
    def cohort(req: CohortRequest, request: Request):
        from ..cohort import build_cohort

        limit(request, "cohort", settings.rate_cohort_per_min, 60.0, "cohort builds per minute")
        if settings.mesh_allowlist:
            blocked = [m for m in req.mesh_ids if m not in settings.mesh_allowlist]
            if blocked:
                raise HTTPException(403, f"MeSH IDs {blocked} are not enabled on this server "
                                         f"(enabled: {list(settings.mesh_allowlist)})")
        if not cohort_sem.acquire(timeout=60):
            raise HTTPException(503, "the server is busy building other cohorts; try again")
        try:
            with chembl_client() as client:
                c = build_cohort(req.mesh_ids, client, min_phase=req.min_phase,
                                 evidence_sources=req.evidence_sources,
                                 exclude_sole_evidence=req.exclude_sole_evidence,
                                 exclude_diagnostics=req.exclude_diagnostics,
                                 disease_label=req.disease_label)
        except HTTPException:
            raise
        except Exception as exc:  # noqa: BLE001
            raise chembl_error(exc) from None
        finally:
            cohort_sem.release()
        token = jb.new_id()
        cdir = settings.cohorts_dir / token
        created = jb.now_utc()
        c.write(cdir / "manifest.csv", timestamp_utc=created)
        side = json.loads((cdir / "manifest.json").read_text(encoding="utf-8"))
        n_inc = int(len(c.included))
        runnable, why = True, None
        if n_inc > settings.max_molecules:
            runnable, why = False, (f"{n_inc} included molecules; the web limit is "
                                    f"{settings.max_molecules} (use the command line)")
        elif n_inc < settings.min_molecules:
            runnable, why = False, (f"only {n_inc} included molecules; the web needs at least "
                                    f"{settings.min_molecules} for 5-fold cross-validation")
        info = dict(disease_label=c.disease_label, mesh_ids=c.parameters["mesh_ids"],
                    mesh_headings=side.get("mesh_headings"), chembl_release=c.release,
                    chembl_release_date=c.release_date, parameters=c.parameters,
                    n_included=n_inc, manifest_sha256=side["manifest_sha256"],
                    created_utc=created, runnable=runnable, not_runnable_reason=why)
        store().add_cohort(token, info, _client_ip(request, settings))
        m = c.manifest
        inc = m[m["decision"] == "include"]
        exc_ = m[m["decision"] != "include"]
        return dict(
            cohort_token=token, **info,
            counts=side["counts"], data_licence=side["data_licence"],
            max_molecules=settings.max_molecules,
            included=_dict_rows(inc, ["parent_chembl_id", "pref_name", "canonical_smiles",
                                      "standard_inchi_key", "phase_for_indication",
                                      "indication_ref_types"]),
            excluded=_dict_rows(exc_, ["molecule_chembl_id", "parent_chembl_id", "pref_name",
                                       "phase_for_indication", "reason"]),
        )

    def cohort_file(token: str, name: str) -> Path:
        if not jb.valid_id(token) or store().get_cohort(token) is None:
            raise HTTPException(404, "unknown cohort")
        p = settings.cohorts_dir / token / name
        if not p.is_file():
            raise HTTPException(404, "cohort files no longer available")
        return p

    @app.get("/api/cohorts/{token}/manifest.csv")
    def cohort_manifest(token: str):
        return FileResponse(cohort_file(token, "manifest.csv"), media_type="text/csv",
                            filename="manifest.csv")

    @app.post("/api/run", status_code=202)
    def run(req: RunRequest, request: Request):
        ip = _client_ip(request, settings)
        if not req.confirmed:
            raise HTTPException(400, "confirm that you have reviewed the cohort manifest first")
        if not jb.valid_id(req.cohort_token):
            raise HTTPException(422, "cohort_token is not a valid token")
        coh = store().get_cohort(req.cohort_token)
        if coh is None:
            raise HTTPException(404, "unknown or expired cohort; build it again")
        info = coh["info"]
        if not info.get("runnable"):
            raise HTTPException(400, info.get("not_runnable_reason") or "cohort cannot be run")
        if settings.mesh_allowlist and any(m not in settings.mesh_allowlist
                                           for m in info["mesh_ids"]):
            raise HTTPException(403, "this cohort's MeSH IDs are not enabled on this server")
        computed = list(dict.fromkeys(req.targets.computed))
        measured = list(dict.fromkeys(req.targets.measured))
        if not computed and not measured:
            raise HTTPException(422, "choose at least one target")
        src = settings.cohorts_dir / req.cohort_token
        if not (src / "manifest.csv").is_file():
            raise HTTPException(404, "cohort files no longer available; build it again")
        if store().count("queued", "running") >= settings.max_queue:
            raise HTTPException(503, "the job queue is full; try again later",
                                headers={"Retry-After": "120"})
        if store().active_for_ip(ip) >= settings.max_active_per_ip:
            raise HTTPException(429, "you already have a job queued or running; wait for it "
                                     "to finish")
        limit(request, "run", settings.rate_runs_per_hour, 3600.0, "runs per hour")
        job_id = jb.new_id()
        created = jb.now_utc()
        jdir = settings.jobs_dir / job_id
        (jdir / "cohort").mkdir(parents=True)
        shutil.copyfile(src / "manifest.csv", jdir / "cohort" / "manifest.csv")
        shutil.copyfile(src / "manifest.json", jdir / "cohort" / "manifest.json")
        targets = []
        if computed:
            targets.append({"provider": "rdkit", "names": computed, "options": {}})
        if measured:
            targets.append({"provider": "expdata", "names": measured, "options": {}})
        config = {
            "name": f"web-exploratory-{job_id[:8]}",
            "cohort": {"source": "manifest", "path": "cohort/manifest.csv"},
            "feature_sets": list(jb.LADDER_SETS),
            "targets": targets,
            "evaluation": {"learners": ["selected"], "families": list(settings.families),
                           "outer_splits": settings.outer_splits, "seeds": list(settings.seeds),
                           "inner_splits": settings.inner_splits},
            "runtime": {"n_jobs": settings.n_jobs, "chembl_cache_dir": settings.chembl_cache_dir,
                        "offline": settings.offline,
                        "expdata_cache_dir": settings.expdata_cache_dir},
        }
        from ..pipeline import resolve_config

        try:
            resolve_config(config)  # fail now, not in the worker
        except (ValueError, KeyError) as exc:
            shutil.rmtree(jdir, ignore_errors=True)
            raise HTTPException(400, str(exc)[:1000]) from None
        cohort_view = dict(token=req.cohort_token, disease_label=info["disease_label"],
                           mesh_ids=info["mesh_ids"], chembl_release=info["chembl_release"],
                           n_included=info["n_included"],
                           manifest_sha256=info["manifest_sha256"])
        request_rec = dict(created_utc=created, mode=req.mode,
                           targets=dict(computed=computed, measured=measured), cohort=cohort_view)
        (jdir / "config.json").write_text(json.dumps(config, sort_keys=True, indent=1) + "\n")
        (jdir / "request.json").write_text(json.dumps(request_rec, sort_keys=True, indent=1) + "\n")
        store().create_job(job_id, req.cohort_token, request_rec, ip, created)
        job = store().get_job(job_id)
        return _job_view(job)

    def get_job_or_404(job_id: str) -> dict:
        if not jb.valid_id(job_id):
            raise HTTPException(404, "unknown job")
        job = store().get_job(job_id)
        if job is None:
            raise HTTPException(404, "unknown job")
        return job

    @app.get("/api/jobs/{job_id}")
    def job_status(job_id: str):
        return _job_view(get_job_or_404(job_id))

    def done_file(job_id: str, name: str) -> Path:
        job = get_job_or_404(job_id)
        if job["status"] != "done":
            raise HTTPException(409, f"job is {job['status']}; results are available when it "
                                     "is done")
        path = {"ladder": settings.jobs_dir / job_id / "ladder.json",
                "zip": settings.jobs_dir / job_id / "results.zip",
                "record": settings.jobs_dir / job_id / "results" / "run_record.json",
                "manifest": settings.jobs_dir / job_id / "results" / "manifest.csv"}[name]
        root = settings.jobs_dir.resolve()
        if not path.resolve().is_relative_to(root) or not path.is_file():
            raise HTTPException(410, "result files are no longer available")
        return path

    @app.get("/api/jobs/{job_id}/results")
    def job_results(job_id: str):
        data = json.loads(done_file(job_id, "ladder").read_text(encoding="utf-8"))
        data["job"] = _job_view(get_job_or_404(job_id))
        return data

    @app.get("/api/jobs/{job_id}/results.zip")
    def job_zip(job_id: str):
        return FileResponse(done_file(job_id, "zip"), media_type="application/zip",
                            filename=f"topoqspr-{job_id[:8]}.zip")

    @app.get("/api/jobs/{job_id}/run_record.json")
    def job_record(job_id: str):
        return FileResponse(done_file(job_id, "record"), media_type="application/json",
                            filename="run_record.json")

    @app.get("/api/jobs/{job_id}/manifest.csv")
    def job_manifest(job_id: str):
        return FileResponse(done_file(job_id, "manifest"), media_type="text/csv",
                            filename="manifest.csv")

    @app.post("/api/indices")
    def indices(req: IndicesRequest, request: Request):
        import pandas as pd
        from rdkit import rdBase

        from ..graph import standardize
        from ..pipeline import compute_index_table

        limit(request, "indices", settings.rate_indices_per_min, 60.0,
              "index calculations per minute")
        try:
            pairs = parse_smiles_text(req.text, settings.max_indices_molecules,
                                      settings.max_smiles_length)
        except ValueError as exc:
            raise HTTPException(422, str(exc)) from None
        if not indices_sem.acquire(timeout=30):
            raise HTTPException(503, "the server is busy; try again in a moment")
        try:
            rows, graphs, ok_ids = [], [], []
            with rdBase.BlockLogs():
                for rid, smi in pairs:
                    res = standardize(smi)
                    rows.append(dict(id=rid, input_smiles=smi,
                                     status="ok" if res.ok else "rejected",
                                     reason="" if res.ok else res.rejection.code,
                                     smiles=res.canonical_smiles or "",
                                     inchikey=res.inchikey or ""))
                    if res.ok:
                        graphs.append(res.graph)
                        ok_ids.append(rid)
                vals, failures = compute_index_table(graphs, ok_ids, exact=True)
        finally:
            indices_sem.release()
        meta_df = pd.DataFrame(rows).set_index("id")
        out = meta_df.join(vals.astype(object).reindex(meta_df.index)).reset_index()
        fail = pd.DataFrame([dict(row=r["id"], feature="*", reason=f"rejected: {r['reason']}")
                             for r in rows if r["status"] != "ok"] + failures,
                            columns=["row", "feature", "reason"])
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            def add(name: str, data: str) -> None:
                info = zipfile.ZipInfo(f"topoqspr-indices/{name}", date_time=(2000, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o644 << 16
                zf.writestr(info, data.encode("utf-8"))
            add("indices.csv", out.to_csv(index=False, lineterminator="\n"))
            add("indices.failures.csv", fail.to_csv(index=False, lineterminator="\n"))
            add("README.txt",
                "TopoQSPR topological indices (v" + __version__ + ")\n\n"
                "indices.csv: standardisation status and every verified index per molecule "
                "(empty cell = undefined; never 0). indices.failures.csv: the reason for every "
                "undefined value or rejected structure. Definitions, references and verification "
                "status: docs/INDEX_VERIFICATION.md in the TopoQSPR repository (MIT licence).\n")
        n_ok = len(ok_ids)
        return Response(buf.getvalue(), media_type="application/zip",
                        headers={"Content-Disposition": 'attachment; filename="topoqspr-indices.zip"',
                                 "X-Molecules": str(len(pairs)), "X-Molecules-OK": str(n_ok),
                                 "X-Undefined-Values": str(len(failures))})

    # ---------------------------------------------------------------- front-end
    @app.get("/", include_in_schema=False)
    def index():
        return FileResponse(STATIC_DIR / "index.html", media_type="text/html")

    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
    return app


def _lazy_app() -> FastAPI:
    return create_app(Settings.from_env())


#: ASGI application for ``uvicorn topoqspr.web.app:app`` (settings from the environment).
app = _lazy_app()
