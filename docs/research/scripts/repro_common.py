"""Shared helpers for the Colab notebooks in notebooks/colab/ and their command-line twins.

Everything here is plumbing: runtime metadata, checksum checks, atomic writes, a PASS/FAIL
collector that fails closed, and a few integrity checks. No modelling code lives here, so the
notebooks and a Linux/HPC run go through the same E3 functions (e3_isomer.py, topoqspr.evaluate).
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
E3_DATA = REPO / "docs/research/data/E3"

#: The dependency set that produced the committed E3 results (run_info.json) plus the ones
#: whose versions decide InChIKeys and therefore the OPERA restore.
PINNED = {"rdkit": "2026.03.6", "chembl_structure_pipeline": "1.2.4", "numpy": "2.2.6",
          "pandas": "2.3.3", "scipy": "1.15.3", "scikit-learn": "1.7.2", "joblib": "1.6.0",
          "xgboost": "3.2.0", "openpyxl": "3.1.5"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_sha(repo: Path = REPO) -> str:
    r = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=repo)
    return r.stdout.strip() or "unknown"


def git_dirty(repo: Path = REPO) -> bool:
    r = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"],
                       capture_output=True, text=True, cwd=repo)
    return bool(r.stdout.strip())


def installed_versions() -> dict:
    from importlib.metadata import PackageNotFoundError, version
    out = {}
    for k in PINNED:
        try:
            out[k] = version(k)
        except PackageNotFoundError:
            out[k] = None
    out["rdkit"] = __import__("rdkit").__version__  # the wheel's own string, e.g. 2026.03.6
    return out


def runtime_metadata() -> dict:
    try:
        import psutil  # present on Colab
        ram_gb = round(psutil.virtual_memory().total / 2**30, 1)
    except ImportError:
        ram_gb = round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**30, 1)
    return dict(utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                python=platform.python_version(), platform=platform.platform(),
                cpus=os.cpu_count(), ram_gb=ram_gb, colab="google.colab" in sys.modules,
                git_sha=git_sha(), git_dirty=git_dirty(), versions=installed_versions(),
                topoqspr_cache=os.environ.get("TOPOQSPR_CACHE", "(default ~/.cache/topoqspr)"))


def version_mismatches() -> dict:
    got = installed_versions()
    norm = lambda v: None if v is None else ".".join(str(int(p)) if p.isdigit() else p
                                                     for p in str(v).split("."))
    return {k: (v, got.get(k)) for k, v in PINNED.items() if norm(got.get(k)) != norm(v)}


# --------------------------------------------------------------------------- checks
@dataclass
class Checks:
    """PASS/FAIL collector. ``require`` records a failure; ``finish`` raises if any failed.

    Failures are never replaced by expected values: the observed value is stored as found.
    """

    name: str
    rows: list = field(default_factory=list)

    def require(self, check: str, ok: bool, observed="", expected="", note="") -> bool:
        self.rows.append(dict(check=check, status="PASS" if ok else "FAIL",
                              observed=str(observed)[:300], expected=str(expected)[:300], note=note))
        print(f"{'PASS' if ok else 'FAIL'}  {check}" + ("" if ok else f"  (observed {observed!s:.120}, expected {expected!s:.120})"),
              flush=True)
        return ok

    def info(self, check: str, observed="", note=""):
        self.rows.append(dict(check=check, status="INFO", observed=str(observed)[:300], expected="", note=note))
        print(f"INFO  {check}: {observed!s:.160}", flush=True)

    def close(self, check: str, a: float, b: float, tol: float, note="") -> bool:
        d = abs(float(a) - float(b)) if np.isfinite(a) and np.isfinite(b) else (0.0 if (np.isnan(a) and np.isnan(b)) else np.inf)
        return self.require(check, d <= tol, f"{a!r} (|diff| {d:.2e})", f"{b!r} ± {tol:g}", note)

    @property
    def table(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows, columns=["check", "status", "observed", "expected", "note"])

    @property
    def n_fail(self) -> int:
        return int((self.table.status == "FAIL").sum())

    def save(self, out: Path) -> Path:
        out = Path(out); out.mkdir(parents=True, exist_ok=True)
        p = out / f"checks_{self.name}.csv"
        atomic_write_csv(self.table, p, index=False)
        atomic_write_text(json.dumps(dict(name=self.name, n_pass=int((self.table.status == "PASS").sum()),
                                          n_fail=self.n_fail, runtime=runtime_metadata()), indent=1),
                          out / f"checks_{self.name}.json")
        return p

    def finish(self, out: Path | None = None) -> None:
        if out is not None:
            self.save(out)
        n_pass = int((self.table.status == "PASS").sum())
        print(f"\n== {self.name}: {n_pass} PASS, {self.n_fail} FAIL ==", flush=True)
        if self.n_fail:
            print(self.table[self.table.status == "FAIL"].to_string(index=False))
            raise SystemExit(f"{self.name}: {self.n_fail} check(s) FAILED; results of this run must not be used")


def check_sha256sums(root: Path, sums_name: str = "SHA256SUMS") -> tuple[int, list]:
    lines = (Path(root) / sums_name).read_text().splitlines()
    bad = []
    for line in lines:
        h, rel = line.split("  ", 1)
        p = Path(root) / rel
        if not p.exists() or sha256_file(p) != h:
            bad.append(rel)
    return len(lines), bad


# --------------------------------------------------------------------------- io
def atomic_write_text(text: str, path: Path) -> None:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp{os.getpid()}")
    tmp.write_text(text)
    os.replace(tmp, path)


def atomic_write_csv(df: pd.DataFrame, path: Path, **kw) -> None:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp{os.getpid()}")
    df.to_csv(tmp, **kw)
    os.replace(tmp, path)


# --------------------------------------------------------------------------- integrity
def check_ids(ck: Checks, label: str, ids, expected=None) -> None:
    s = pd.Index(ids)
    ck.require(f"{label}: no missing ids", not s.isna().any(), int(s.isna().sum()), 0)
    ck.require(f"{label}: ids unique", s.is_unique, int(s.duplicated().sum()), 0)
    if expected is not None:
        e = pd.Index(expected)
        ck.require(f"{label}: id set equals the archived set", set(s) == set(e),
                   f"{len(set(s) - set(e))} extra, {len(set(e) - set(s))} missing", "0 extra, 0 missing")


def check_group_disjoint(ck: Checks, label: str, train_groups, test_groups) -> None:
    inter = set(train_groups) & set(test_groups)
    ck.require(f"{label}: train/test groups disjoint", not inter, len(inter), 0)


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (independent implementation)."""
    p = np.asarray(p, float); m = p.size; o = np.argsort(p); q = np.empty(m); run = 1.0
    for j in range(m - 1, -1, -1):
        run = min(run, p[o[j]] * m / (j + 1)); q[o[j]] = run
    return q


def prefetch_sources(keys=("bradley", "opera")) -> dict:
    """Download the pinned public source files (sha256-checked by topoqspr.expdata) into the
    cache, so the E3 builders, which read the cache offline, can run on a fresh runtime."""
    from topoqspr import expdata
    out = {}
    for k in keys:
        src = expdata.SOURCES[k]
        p = expdata.fetch(k)
        out[k] = dict(path=str(p), url=src.url, sha256=src.sha256, licence=src.licence,
                      redistribution=src.redistribution)
    return out
