"""Fixtures: a throw-away Git repository with small synthetic (non-sensitive) E3-shaped artifacts."""
from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

import pytest

from researchos import config
from researchos.store import Store

SUMMARY = """test,set,n,n_formulas,family,r2,mae
unseen_formula,composition,100,40,svr,0.20,70.0
unseen_formula,A_comp,100,40,extratrees,0.47,58.0
unseen_formula,B_comp,100,40,svr,0.42,60.0
"""
PAIRED = """test,set,reference,d_r2,d_r2_ci_low,d_r2_ci_high
unseen_formula,A_comp,composition,0.26999999999999996,0.2,0.3
unseen_formula,B_comp,composition,0.21999999999999997,0.15,0.28
"""
RECOMP = """test,n,set,reference,d_r2,lo,hi
unseen_formula,100,A_comp,composition,0.26999999999999996,0.2,0.3
unseen_formula,100,B_comp,composition,0.21999999999999997,0.15,0.28
"""
INTERNAL = """endpoint,estimand,set,n,p,r2_fold_mean,family_counts
mp,E3b,A_comp,100,10,0.60,"{""svr"": 5}"
mp,E3b,B_comp,100,10,0.61,"{""svr"": 5}"
"""
FOLDS = "set,repeat,fold,n_train,n_test,r2\n" + "".join(
    f"{s},{r},{f},80,20,{0.6 + (0.01 if s == 'A_comp' else 0) + 0.001 * (r * 5 + f)}\n"
    for s in ("A_comp", "B_comp") for r in range(3) for f in range(5))


def _git(repo, *a):
    subprocess.run(["git", "-C", str(repo), *a], check=True, capture_output=True)


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


@pytest.fixture
def repo(tmp_path) -> Path:
    r = tmp_path / "repo"
    d = r / "out" / "ext"
    d.mkdir(parents=True)
    (d / "summary.csv").write_text(SUMMARY)
    (d / "paired.csv").write_text(PAIRED)
    (d / "SHA256SUMS").write_text(f"{sha(SUMMARY)}  summary.csv\n{sha(PAIRED)}  paired.csv\n"
                                  f"{'0' * 64}  predictions.csv.gz\n")
    (r / "out" / "recomputed.csv").write_text(RECOMP)
    (r / "out" / "internal.csv").write_text(INTERNAL)
    (r / "out" / "folds.csv").write_text(FOLDS)
    (r / ".gitignore").write_text(".researchos/\n")
    _git(r, "init", "-q", "-b", "main")
    _git(r, "-c", "user.name=t", "-c", "user.email=t@example.org", "add", ".")
    _git(r, "-c", "user.name=t", "-c", "user.email=t@example.org", "commit", "-q", "-m", "fixture")
    return r


@pytest.fixture
def cfg(repo, tmp_path):
    return config.load(repo, tmp_path / "state")


@pytest.fixture
def store(cfg):
    return Store(cfg.db_path)
