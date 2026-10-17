"""Notebook 00: environment and data audit. PASS/FAIL with reasons; fails closed.

    python -I docs/research/scripts/restore_opera_values.py docs/research/data/E3 RESTORED
    python -I docs/research/scripts/audit_environment.py docs/research/data/E3 RESTORED OUT
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import repro_common as RC  # noqa: E402
from restore_opera_values import col_sha256  # noqa: E402

PRIVATE = "|".join(["/home/" + "mpcrlab", "/mnt/" + "beegfs", r"\." + "ven" + "vs/"])  # assembled so this line does not match itself
PLAUSIBLE = {"mp": (-200.0, 700.0), "logp": (-10.0, 15.0)}  # degC, log10 units


def main(data: Path, restored: Path, out: Path) -> RC.Checks:
    ck = RC.Checks("00_environment_and_data_audit")
    meta = RC.runtime_metadata()
    ck.info("runtime", json.dumps(meta))
    ck.require("git checkout clean (tracked files)", not meta["git_dirty"], meta["git_dirty"], False)
    mm = RC.version_mismatches()
    ck.require("installed versions equal the pinned set", not mm, mm or "none", "none")
    ck.require("Python >= 3.10", tuple(map(int, meta["python"].split(".")[:2])) >= (3, 10), meta["python"], ">= 3.10")
    ck.require("RAM >= 4 GB (E3 peak was about 1.1 GB per worker)", meta["ram_gb"] >= 4, meta["ram_gb"], ">= 4")
    # public sources: pinned URL + sha256 (topoqspr.expdata raises on mismatch)
    try:
        src = RC.prefetch_sources()
        for k, v in src.items():
            ck.require(f"source {k}: downloaded/cached file matches pinned sha256", True, v["sha256"][:16], v["url"])
            ck.info(f"source {k} licence", f"{v['licence']} (redistribution={v['redistribution']})")
    except Exception as exc:  # noqa: BLE001 - reported as FAIL
        ck.require("public sources reachable and sha256-verified", False, repr(exc)[:200], "pinned files")
    # published dataset
    n, bad = RC.check_sha256sums(data)
    ck.require(f"published E3 files match SHA256SUMS ({n} files)", not bad, bad, "[]")
    rec = json.loads((data / "removed_columns.json").read_text())["removed"]
    for rel, m in rec.items():
        P = pd.read_csv(restored / rel.replace("predictions_no_y", "predictions"), float_precision="round_trip",
                        index_col=0 if m["column"] == "y" else None)
        col = P[m["column"]]
        ck.require(f"restored {rel}: row count", len(P) == m["n_rows"], len(P), m["n_rows"])
        ck.require(f"restored {rel}: no unmatched identifiers", not col.isna().any(), int(col.isna().sum()), 0)
        ck.require(f"restored {rel}: value hash equals the recorded hash", col_sha256(col) == m["sha256_float64_in_row_order"],
                   col_sha256(col)[:16], m["sha256_float64_in_row_order"][:16])
        lo, hi = PLAUSIBLE["mp" if rel.startswith("external") else "logp"]
        ck.require(f"restored {rel}: values within plausible unit range [{lo}, {hi}]",
                   bool(col.between(lo, hi).all()), f"[{col.min():.2f}, {col.max():.2f}]", f"[{lo}, {hi}]")
    for ep, kind in (("mp_bradley", "mp"), ("logp_opera", "logp")):
        inc = pd.read_csv(data / f"inclusion_{ep}.csv", index_col=0)
        RC.check_ids(ck, f"inclusion {ep}", inc.index)
        ck.require(f"inclusion {ep}: constitution = InChIKey first block", bool((inc.constitution == inc.index.str[:14]).all()), "", "")
        ck.require(f"inclusion {ep}: every formula has >= 5 constitutions (E3 rule)",
                   bool((inc.drop_duplicates("constitution").groupby("formula").size() >= 5).all()), "", "")
    P = pd.read_csv(data / "cells/mp_bradley_E3a/predictions.csv.gz", float_precision="round_trip")
    lo, hi = PLAUSIBLE["mp"]
    ck.require("Bradley y_true within plausible degC range", bool(P.y_true.between(lo, hi).all()),
               f"[{P.y_true.min()}, {P.y_true.max()}]", f"[{lo}, {hi}]")
    split = pd.read_csv(data / "external/test_split.csv")
    RC.check_ids(ck, "external test split", split.id)
    sizes = split.test.value_counts().to_dict()
    ck.require("external split sizes 1,014 / 7,327", sizes == {"unseen_formula": 7327, "known_formula": 1014}, sizes, "1014/7327")
    inc = pd.read_csv(data / "inclusion_mp_bradley.csv", index_col=0)
    ck.info("external test constitutions also in the internal MP inclusion set (expected: the external "
            "train set is the PHYSPROP-free subset only)", int(split.id.str[:14].isin(inc.constitution).sum()))
    # no private paths in the reproducibility code
    hits = []
    for f in list((RC.REPO / "docs/research/scripts").glob("*.py")) + list((RC.REPO / "notebooks/colab").glob("*.ipynb")):
        for i, line in enumerate(f.read_text(errors="replace").splitlines(), 1):
            if re.search(PRIVATE, line) and not line.lstrip().startswith("#"):
                hits.append(f"{f.name}:{i}")
    ck.require("no workstation paths in scripts or notebooks", not hits, hits[:5], "[]")
    ck.finish(out)
    return ck


if __name__ == "__main__":
    main(*(Path(p).expanduser() for p in sys.argv[1:4]))
