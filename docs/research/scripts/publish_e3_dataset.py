"""Build the public E3 dataset under docs/research/data/E3 from the local run archive.

Licensing (docs/research/01b_DATA_LICENCES.md, 10_LICENSING_RISK_REVIEW.md):
- Bradley Open Melting Point (CC0): published in full, measured values included.
- OPERA / PHYSPROP (logP cells, external MP test): data terms upstream are unknown, so the measured
  values (y_true / y) are REMOVED. Identifiers, formulas, folds and model predictions are kept, and
  sha256 of the removed value column is recorded so `restore_opera_values.py` can verify a restore.

    python -I docs/research/scripts/publish_e3_dataset.py ~/topoqspr_runs/research_2026-10-13_E3
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def col_sha256(s: pd.Series) -> str:
    return hashlib.sha256(np.ascontiguousarray(s.to_numpy(dtype=np.float64)).tobytes()).hexdigest()


def file_sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main(run: Path) -> None:
    repo = Path(__file__).resolve().parents[3]
    out = repo / "docs" / "research" / "data" / "E3"
    out.mkdir(parents=True, exist_ok=True)
    removed = {}
    originals = {}
    for name in ("run_info.json", "inclusion_mp_bradley.csv", "inclusion_logp_opera.csv"):
        shutil.copy2(run / name, out / name)
    for cdir in sorted((run / "cells").glob("*_E3*")):
        d = out / "cells" / cdir.name
        d.mkdir(parents=True, exist_ok=True)
        shutil.copy2(cdir / "folds.csv", d / "folds.csv")
        originals[f"cells/{cdir.name}/predictions.csv.gz"] = file_sha256(cdir / "predictions.csv.gz")
        if cdir.name.startswith("mp_bradley"):
            shutil.copy2(cdir / "predictions.csv.gz", d / "predictions.csv.gz")
        else:
            P = pd.read_csv(cdir / "predictions.csv.gz", float_precision="round_trip")
            removed[f"cells/{cdir.name}/predictions_no_y.csv.gz"] = dict(
                column="y_true", key=["set", "repeat", "fold", "row"], n_rows=len(P),
                sha256_float64_in_row_order=col_sha256(P.y_true))
            P.drop(columns="y_true").to_csv(d / "predictions_no_y.csv.gz", index=False)
    ext = out / "external"
    ext.mkdir(exist_ok=True)
    shutil.copy2(run / "external" / "info.json", ext / "info.json")
    shutil.copy2(run / "external" / "summary.csv", ext / "summary.csv")
    originals["external/predictions.csv.gz"] = file_sha256(run / "external" / "predictions.csv.gz")
    pr = pd.read_csv(run / "external" / "predictions.csv.gz", float_precision="round_trip")
    removed["external/predictions_no_y.csv.gz"] = dict(
        column="y", key=["id"], n_rows=len(pr), sha256_float64_in_row_order=col_sha256(pr.y))
    pr.drop(columns="y").to_csv(ext / "predictions_no_y.csv.gz", index=False)
    # known/unseen split (formula seen in the PHYSPROP-free Bradley training set), so the external
    # analysis does not need the training set rebuilt. Taken from the P0 run, which rebuilt it.
    split = Path(__file__).resolve().parents[1] / "out" / "P0" / "external_split.csv"
    if split.exists():
        shutil.copy2(split, ext / "test_split.csv")
    (out / "removed_columns.json").write_text(json.dumps(dict(
        note="Measured OPERA/PHYSPROP values removed (upstream data terms unknown). "
             "Restore with docs/research/scripts/restore_opera_values.py.",
        removed=removed, sha256_of_original_archive_files=originals), indent=1))
    lines = [f"{file_sha256(p)}  {p.relative_to(out).as_posix()}"
             for p in sorted(out.rglob("*")) if p.is_file() and p.name != "SHA256SUMS"]
    (out / "SHA256SUMS").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(Path(sys.argv[1]).expanduser())
