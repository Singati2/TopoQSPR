"""Restore the OPERA/PHYSPROP measured values removed from the public E3 dataset, and verify them.

Downloads the pinned OPERA_Data.zip (github.com/kmansouri/OPERA, MIT) through topoqspr.expdata,
which checks its sha256, re-aggregates per parent InChIKey exactly as E3 did, joins by InChIKey,
and checks the restored column against the sha256 recorded in removed_columns.json. Writes a
directory with the same layout as the private run archive, so any E3 / P0 script can run on it.
Bradley (CC0) files are copied unchanged.

    pip install -e .            # from the repository root (Colab: see notebooks/)
    python -I docs/research/scripts/restore_opera_values.py docs/research/data/E3 /content/e3_restored
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from topoqspr import expdata


def col_sha256(s: pd.Series) -> str:
    return hashlib.sha256(np.ascontiguousarray(s.to_numpy(dtype=np.float64)).tobytes()).hexdigest()


def main(src: Path, dst: Path) -> None:
    rec = json.loads((src / "removed_columns.json").read_text())["removed"]
    shutil.copytree(src, dst, dirs_exist_ok=True)
    values = {}
    for name in ("logp_opera", "mp_opera"):
        agg = expdata.get_target(name, drop_disagreeing=True)
        values[name] = agg.set_index("inchikey")["value"].astype(float)
    ok = True
    for rel, meta in rec.items():
        P = pd.read_csv(src / rel, float_precision="round_trip")
        endpoint = "mp_opera" if rel.startswith("external/") else "logp_opera"
        keycol = "id" if meta["column"] == "y" else "row"
        y = P[keycol].map(values[endpoint])
        P.insert(len(P.columns) - (1 if "formula" in P.columns else 0) if meta["column"] == "y"
                 else list(P.columns).index("y_pred"), meta["column"], y)
        got = col_sha256(P[meta["column"]])
        match = got == meta["sha256_float64_in_row_order"] and len(P) == meta["n_rows"]
        ok &= match
        print(f"{'OK  ' if match else 'FAIL'} {rel}: {y.isna().sum()} unmatched, sha256 {got[:16]}…")
        out = dst / rel.replace("predictions_no_y", "predictions")
        if meta["column"] == "y":
            P.set_index("id").to_csv(out)
        else:
            P.to_csv(out, index=False)
        (dst / rel).unlink()
    if not ok:
        sys.exit("restored values do not match the recorded hashes (different OPERA version or "
                 "aggregation); do not use them to verify E3")
    print(f"all restored columns verified; archive layout written to {dst}")


if __name__ == "__main__":
    main(Path(sys.argv[1]).expanduser(), Path(sys.argv[2]).expanduser())
