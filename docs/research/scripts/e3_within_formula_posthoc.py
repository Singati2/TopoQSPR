"""POST HOC, DESCRIPTIVE (not in AMD-2026-10-13-E3): within-formula skill of each feature set
on MP E3a. Within each (repeat, outer fold, test formula with >= 2 test members), centre both
y and y_hat on their means, then pool: R2_c = 1 - SS(yc - pc) / SS(yc), and corr(yc, pc).
Composition is constant within a formula, so its R2_c is 0 by construction. The pre-declared
'within-formula R2' (06 §5) mixes formula-level bias with within-formula ranking; this metric
isolates the ranking. Uses test labels only for centring the metric, never for prediction.

    python -I docs/research/scripts/e3_within_formula_posthoc.py ~/topoqspr_runs/research_2026-10-13_E3
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

run = Path(sys.argv[1]).expanduser()
out = Path(__file__).resolve().parents[1] / "out" / "E3"
out.mkdir(parents=True, exist_ok=True)
for cell in ("mp_bradley_E3a", "mp_bradley_E3b"):
    f = run / "cells" / cell / "predictions.csv.gz"
    if not f.exists():
        continue
    P = pd.read_csv(f, float_precision="round_trip")
    meta = pd.read_csv(run / "inclusion_mp_bradley.csv", index_col=0)
    P["f"] = meta.loc[P.row, "formula"].values
    P = P[P.groupby(["set", "repeat", "fold", "f"]).row.transform("size") >= 2]
    rows = []
    for s, d in P.groupby("set"):
        g = d.groupby(["repeat", "fold", "f"])
        yc = d.y_true - g.y_true.transform("mean"); pc = d.y_pred - g.y_pred.transform("mean")
        rows.append(dict(cell=cell, set=s, n=len(d), within_R2_centred=1 - ((yc - pc) ** 2).sum() / (yc ** 2).sum(),
                         within_corr=float(np.corrcoef(yc, pc)[0, 1]),
                         within_MAE_centred=float((yc - pc).abs().mean())))
    df = pd.DataFrame(rows).sort_values("within_R2_centred")
    df.to_csv(out / f"within_formula_posthoc_{cell}.csv", index=False)
    print(df.round(3).to_string(index=False))
