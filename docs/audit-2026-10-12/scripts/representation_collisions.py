"""Audit 2026-10-12, section C: representation-limited empirical error floor.

For a feature map phi, molecules with identical phi(G) form one equivalence class C.
Any deterministic predictor f(phi(G)) gives one value per class, so on a finite set

    min_f (1/n) sum_i (y_i - f(phi_i))^2 = (1/n) sum_C sum_{i in C} (y_i - ybar_C)^2,

attained by f = class mean (squared-error minimiser of each class separately). Hence
R^2 <= 1 - floor / Var_n(y) for *every* model on *this* dataset (in-sample, empirical;
not an out-of-sample bound). This script computes, per target and per feature map
(TI = the 102 archived indices, TI+composition, composition, size), the classes,
the floor and the implied R^2 ceiling, and lists every multi-member class.

Equality is tested on the archived floating-point values rounded to 9 significant
digits (indices computed from the same graph are bitwise equal; rounding only guards
against summation-order noise). NaN compares equal to NaN.

    python -I docs/audit-2026-10-12/scripts/representation_collisions.py \
        ~/topoqspr_runs/main_confirmatory docs/audit-2026-10-12/out
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from topoqspr.experiments import index_table, prepare_structures

TARGETS = ["mp_bradley", "logp_opera", "logd74_az", "logs_aqsoldb"]


def keys(X: pd.DataFrame) -> pd.Series:
    r = X.apply(lambda c: c.map(lambda v: "nan" if pd.isna(v) else f"{v:.9g}"))
    return r.agg("|".join, axis=1)


def floor(y: pd.Series, k: pd.Series) -> dict:
    g = y.groupby(k)
    sse = float(((y - g.transform("mean")) ** 2).sum())
    var = float(((y - y.mean()) ** 2).sum())
    sizes = g.size()
    return dict(n=len(y), n_classes=int(len(sizes)), n_in_multi=int(sizes[sizes > 1].sum()),
                mse_floor=sse / len(y), r2_ceiling=1 - sse / var if var > 0 else np.nan)


def ethanol_vs_dimethyl_ether() -> pd.DataFrame:
    tab, _ = index_table(prepare_structures(["ethanol", "dimethyl_ether"], ["CCO", "COC"]))
    a, b = tab.loc["ethanol"], tab.loc["dimethyl_ether"]
    names = list(tab.columns)
    same = [n for n in names if (pd.isna(a[n]) and pd.isna(b[n])) or a[n] == b[n]]
    return pd.DataFrame(dict(n_indices=[len(names)], n_identical=[len(same)],
                             differing=[",".join(n for n in names if n not in same)]))


def main(run: Path, out: Path) -> None:
    F = pd.read_csv(run / "features.csv", index_col="id", float_precision="round_trip")
    Y = pd.read_csv(run / "targets.csv", index_col="id", float_precision="round_trip")
    M = pd.read_csv(run / "molecules.csv", index_col="id", float_precision="round_trip")
    blocks = {p: [c for c in F.columns if c.startswith(p)] for p in ("ti.", "composition.", "size.")}
    maps = {"size": blocks["size."][:2], "composition": blocks["composition."],
            "TI": blocks["ti."], "TI_composition": blocks["ti."] + blocks["composition."]}
    rows, members = [], []
    for t in TARGETS:
        ids = Y.index[Y[t].notna() & (M.loc[Y.index, "status"] == "ok")]
        y = Y.loc[ids, t].astype(float)
        for m, cols in maps.items():
            k = keys(F.loc[ids, cols])
            rows.append(dict(target=t, map=m, **floor(y, k)))
            if m in ("TI", "TI_composition"):
                for key, grp in k.groupby(k):
                    if len(grp) > 1:
                        ids_ = list(grp.index)
                        members.append(dict(target=t, map=m, size=len(ids_),
                                            ids=";".join(ids_),
                                            smiles=";".join(M.loc[ids_, "smiles"]),
                                            y=";".join(f"{v:g}" for v in y.loc[ids_]),
                                            y_range=float(y.loc[ids_].max() - y.loc[ids_].min())))
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "collision_floor.csv", index=False)
    pd.DataFrame(members).sort_values(["target", "map", "y_range"], ascending=[True, True, False]) \
        .to_csv(out / "collision_classes.csv", index=False)
    e = ethanol_vs_dimethyl_ether()
    e.to_csv(out / "ethanol_dimethyl_ether.csv", index=False)
    with pd.option_context("display.width", 200):
        print(df.round(4).to_string(index=False))
        print(e.to_string(index=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]).expanduser(), Path(sys.argv[2]))
