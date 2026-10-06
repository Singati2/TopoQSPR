"""Why linear learners exploded on a few molecules: how far each flagged molecule lies from the
other molecules of its own target set, feature by feature (z-scores computed without the molecule,
as a training fold that excludes it would standardise).

Usage: python tools/check_extrapolation.py RESULTS_DIR OUT_DIR [--sd-factor 10]

Reads features.csv, targets.csv, predictions.csv (learner ``selected``) of a run; a molecule is
"flagged" when one of its out-of-fold predictions is more than ``sd_factor`` standard deviations of
the target away from the truth (scale-free, so melting points in degrees and log values are treated alike).  Writes extrapolation.csv and extrap_macros.tex (\\extrapMaxZ, \\extrapIndex, ...).
Nothing is typed by hand.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("out")
    ap.add_argument("--sd-factor", type=float, default=10.0,
                    help="a prediction is 'exploded' when |residual| exceeds this many SDs of the target")
    a = ap.parse_args()
    res, out = Path(a.results), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    f = pd.read_csv(res / "features.csv", keep_default_na=False, na_values=[""]).set_index("id")
    t = pd.read_csv(res / "targets.csv", keep_default_na=False, na_values=[""]).set_index("id")
    p = pd.read_csv(res / "predictions.csv", keep_default_na=False, na_values=[""])
    man = pd.read_csv(res / "manifest.csv", keep_default_na=False, na_values=[""])
    names = man.drop_duplicates("parent_chembl_id").set_index("parent_chembl_id")["pref_name"].to_dict()
    p = p[p["learner"] == "selected"].copy()
    p["res"] = (p["y_pred"] - p["y_true"]).abs()
    sd = p.groupby("target")["y_true"].transform("std")
    flagged = p[p["res"] > a.sd_factor * sd].groupby(["target", "row"])["res"].max().reset_index()
    ti_cols = [c for c in f.columns if c.startswith("ti.")]
    rows = []
    for _, r in flagged.iterrows():
        target, cid = r["target"], r["row"]
        ids = t.index[t[target].notna()]
        sub = f.loc[ids, ti_cols].astype(float)
        rest = sub.drop(index=cid)
        z = ((sub.loc[cid] - rest.mean()) / rest.std(ddof=0)).abs().sort_values(ascending=False)
        n_heavy = float(f.loc[cid, "size.n_heavy"])
        rows.append(dict(target=target, molecule=cid, name=str(names.get(cid, "")).title(),
                         n_heavy=int(n_heavy), size_rank=int((f.loc[ids, "size.n_heavy"] >= n_heavy).sum()), n=int(len(ids)),
                         worst_residual=float(r["res"]), top_index=z.index[0][3:], top_z=float(z.iloc[0]),
                         z_Z=float(z.get("ti.Z", float("nan"))), n_indices_over_5=int((z > 5).sum())))
    tab = pd.DataFrame(rows).sort_values(["target", "size_rank"])
    tab.to_csv(out / "extrapolation.csv", index=False)
    if len(tab):
        worst = tab.sort_values("top_z", ascending=False).iloc[0]
        macros = {
            "extrapNMolecules": int(tab["molecule"].nunique()),
            "extrapNames": ", ".join(dict.fromkeys(tab.sort_values("n_heavy", ascending=False)["name"])),
            "extrapMinHeavy": int(tab["n_heavy"].min()), "extrapMaxHeavy": int(tab["n_heavy"].max()),
            "extrapAllLargest": "yes" if (tab["size_rank"] <= 2).all() else "no",
            "extrapMaxZ": f"{worst['top_z']:,.0f}", "extrapMaxZIndex": worst["top_index"].replace("_", "\\_"),
            "extrapMinZZ": f"{tab['z_Z'].min():.0f}", "extrapMaxResidual": f"{tab['worst_residual'].max():,.0f}",
        }
        (out / "extrap_macros.tex").write_text("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in macros.items()))
        print(tab.to_string(index=False))
        print(macros)
    else:
        print("no flagged molecules")


if __name__ == "__main__":
    main()
