"""Audit 2026-10-12, B2: trace every catastrophic out-of-fold prediction to a molecule.

Reads an archived run directory (predictions.csv, features.csv, molecules.csv,
folds.csv) and, for each (target, set) whose selected-learner pooled R^2 is
below -1 in any fold, lists the test molecules with the largest absolute
errors and, for each, the feature columns that lie furthest outside the
training fold's range (robust z = (x - median_train) / MAD_train).

Read-only: nothing in the run directory is modified.

    python -I docs/audit-2026-10-12/scripts/catastrophic_predictions.py \
        ~/topoqspr_runs/main_confirmatory docs/audit-2026-10-12/out
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PREFIX = {"TI": ("ti.",), "TI_composition": ("ti.", "composition."),
          "rdkit2d": ("rdkit2d.",), "composition": ("composition.",), "size": ("size.",)}


def main(run: Path, out: Path) -> None:
    folds = pd.read_csv(run / "folds.csv")
    pred = pd.read_csv(run / "predictions.csv")
    feats = pd.read_csv(run / "features.csv", index_col="id")
    mols = pd.read_csv(run / "molecules.csv", index_col="id")
    sel = folds[(folds.learner == "selected") & (folds.r2 < -1)]
    rows = []
    for f in sel.itertuples(index=False):
        p = pred[(pred.target == f.target) & (pred.learner == "selected") & (pred.set == f.set)
                 & (pred.repeat == f.repeat)]
        test = p[p.fold == f.fold].copy()
        train_ids = p[p.fold != f.fold].row.values  # every non-test row of this repeat
        test["abs_err"] = (test.y_pred - test.y_true).abs()
        cols = [c for c in feats.columns if c.startswith(PREFIX[f.set])]
        tr = feats.loc[train_ids, cols]
        med = tr.median()
        mad = (tr - med).abs().median().replace(0, np.nan)
        sse = float(((test.y_pred - test.y_true) ** 2).sum())
        for _, r in test.nlargest(3, "abs_err").iterrows():
            z = ((feats.loc[r.row, cols] - med) / mad).abs().sort_values(ascending=False)
            above = (feats.loc[r.row, cols] > tr.max()).sum()
            rows.append(dict(
                target=f.target, set=f.set, repeat=f.repeat, fold=f.fold, family=f.family,
                params=f.params, fold_r2=f.r2, molecule=r.row, smiles=mols.loc[r.row, "smiles"],
                n_heavy=mols.loc[r.row, "n_heavy"], y_true=r.y_true, y_pred=r.y_pred,
                abs_err=r.abs_err, share_of_fold_sse=float((r.y_pred - r.y_true) ** 2 / sse),
                in_domain=r.in_domain, n_cols_above_train_max=int(above),
                top_feature=z.index[0], top_robust_z=float(z.iloc[0]),
                second_feature=z.index[1], second_robust_z=float(z.iloc[1])))
    df = pd.DataFrame(rows)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / f"catastrophic_{run.name}.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 40):
        print(df[["target", "set", "repeat", "fold", "fold_r2", "molecule", "n_heavy", "y_true",
                  "y_pred", "share_of_fold_sse", "in_domain", "top_feature", "top_robust_z"]])


if __name__ == "__main__":
    main(Path(sys.argv[1]).expanduser(), Path(sys.argv[2]))
