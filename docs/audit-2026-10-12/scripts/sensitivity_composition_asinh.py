"""Audit 2026-10-12, B1/B2/B4 sensitivity study (EXPLORATORY, not preregistered).

Question: do unlabeled-graph topological indices (TI) add out-of-fold accuracy over
elemental composition, and does the answer survive (i) a monotone, fit-free
transform that removes the unbounded-feature extrapolation failures and (ii) a
Bemis-Murcko scaffold split with the *same* 102-index registry?

Decisions fixed before any result of this script was seen (2026-10-12):
  * Features: the archived matrix of the frozen run (features.csv of
    main_confirmatory, 102-index registry at plan-freeze-v1 + 1 commit). Nothing is
    recomputed, so the raw sets must reproduce the archived summary exactly; the
    script asserts this for the random split.
  * Transform: asinh(x) applied elementwise to every TI column and every rdkit2d
    column; size and composition columns stay raw. asinh is monotone, sign-preserving,
    has no fitted parameters (so nothing leaks across folds) and is ~log(2x) for large
    x. It changes no descriptor definition; it is a model-input transform only.
  * Reference for the paired test: composition. Secondary paired comparisons against
    size are also written.
  * Splits: "random" (5-fold x seeds 0-4 grouped by InChIKey connectivity group, as
    archived) and "scaffold" (Bemis-Murcko groups, acyclic molecules singletons,
    connectivity-group members on their first member's scaffold, as in main_scaffold).
  * Learners: the frozen family list with in-fold selection ("selected"), inner 3-fold.
  * All rows are kept; no applicability-domain filtering of headline numbers.

    OMP_NUM_THREADS=1 python -I docs/audit-2026-10-12/scripts/sensitivity_composition_asinh.py \
        ~/topoqspr_runs/main_confirmatory ~/topoqspr_runs/audit_2026-10-12_sensitivity
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

from topoqspr import evaluate as ev
from topoqspr import stats

FAMILIES = ["ridge", "elasticnet", "pls", "knn", "svr", "extratrees"]
TARGETS = ["mp_bradley", "logp_opera", "logd74_az", "logs_aqsoldb"]


def feature_sets(F: pd.DataFrame) -> dict[str, pd.DataFrame]:
    def block(prefix):
        return F[[c for c in F.columns if c.startswith(prefix)]]

    size, comp, ti, rd = block("size."), block("composition."), block("ti."), block("rdkit2d.")
    size = size[["size.n_heavy", "size.n_bonds"]]  # the archived size set (p = 2)
    a_ti, a_rd = np.arcsinh(ti), np.arcsinh(rd)
    return {
        "size": size, "composition": comp,
        "TI": ti, "TI_composition": pd.concat([ti, comp], axis=1), "rdkit2d": rd,
        "TI_asinh": a_ti, "TI_asinh_composition": pd.concat([a_ti, comp], axis=1),
        "rdkit2d_asinh": a_rd,
    }


def paired(results, a, b, alpha=0.05):
    fa = results[a].folds.set_index(["repeat", "fold"])
    fb = results[b].folds.set_index(["repeat", "fold"])
    d_r2, d_mae = fa.r2 - fb.r2, fb.mae - fa.mae
    out = dict(set=a, reference=b)
    for k, d in (("d_r2", d_r2), ("d_mae", d_mae)):
        t = stats.corrected_ttest(d, fa.n_train, fa.n_test, alpha=alpha)
        out.update({f"{k}_mean": t["mean"], f"{k}_ci_low": t["ci_low"],
                    f"{k}_ci_high": t["ci_high"], f"{k}_p": t["p"]})
        out[f"{k}_n_folds_positive"] = int((d > 0).sum())
    out["d_mae_median"] = float(d_mae.median())
    return out


def main(run: Path, out: Path, n_jobs: int = 14) -> None:
    out.mkdir(parents=True, exist_ok=True)
    F = pd.read_csv(run / "features.csv", index_col="id")
    Y = pd.read_csv(run / "targets.csv", index_col="id")
    M = pd.read_csv(run / "molecules.csv", index_col="id")
    archived = pd.read_csv(run / "summary.csv")
    scaffold_run = run.parent / "main_scaffold" / "summary.csv"
    archived_scaffold = pd.read_csv(scaffold_run) if scaffold_run.exists() else None
    ok = M[M.status == "ok"]
    conn_all = ok.connectivity_group.fillna(pd.Series(ok.index, index=ok.index))
    scaffold_of = pd.Series(ev.scaffold_groups([Chem.MolFromSmiles(s) for s in ok.smiles]),
                            index=ok.index)
    summ, pair_rows, fold_rows, pred_rows = [], [], [], []
    for split in ("random", "scaffold"):
        for t in TARGETS:
            ids = Y.index[Y[t].notna() & (M.loc[Y.index, "status"] == "ok")]
            y = Y.loc[ids, t].astype(float)
            sets = {k: v.loc[ids] for k, v in feature_sets(F).items()}
            # exactly as experiments.cv_groups: labels are computed over all accepted
            # molecules (so acyclic singleton names match), then restricted to the rows
            # with this target; connectivity-group members share their first member's label
            groups = scaffold_of[conn_all.loc[ids]].to_numpy(dtype=object) \
                if split == "scaffold" else conn_all.loc[ids].to_numpy(dtype=object)
            t0 = time.time()
            res = {k: ev.nested_cv(X, y, FAMILIES, groups=groups, n_jobs=n_jobs)
                   for k, X in sets.items()}
            print(f"{split} {t} n={len(y)} {time.time() - t0:.0f}s", flush=True)
            for k, r in res.items():
                s = r.summary()
                summ.append(dict(split=split, target=t, set=k, n=len(y), p=r.meta["p"],
                                 r2_fold_mean=s["r2_fold_mean"], r2_fold_median=s["r2_fold_median"],
                                 r2_pooled_mean=s["r2_pooled_mean"], mae_fold_mean=s["mae_fold_mean"],
                                 rmse_fold_mean=s["rmse_fold_mean"], min_fold_r2=float(r.folds.r2.min()),
                                 family_counts=json.dumps(s["family_counts"]), reason=s["reason"]))
                fold_rows.append(r.folds.assign(split=split, target=t, set=k))
                pred_rows.append(r.predictions.assign(split=split, target=t, set=k))
            for a, b in [("TI_composition", "composition"), ("TI_asinh_composition", "composition"),
                         ("TI", "size"), ("TI_asinh", "size"), ("composition", "size"),
                         ("TI_asinh", "TI"), ("TI_asinh_composition", "TI_composition"),
                         ("rdkit2d_asinh", "rdkit2d"), ("rdkit2d_asinh", "composition")]:
                pair_rows.append(dict(split=split, target=t, **paired(res, a, b)))
            # the raw sets must reproduce the archived runs exactly: all raw sets of the
            # random split (main_confirmatory), and the registry-independent sets of the
            # scaffold split (main_scaffold used the 104-index registry for TI)
            ref, check = (archived, ("size", "composition", "TI", "TI_composition", "rdkit2d")) \
                if split == "random" else (archived_scaffold, ("size", "composition", "rdkit2d"))
            if ref is not None:
                for k in check:
                    a = ref[(ref.target == t) & (ref.learner == "selected")
                            & (ref.set == k)].r2_fold_mean.item()
                    mine = [r for r in summ if r["split"] == split and r["target"] == t
                            and r["set"] == k][0]["r2_fold_mean"]
                    assert np.isclose(a, mine, rtol=1e-9, atol=1e-9), (split, t, k, a, mine)
            pd.DataFrame(summ).to_csv(out / "summary.csv", index=False)
            pd.DataFrame(pair_rows).to_csv(out / "paired_tests.csv", index=False)
    pd.concat(fold_rows).to_csv(out / "folds.csv", index=False)
    pd.concat(pred_rows).to_csv(out / "predictions.csv.gz", index=False)
    print("reproduction of archived raw sets (random: all; scaffold: size/composition/rdkit2d): OK")


if __name__ == "__main__":
    main(Path(sys.argv[1]).expanduser(), Path(sys.argv[2]).expanduser())
