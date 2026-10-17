"""P0 step 1 (doc 12): paired contrast TIa_composition - edgetypes_composition from archived E3 outputs.

No fitting. Uses the archived fold scores and predictions only:
- internal cells: corrected resampled t (Nadeau-Bengio, as in E3) over the 15 seed x fold scores,
  recomputed from predictions and checked against the archived fold R2;
- external test: formula-cluster bootstrap (2000 resamples) on the archived predictions, with the
  run's seed (2026) and an independent seed (99) to show Monte Carlo stability.
Also tabulates the selected learner family per fold for both sets. The contrast compares the
best-selected model per representation, so it does NOT remove learner confounding.

    python -I docs/research/scripts/p0_full_vs_edgetypes.py ~/topoqspr_runs/research_2026-10-13_E3
    python -I docs/research/scripts/p0_full_vs_edgetypes.py <restored dataset dir> <output dir>   # Colab

The second form runs on the public dataset after restore_opera_values.py has put the OPERA values back.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

A_SET, B_SET = "TIa_composition", "edgetypes_composition"


def r2(y, p):
    s = ((y - y.mean()) ** 2).sum()
    return np.nan if s == 0 else 1 - ((y - p) ** 2).sum() / s


def corrected_t(d, n1, n2, alpha=0.05):
    d = np.asarray(d, float); J = d.size; s2 = d.var(ddof=1)
    se = np.sqrt((1 / J + np.mean(n2) / np.mean(n1)) * s2); c = st.t.ppf(1 - alpha / 2, J - 1)
    return d.mean(), d.mean() - c * se, d.mean() + c * se, 2 * st.t.sf(abs(d.mean() / se), J - 1)


def internal(run: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, fams = [], []
    for cdir in sorted((run / "cells").glob("*_E3*")):
        P = pd.read_csv(cdir / "predictions.csv.gz", float_precision="round_trip")
        F = pd.read_csv(cdir / "folds.csv", float_precision="round_trip")
        P = P[P.set.isin([A_SET, B_SET])]
        fm = P.groupby(["set", "repeat", "fold"]).apply(
            lambda d: pd.Series(dict(r2=r2(d.y_true.to_numpy(), d.y_pred.to_numpy()),
                                     mae=float(np.abs(d.y_true - d.y_pred).mean()))),
            include_groups=False).reset_index()
        fm = fm.merge(F[["set", "repeat", "fold", "n_train", "n_test", "family", "r2"]],
                      on=["set", "repeat", "fold"], suffixes=("", "_arch"))
        a = fm[fm.set == A_SET].sort_values(["repeat", "fold"]).reset_index(drop=True)
        b = fm[fm.set == B_SET].sort_values(["repeat", "fold"]).reset_index(drop=True)
        d = a.r2.to_numpy() - b.r2.to_numpy()
        m, lo, hi, p = corrected_t(d, a.n_train, a.n_test)
        _, lo90, hi90, _ = corrected_t(d, a.n_train, a.n_test, alpha=0.10)
        dm, mlo, mhi, _ = corrected_t(b.mae.to_numpy() - a.mae.to_numpy(), a.n_train, a.n_test)
        same = (a.family.to_numpy() == b.family.to_numpy())
        rows.append(dict(context=cdir.name, n_folds=d.size, r2_full=a.r2.mean(), r2_edge=b.r2.mean(),
                         d_r2=m, ci_low=lo, ci_high=hi, ci90_low=lo90, ci90_high=hi90, p=p,
                         folds_full_better=int((d > 0).sum()), d_mae=dm, mae_ci_low=mlo, mae_ci_high=mhi,
                         folds_same_family=int(same.sum()),
                         d_r2_same_family_folds=float(d[same].mean()) if same.any() else np.nan,
                         max_fold_r2_diff_vs_archived=float(np.nanmax(np.abs(
                             np.r_[a.r2 - a.r2_arch, b.r2 - b.r2_arch])))))
        for s, x in ((A_SET, a), (B_SET, b)):
            for fam, n in x.family.value_counts().items():
                fams.append(dict(context=cdir.name, set=s, family=fam, folds=int(n)))
    return pd.DataFrame(rows), pd.DataFrame(fams)


def external(run: Path, out: Path) -> pd.DataFrame:
    pr = pd.read_csv(run / "external" / "predictions.csv.gz", index_col=0, float_precision="round_trip")
    summ = pd.read_csv(run / "external" / "summary.csv")
    split_file = run / "external" / "test_split.csv"
    if split_file.exists():  # published dataset: split precomputed, no training-set rebuild needed
        split = pd.read_csv(split_file, index_col=0)
        known = split.loc[pr.index, "test"].eq("known_formula")
    else:  # local archive: rebuild the training formulas (data building only, no fitting)
        sys.path.insert(0, str(Path(__file__).parent))
        import e3_isomer as E3
        _, tmeta_train, _, _ = E3.build("mp_bradley", restrict=E3.physprop_free_bradley_keys())
        known = pr.formula.isin(set(tmeta_train.formula))
        pd.DataFrame(dict(formula=pr.formula, test=np.where(known, "known_formula", "unseen_formula")),
                     index=pr.index).to_csv(out / "external_split.csv")
    sizes = json.loads((run / "external" / "info.json").read_text())["test_sizes"]
    assert int(known.sum()) == sizes["known_formula"] and int((~known).sum()) == sizes["unseen_formula"], \
        "known/unseen split does not reproduce the archived test sizes"
    rows = []
    for test, mask in (("known_formula", known), ("unseen_formula", ~known)):
        d = pr[mask]; yy = d.y.to_numpy(); Fm = d.formula.to_numpy()
        pa, pb = d[A_SET].to_numpy(), d[B_SET].to_numpy()
        uf = np.unique(Fm); mem = [np.flatnonzero(Fm == f) for f in uf]
        dr = r2(yy, pa) - r2(yy, pb)
        dmae = np.abs(yy - pb).mean() - np.abs(yy - pa).mean()
        fam = summ[summ.test == test].set_index("set").family
        row = dict(context=f"external_{test}", n=len(d), n_formulas=uf.size, r2_full=r2(yy, pa),
                   r2_edge=r2(yy, pb), d_r2=dr, d_mae=dmae, family_full=fam[A_SET], family_edge=fam[B_SET])
        for seed in (2026, 99):
            rng = np.random.default_rng(seed)
            boots = [np.concatenate([mem[j] for j in rng.integers(0, uf.size, uf.size)]) for _ in range(2000)]
            bs = np.array([r2(yy[i], pa[i]) - r2(yy[i], pb[i]) for i in boots])
            bm = np.array([np.abs(yy[i] - pb[i]).mean() - np.abs(yy[i] - pa[i]).mean() for i in boots])
            row.update({f"ci_low_s{seed}": np.percentile(bs, 2.5), f"ci_high_s{seed}": np.percentile(bs, 97.5),
                        f"ci90_low_s{seed}": np.percentile(bs, 5), f"ci90_high_s{seed}": np.percentile(bs, 95),
                        f"p_boot_s{seed}": float(min(1, 2 * min((bs <= 0).mean(), (bs >= 0).mean()))),
                        f"mae_ci_low_s{seed}": np.percentile(bm, 2.5), f"mae_ci_high_s{seed}": np.percentile(bm, 97.5)})
        rows.append(row)
    return pd.DataFrame(rows)


def main(run: Path, out: Path | None = None) -> None:
    out = out or Path(__file__).resolve().parents[1] / "out" / "P0"
    out.mkdir(parents=True, exist_ok=True)
    I, fams = internal(run)
    E = external(run, out)
    I.to_csv(out / "internal_full_vs_edgetypes.csv", index=False)
    fams.to_csv(out / "learner_families_full_vs_edgetypes.csv", index=False)
    E.to_csv(out / "external_full_vs_edgetypes.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print(I.round(4).to_string(index=False)); print(); print(E.round(4).T.to_string()); print()
    print(fams.pivot_table(index=["context", "set"], columns="family", values="folds", fill_value=0))


if __name__ == "__main__":
    main(Path(sys.argv[1]).expanduser(), Path(sys.argv[2]).expanduser() if len(sys.argv) > 2 else None)
