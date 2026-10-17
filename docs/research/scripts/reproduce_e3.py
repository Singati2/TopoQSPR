"""Notebook 01: recompute the committed E3 and P0 tables from the public E3 artifacts.

Input is a *restored* E3 directory (docs/research/data/E3 after restore_opera_values.py). Nothing is
trained. Independent code (numpy/scipy only) recomputes:
  * fold R2/MAE from per-molecule predictions, against the archived folds.csv;
  * the per-cell summary (fold-mean R2/MAE, pooled R2, family counts, pre-declared
    within-formula R2), against docs/research/out/E3/<cell>/summary.csv;
  * every paired contrast with the Nadeau-Bengio corrected resampled t, against paired.csv;
  * the external test (R2/MAE per set; paired contrasts with the seed-2026 formula-cluster
    bootstrap used by E3), against out/E3/external/{summary,paired}.csv;
  * the BH family (internal secondary contrasts + external, seed-99 bootstrap) against
    out/E3/recomputed_*.csv;
  * the post hoc centred within-formula metric, against out/E3/within_formula_posthoc_*.csv;
  * P0 (full vs edge types), by re-running p0_full_vs_edgetypes.py, against out/P0/.
Every comparison is a PASS/FAIL with an explicit tolerance. Not recomputable from artifacts
(they need re-training) and therefore only listed: the negative control and S-noflt.

    python -I docs/research/scripts/reproduce_e3.py RESTORED_DIR OUT_DIR
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
import repro_common as RC  # noqa: E402

COMMITTED = RC.REPO / "docs/research/out"
TOL = 1e-9          # same arithmetic on the same numbers: only float noise allowed
TOL_CSV = 1e-12     # relative slack for values that went through CSV round trips is not needed
                    # (float_precision="round_trip"), so TOL is used everywhere


def r2(y, p):
    s = float(((y - y.mean()) ** 2).sum())
    return np.nan if s == 0 else 1.0 - float(((y - p) ** 2).sum()) / s


def corrected_t(d, n1, n2, alpha=0.05):
    d = np.asarray(d, float); J = d.size; s2 = d.var(ddof=1)
    if s2 == 0:
        return d.mean(), np.nan, np.nan, np.nan
    se = np.sqrt((1 / J + np.mean(n2) / np.mean(n1)) * s2)
    c = st.t.ppf(1 - alpha / 2, J - 1)
    return d.mean(), d.mean() - c * se, d.mean() + c * se, 2 * st.t.sf(abs(d.mean() / se), J - 1)


def read(p: Path, **kw) -> pd.DataFrame:
    return pd.read_csv(p, float_precision="round_trip", **kw)


def predictions(cdir: Path) -> pd.DataFrame:
    f = cdir / "predictions.csv.gz"
    if not f.exists():
        raise SystemExit(f"{f} missing: run restore_opera_values.py first (OPERA cells need it)")
    return read(f)


def fold_metrics(P: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (s, r, k), d in P.groupby(["set", "repeat", "fold"], sort=True):
        y, p = d.y_true.to_numpy(), d.y_pred.to_numpy()
        rows.append(dict(set=s, repeat=r, fold=k, r2=r2(y, p), mae=float(np.abs(y - p).mean()), n_test=len(d)))
    return pd.DataFrame(rows)


def within_formula_predeclared(P: pd.DataFrame, formula: pd.Series) -> float:
    """06 §5 as implemented in e3_isomer.within_formula_r2 (independent re-implementation)."""
    vals = []
    for _, d in P.groupby("repeat"):
        d = d.assign(f=formula.loc[d.row].values)
        d = d[d.groupby(["fold", "f"]).row.transform("size") >= 2]
        den = ((d.y_true - d.groupby(["fold", "f"]).y_true.transform("mean")) ** 2).sum()
        vals.append(1 - ((d.y_true - d.y_pred) ** 2).sum() / den if den > 0 else np.nan)
    return float(np.mean(vals))


def within_formula_centred(P: pd.DataFrame, formula: pd.Series) -> pd.DataFrame:
    """POST HOC metric (07 §5): centre y and y_hat within (repeat, fold, formula)."""
    P = P.assign(f=formula.loc[P.row].values)
    P = P[P.groupby(["set", "repeat", "fold", "f"]).row.transform("size") >= 2]
    rows = []
    for s, d in P.groupby("set"):
        g = d.groupby(["repeat", "fold", "f"])
        yc = d.y_true - g.y_true.transform("mean"); pc = d.y_pred - g.y_pred.transform("mean")
        rows.append(dict(set=s, n=len(d), within_R2_centred=1 - ((yc - pc) ** 2).sum() / (yc ** 2).sum(),
                         within_corr=float(np.corrcoef(yc, pc)[0, 1]),
                         within_MAE_centred=float((yc - pc).abs().mean())))
    return pd.DataFrame(rows)


def internal(ck: RC.Checks, run: Path) -> pd.DataFrame:
    rows = []
    for cdir in sorted((run / "cells").glob("*_E3*")):
        cell = cdir.name
        endpoint = cell.rsplit("_", 1)[0]
        P = predictions(cdir)
        F = read(cdir / "folds.csv")
        meta = read(run / f"inclusion_{endpoint}.csv", index_col=0)
        S_c = read(COMMITTED / "E3" / cell / "summary.csv")
        A_c = read(COMMITTED / "E3" / cell / "paired.csv")
        # integrity
        RC.check_ids(ck, f"{cell} inclusion", meta.index)
        once = P.groupby(["set", "repeat", "row"]).size()
        ck.require(f"{cell}: each molecule predicted exactly once per set x repeat", bool((once == 1).all()),
                   int((once != 1).sum()), 0)
        ck.require(f"{cell}: every included molecule predicted", set(P.row) == set(meta.index),
                   len(set(meta.index) ^ set(P.row)), 0)
        same = P.pivot_table(index=["repeat", "row"], columns="set", values="fold").nunique(axis=1)
        ck.require(f"{cell}: identical fold membership across representations", bool((same == 1).all()),
                   int((same != 1).sum()), 0)
        gcol = "constitution" if cell.endswith("E3a") else "formula"
        g = meta[gcol]
        straddle = P.drop_duplicates(["repeat", "row"]).assign(g=g.loc[P.drop_duplicates(["repeat", "row"]).row].values) \
                    .groupby(["repeat", "g"]).fold.nunique().gt(1).sum()
        ck.require(f"{cell}: no {gcol} group in more than one outer fold (train/test disjoint)", straddle == 0, straddle, 0)
        ck.require(f"{cell}: predictions finite", bool(np.isfinite(P.y_pred).all() and np.isfinite(P.y_true).all()),
                   int((~np.isfinite(P.y_pred)).sum()), 0)
        # fold metrics
        fm = fold_metrics(P).merge(F[["set", "repeat", "fold", "n_train", "r2", "mae", "family"]],
                                   on=["set", "repeat", "fold"], suffixes=("", "_arch"))
        ck.require(f"{cell}: 15 folds per set", bool((fm.groupby("set").size() == 15).all()),
                   fm.groupby("set").size().to_dict(), 15)
        ck.close(f"{cell}: fold R2 from predictions = archived folds.csv (max |diff|)",
                 float(np.nanmax(np.abs(fm.r2 - fm.r2_arch))), 0.0, TOL)
        # summary
        for _, s in S_c.iterrows():
            f = fm[fm.set == s.set]
            pooled = np.mean([r2(d.y_true.to_numpy(), d.y_pred.to_numpy())
                              for _, d in P[P.set == s.set].groupby("repeat")])
            fam = json.dumps(f.family.value_counts().sort_index().to_dict())
            wf = within_formula_predeclared(P[P.set == s.set], meta.formula)
            ck.close(f"{cell} {s.set}: R2 fold mean", f.r2.mean(), s.r2_fold_mean, TOL)
            ck.close(f"{cell} {s.set}: R2 pooled", pooled, s.r2_pooled_mean, TOL)
            ck.close(f"{cell} {s.set}: MAE fold mean", f.mae.mean(), s.mae_fold_mean, TOL)
            ck.close(f"{cell} {s.set}: within-formula R2 (PRE-DECLARED metric)", wf, s.within_formula_r2, TOL)
            ck.require(f"{cell} {s.set}: selected families", json.loads(fam) == json.loads(s.family_counts),
                       fam, s.family_counts)
        # paired contrasts
        for _, a in A_c.iterrows():
            A = fm[fm.set == a.set].sort_values(["repeat", "fold"])
            B = fm[fm.set == a.reference].sort_values(["repeat", "fold"])
            m, lo, hi, p = corrected_t(A.r2.to_numpy() - B.r2.to_numpy(), A.n_train, A.n_test)
            m90, lo90, hi90, _ = corrected_t(A.r2.to_numpy() - B.r2.to_numpy(), A.n_train, A.n_test, 0.10)
            dm, dlo, dhi, _ = corrected_t(B.mae.to_numpy() - A.mae.to_numpy(), A.n_train, A.n_test)
            lab = f"{cell} {a.set} - {a.reference}"
            for nm, x, y in (("dR2", m, a.d_r2_mean), ("CI low", lo, a.d_r2_ci_low), ("CI high", hi, a.d_r2_ci_high),
                             ("p", p, a.d_r2_p), ("90% CI low", lo90, a.d_r2_ci90_low),
                             ("90% CI high", hi90, a.d_r2_ci90_high), ("dMAE", dm, a.d_mae_mean)):
                ck.close(f"{lab}: {nm}", x, y, TOL)
            rows.append(dict(cell=cell, set=a.set, reference=a.reference, d_r2=m, lo=lo, hi=hi, p=p,
                             ci90_low=lo90, ci90_high=hi90, d_mae=dm, mae_lo=dlo, mae_hi=dhi,
                             r2_set=A.r2.mean(), r2_reference=B.r2.mean(),
                             mae_set=A.mae.mean(), mae_reference=B.mae.mean()))
        # post hoc within-formula (MP only, as committed)
        cf = COMMITTED / "E3" / f"within_formula_posthoc_{cell}.csv"
        if cf.exists():
            W = within_formula_centred(P, meta.formula).set_index("set")
            Wc = read(cf).set_index("set")
            for s in Wc.index:
                ck.close(f"{cell} {s}: within-formula centred R2 (POST HOC metric)",
                         W.loc[s, "within_R2_centred"], Wc.loc[s, "within_R2_centred"], TOL)
    return pd.DataFrame(rows)


def cluster_boot(F: np.ndarray, rng, B=2000):
    uf = np.unique(F); mem = [np.flatnonzero(F == f) for f in uf]
    return [np.concatenate([mem[j] for j in rng.integers(0, uf.size, uf.size)]) for _ in range(B)]


def external(ck: RC.Checks, run: Path, seed: int) -> pd.DataFrame:
    """Recompute external R2/MAE and paired bootstrap CIs (seed 2026 = E3's own; 99 = recompute)."""
    pr = read(run / "external/predictions.csv.gz", index_col=0)
    split = read(run / "external/test_split.csv").set_index("id")
    RC.check_ids(ck, f"external (seed {seed}) predictions", pr.index, split.index)
    ck.require(f"external (seed {seed}): split order = prediction order", list(split.index) == list(pr.index),
               "differs", "identical")
    ck.require(f"external (seed {seed}): y finite after restore", bool(np.isfinite(pr.y).all()),
               int((~np.isfinite(pr.y)).sum()), 0)
    S_c = read(COMMITTED / "E3/external/summary.csv")
    A_c = read(COMMITTED / "E3/external/paired.csv")
    rng = np.random.default_rng(seed)
    rows = []
    for test in ("known_formula", "unseen_formula"):  # E3's order: the rng stream depends on it
        d = pr[split.test == test]; yy = d.y.to_numpy(); Fm = d.formula.to_numpy()
        boots = cluster_boot(Fm, rng)
        if seed == 2026:
            for _, s in S_c[S_c.test == test].iterrows():
                ck.close(f"external {test} {s.set}: R2", r2(yy, d[s.set].to_numpy()), s.r2, TOL)
                ck.close(f"external {test} {s.set}: MAE", float(np.abs(yy - d[s.set]).mean()), s.mae, TOL)
        for _, a in A_c[A_c.test == test].iterrows():
            pa, pb = d[a.set].to_numpy(), d[a.reference].to_numpy()
            dr = r2(yy, pa) - r2(yy, pb)
            bs = np.array([r2(yy[i], pa[i]) - r2(yy[i], pb[i]) for i in boots])
            row = dict(test=test, n=len(d), n_formulas=int(np.unique(Fm).size), set=a.set, reference=a.reference,
                       d_r2=dr, lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                       ci90_low=np.percentile(bs, 5), ci90_high=np.percentile(bs, 95), boot_se=float(bs.std(ddof=1)),
                       p_boot=float(min(1, 2 * min((bs <= 0).mean(), (bs >= 0).mean()))))
            if seed == 2026:
                lab = f"external {test} {a.set} - {a.reference}"
                ck.close(f"{lab}: dR2", dr, a.d_r2, TOL)
                ck.close(f"{lab}: CI low (seed 2026)", row["lo"], a.d_r2_ci_low, TOL)
                ck.close(f"{lab}: CI high (seed 2026)", row["hi"], a.d_r2_ci_high, TOL)
                ck.close(f"{lab}: bootstrap p (seed 2026)", row["p_boot"], a.boot_p_two_sided, TOL)
            rows.append(row)
    return pd.DataFrame(rows)


def bh_family(ck: RC.Checks, R: pd.DataFrame, E99: pd.DataFrame) -> None:
    """BH over internal secondary + external (seed 99) contrasts, as e3_recompute_and_failures.py."""
    Rc = read(COMMITTED / "E3/recomputed_paired.csv")
    Ec = read(COMMITTED / "E3/recomputed_external.csv")
    R = R.merge(Rc[["cell", "set", "reference", "role", "q_bh_family"]], on=["cell", "set", "reference"], how="left")
    ck.require("BH family: every committed internal contrast recomputed", R.q_bh_family.notna().sum() + (R.role == "primary").sum() == len(Rc),
               len(R), len(Rc))
    # committed row order (cell sorted, CONTRASTS order) is the order p-values entered BH
    R = Rc[["cell", "set", "reference"]].merge(R, on=["cell", "set", "reference"])
    E = Ec[["test", "set", "reference"]].merge(E99, on=["test", "set", "reference"])
    sec = R.role == "secondary"
    q = RC.bh(np.r_[R.loc[sec, "p"].to_numpy(), E.p_boot.to_numpy()])
    ck.close("BH family size (internal secondary + external)", float(q.size), float(sec.sum() + len(Ec)), 0)
    ck.close("BH q: max |diff| vs committed, internal",
             float(np.nanmax(np.abs(q[: sec.sum()] - R.loc[sec, "q_bh_family"].to_numpy()))), 0.0, TOL)
    ck.close("BH q: max |diff| vs committed, external (seed 99)",
             float(np.nanmax(np.abs(q[sec.sum():] - Ec.q_bh_family.to_numpy()))), 0.0, TOL)
    for _, (c, s, r_) in enumerate([("mp_bradley_E3b", "ECFP_TIa", "ECFP"), ("logp_opera_E3a", "TIa_composition", "composition")]):
        i = np.flatnonzero(((R.cell == c) & (R.set == s) & (R.reference == r_))[sec].to_numpy())
        if i.size:
            ck.info(f"BH q {c} {s} - {r_} (reported in 07)", round(float(q[i[0]]), 3))


def p0(ck: RC.Checks, run: Path, out: Path) -> None:
    import p0_full_vs_edgetypes as P0
    P0.main(run, out / "P0")
    for name in ("internal_full_vs_edgetypes.csv", "external_full_vs_edgetypes.csv",
                 "learner_families_full_vs_edgetypes.csv"):
        a = read(COMMITTED / "P0" / name); b = read(out / "P0" / name)
        ck.require(f"P0 {name}: same shape and columns", a.shape == b.shape and list(a.columns) == list(b.columns),
                   b.shape, a.shape)
        num = a.select_dtypes("number").columns
        if len(num):
            ck.close(f"P0 {name}: max |diff| vs committed", float(np.nanmax(np.abs(a[num].to_numpy(float) - b[num].to_numpy(float)))), 0.0, TOL)
        txt = a.select_dtypes(exclude="number").columns
        ck.require(f"P0 {name}: text columns identical", a[txt].astype(str).equals(b[txt].astype(str)), "", "")


def main(run: Path, out: Path) -> RC.Checks:
    out.mkdir(parents=True, exist_ok=True)
    ck = RC.Checks("01_E3_P0_reproduction")
    for k, v in RC.version_mismatches().items():
        ck.info(f"version differs from the pinned set: {k}", f"pinned {v[0]}, installed {v[1]}")
    R = internal(ck, run)
    E26 = external(ck, run, 2026)
    E99 = external(ck, run, 99)
    bh_family(ck, R, E99)
    p0(ck, run, out)
    ck.info("negative control, S-noflt", "not recomputable from artifacts (need re-training); committed values only",
            "docs/research/out/E3/negative_control_paired.csv, sensitivity_noflt_paired.csv")
    RC.atomic_write_csv(R, out / "recomputed_internal_contrasts.csv", index=False)
    RC.atomic_write_csv(pd.concat([E26.assign(seed=2026), E99.assign(seed=99)]), out / "recomputed_external_contrasts.csv", index=False)
    ck.finish(out)
    return ck


if __name__ == "__main__":
    main(Path(sys.argv[1]).expanduser(), Path(sys.argv[2]).expanduser())
