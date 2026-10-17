"""Notebook 04: statistical robustness of the E3 inference (exploratory, post hoc checks).

  (a) bootstrap unit: formula-cluster vs molecule-level bootstrap on the archived external
      predictions (training held fixed in both: test-set sampling uncertainty only);
  (b) coverage of the corrected resampled t (Nadeau & Bengio 2003) under formula-grouped
      repeated CV, by simulation on the MP isomer data (finite-population design below);
  (c) training-seed variability of the archived external pipeline (E3 select-and-refit over the
      six families): seed 0 must reproduce the archived predictions exactly; seeds 1-4 measure
      training variability (FULL mode only, needs owner approval like doc 14);
  (d) within-formula metrics: the PRE-DECLARED within-formula R2 (06 §5) next to the POST HOC
      centred variant (07 §5), for every cell, both labelled.

Coverage design (b). The 758 MP formulas are the population. A replicate draws half of the
formulas (the "study sample"), runs E3b-type grouped 3 x 5 CV of C+T vs C with a fixed learner
(ridge, alpha 1; SVR C=3 in FULL mode), and forms the corrected-t 95% CI of the fold-mean
dR2, plus the naive t CI (no correction) for comparison. The target is the expected fold-level
dR2 for the same train/test sizes, estimated by Monte Carlo over random formula-disjoint
train/test draws from the whole population. Study samples overlap between replicates (finite
population), so this is an approximation of coverage, not an exact frequentist guarantee.

    python -I docs/research/scripts/robustness.py RESTORED_DIR FIXED_LEARNER_OUT OUT --mode smoke
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
import repro_common as RC  # noqa: E402
import reproduce_e3 as R3  # noqa: E402

COMMITTED = RC.REPO / "docs/research/out/E3"


# ----------------------------------------------------------------------------- (a)
def bootstrap_unit(restored: Path, B: int = 2000) -> pd.DataFrame:
    pr = pd.read_csv(restored / "external/predictions.csv.gz", index_col=0, float_precision="round_trip")
    split = pd.read_csv(restored / "external/test_split.csv").set_index("id")
    rows = []
    for test in ("known_formula", "unseen_formula"):
        d = pr[split.test == test]; yy = d.y.to_numpy(); Fm = d.formula.to_numpy()
        n_per = pd.Series(Fm).value_counts()
        rng_c, rng_m = np.random.default_rng(2026), np.random.default_rng(2027)
        bc = R3.cluster_boot(Fm, rng_c, B)
        bm = [rng_m.integers(0, len(yy), len(yy)) for _ in range(B)]
        for a, b in (("TIa_composition", "composition"), ("TIa_composition", "edgetypes_composition"),
                     ("rdkit2d_asinh_TIa", "rdkit2d_asinh"), ("ECFP_TIa", "ECFP")):
            pa, pb = d[a].to_numpy(), d[b].to_numpy()
            f = lambda i: R3.r2(yy[i], pa[i]) - R3.r2(yy[i], pb[i])
            sc = np.array([f(i) for i in bc]); sm = np.array([f(i) for i in bm])
            rows.append(dict(test=test, set=a, reference=b, n=len(yy), n_formulas=int(n_per.size),
                             mean_molecules_per_formula=round(float(n_per.mean()), 2),
                             d_r2=f(np.arange(len(yy))), se_formula_cluster=sc.std(ddof=1), se_molecule=sm.std(ddof=1),
                             se_ratio_cluster_over_molecule=sc.std(ddof=1) / sm.std(ddof=1),
                             ci95_cluster=f"[{np.percentile(sc, 2.5):.4f}, {np.percentile(sc, 97.5):.4f}]",
                             ci95_molecule=f"[{np.percentile(sm, 2.5):.4f}, {np.percentile(sm, 97.5):.4f}]",
                             note="training fixed (single E3 fit): test-sampling uncertainty only"))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- (b)
def _pipe(learner: str):
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import Ridge
    from sklearn.svm import SVR
    from sklearn.compose import TransformedTargetRegressor
    from sklearn.impute import SimpleImputer
    imp = SimpleImputer(strategy="median", keep_empty_features=True)  # fitted on training rows only
    if learner == "ridge":
        return make_pipeline(imp, StandardScaler(), Ridge(alpha=1.0))
    return TransformedTargetRegressor(make_pipeline(imp, StandardScaler(), SVR(C=3.0, epsilon=0.1, gamma="scale")),
                                      transformer=StandardScaler())


def _delta(XA, XB, y, tr, te, learner):
    from sklearn.base import clone
    pa = clone(_pipe(learner)).fit(XA[tr], y[tr]).predict(XA[te])
    pb = clone(_pipe(learner)).fit(XB[tr], y[tr]).predict(XB[te])
    return R3.r2(y[te], pa) - R3.r2(y[te], pb)


def coverage(fl_out: Path, learner: str, n_rep: int, n_truth: int, n_jobs: int, seed: int = 4242) -> dict:
    """Corrected-t vs naive-t coverage under formula-grouped repeated CV (see module doc)."""
    from joblib import Parallel, delayed
    from topoqspr.evaluate import balanced_group_kfold
    with np.load(fl_out / "data/internal.npz", allow_pickle=False) as z:
        I = dict(z)
    XA, XB, y, form = I["C+T"], I["C"], I["y"], I["formula"]
    uf = np.unique(form)
    rng = np.random.default_rng(seed)
    samples = [rng.choice(uf, uf.size // 2, replace=False) for _ in range(n_rep)]
    n_s = int(np.mean([np.isin(form, s).sum() for s in samples]))
    n1, n2 = int(round(n_s * 0.8)), int(round(n_s * 0.2))

    def one_rep(s):
        idx = np.flatnonzero(np.isin(form, s))
        d, ntr, nte = [], [], []
        for r in range(3):
            for tr, te in balanced_group_kfold(form[idx], 5, r):
                d.append(_delta(XA, XB, y, idx[tr], idx[te], learner)); ntr.append(tr.size); nte.append(te.size)
        d = np.array(d)
        m, lo, hi, _ = R3.corrected_t(d, ntr, nte)
        se_n = d.std(ddof=1) / np.sqrt(d.size); c = st.t.ppf(0.975, d.size - 1)
        return m, lo, hi, m - c * se_n, m + c * se_n

    def one_truth(j):
        r = np.random.default_rng(seed + 1 + j)
        perm = r.permutation(uf)
        cnt = pd.Series(form).value_counts()
        cum = np.cumsum(cnt.loc[perm].to_numpy())
        te_f = perm[: np.searchsorted(cum, n2) + 1]
        rest = perm[np.searchsorted(cum, n2) + 1:]
        cum2 = np.cumsum(cnt.loc[rest].to_numpy())
        tr_f = rest[: np.searchsorted(cum2, n1) + 1]
        return _delta(XA, XB, y, np.flatnonzero(np.isin(form, tr_f)), np.flatnonzero(np.isin(form, te_f)), learner)

    t0 = time.time()
    reps = np.array(Parallel(n_jobs=n_jobs)(delayed(one_rep)(s) for s in samples))
    truth = np.array(Parallel(n_jobs=n_jobs)(delayed(one_truth)(j) for j in range(n_truth)))
    theta = truth.mean(); theta_se = truth.std(ddof=1) / np.sqrt(truth.size)
    cov_c = float(np.mean((reps[:, 1] <= theta) & (theta <= reps[:, 2])))
    cov_n = float(np.mean((reps[:, 3] <= theta) & (theta <= reps[:, 4])))
    k = len(reps)
    return dict(learner=learner, n_replicates=k, n_truth_draws=n_truth, study_sample_rows_mean=n_s,
                train_rows=n1, test_rows=n2, theta_target=theta, theta_mc_se=theta_se,
                coverage_corrected_t=cov_c, coverage_corrected_t_mc_se=float(np.sqrt(cov_c * (1 - cov_c) / k)),
                coverage_naive_t=cov_n, mean_width_corrected=float(np.mean(reps[:, 2] - reps[:, 1])),
                mean_width_naive=float(np.mean(reps[:, 4] - reps[:, 3])),
                mean_estimate=float(reps[:, 0].mean()), sd_estimate=float(reps[:, 0].std(ddof=1)),
                seconds=round(time.time() - t0, 1),
                caveat="finite-population approximation: study samples overlap; target estimated by Monte Carlo")


# ----------------------------------------------------------------------------- (c)
def external_seeds(fl_out: Path, restored: Path, sets: list, seeds: list, n_jobs: int) -> pd.DataFrame:
    """Re-run the E3 external select-and-refit (all six families) with random_state = seed."""
    from joblib import Parallel, delayed
    import e3_isomer as E3
    from topoqspr import evaluate as ev
    with np.load(fl_out / "data/external.npz", allow_pickle=False) as z:
        X = dict(z)  # picklable for joblib workers
    name = {"composition": "C", "edgetypes_composition": "C+E", "TIa_composition": "C+T"}
    pr = pd.read_csv(restored / "external/predictions.csv.gz", index_col=0, float_precision="round_trip")
    split = pd.read_csv(restored / "external/test_split.csv").set_index("id")
    assert list(X["test_ids"]) == list(split.index) == list(pr.index)

    def fit(k, s):
        sel = ev.select_and_refit(X[f"train_{name[k]}"], X["y_train"], E3.FAMILIES, inner=3, random_state=s,
                                  groups=X["train_constitution"].astype(object))
        return k, s, sel.family, sel.estimator.predict(X[f"test_{name[k]}"])
    fits = Parallel(n_jobs=n_jobs)(delayed(fit)(k, s) for k in sets for s in seeds)
    rows = []
    for k, s, fam, p in fits:
        for test in ("known_formula", "unseen_formula"):
            m = (split.test == test).to_numpy(); yy = pr.y.to_numpy()[m]
            rows.append(dict(set=k, seed=s, family=fam, test=test, r2=R3.r2(yy, p[m]),
                             mae=float(np.abs(yy - p[m]).mean()),
                             max_abs_diff_vs_archived_prediction=float(np.abs(p - pr[k].to_numpy()).max()) if s == 0 else np.nan))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- (d)
def within_formula(restored: Path) -> pd.DataFrame:
    rows = []
    for cdir in sorted((restored / "cells").glob("*_E3*")):
        endpoint = cdir.name.rsplit("_", 1)[0]
        P = R3.predictions(cdir)
        meta = pd.read_csv(restored / f"inclusion_{endpoint}.csv", index_col=0)
        W = R3.within_formula_centred(P, meta.formula).set_index("set")
        for s in sorted(P.set.unique()):
            rows.append(dict(cell=cdir.name, set=s,
                             within_formula_r2_PREDECLARED=R3.within_formula_predeclared(P[P.set == s], meta.formula),
                             within_formula_r2_centred_POSTHOC=W.loc[s, "within_R2_centred"],
                             within_corr_centred_POSTHOC=W.loc[s, "within_corr"]))
    return pd.DataFrame(rows)


def main(restored: Path, fl_out: Path, out: Path, mode: str, n_jobs: int) -> RC.Checks:
    out.mkdir(parents=True, exist_ok=True)
    ck = RC.Checks(f"04_statistical_robustness_{mode}")
    A = bootstrap_unit(restored)
    RC.atomic_write_csv(A, out / "a_bootstrap_unit.csv", index=False)
    ck.require("(a) bootstrap SEs finite and positive", bool((A[["se_formula_cluster", "se_molecule"]] > 0).all().all()), "", "")
    ck.info("(a) formula-cluster / molecule-level bootstrap SE ratio per contrast (finding, not a gate)",
            "; ".join(f"{r.test} {r.set}-{r.reference}: {r.se_ratio_cluster_over_molecule:.2f}" for r in A.itertuples()))
    cov = [coverage(fl_out, "ridge", 40 if mode == "smoke" else 400, 100 if mode == "smoke" else 1000, n_jobs)]
    if mode == "full":
        cov.append(coverage(fl_out, "svr", 100, 300, n_jobs))
    Cv = pd.DataFrame(cov)
    if mode == "smoke":
        Cv.insert(0, "WARNING", "SMOKE: too few replicates for a coverage estimate")
    RC.atomic_write_csv(Cv, out / "b_coverage_corrected_t.csv", index=False)
    for r in cov:
        ck.info(f"(b) coverage {r['learner']}: corrected t {r['coverage_corrected_t']:.3f} "
                f"(MC SE {r['coverage_corrected_t_mc_se']:.3f}), naive t {r['coverage_naive_t']:.3f}",
                json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}))
    seeds = [0] if mode == "smoke" else [0, 1, 2, 3, 4]
    sets = ["composition"] if mode == "smoke" else ["composition", "edgetypes_composition", "TIa_composition"]
    Sx = external_seeds(fl_out, restored, sets, seeds, n_jobs)
    RC.atomic_write_csv(Sx, out / "c_external_training_seeds.csv", index=False)
    s0 = Sx[Sx.seed == 0]
    ck.close("(c) seed 0 re-fit reproduces archived external predictions (max |diff|, degC)",
             float(s0.max_abs_diff_vs_archived_prediction.max()), 0.0, 1e-6,
             note="cross-machine refit; BLAS/libsvm may differ in the last bits")
    W = within_formula(restored)
    RC.atomic_write_csv(W, out / "d_within_formula_both_metrics.csv", index=False)
    Sc = {c: pd.read_csv(COMMITTED / c / "summary.csv", float_precision="round_trip").set_index("set")
          for c in W.cell.unique()}
    ck.close("(d) pre-declared within-formula R2 matches committed summaries (max |diff|)",
             float(max(abs(r.within_formula_r2_PREDECLARED - Sc[r.cell].loc[r.set, "within_formula_r2"]) for r in W.itertuples())), 0.0, 1e-9)
    ck.finish(out)
    return ck


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("restored"); ap.add_argument("fl_out"); ap.add_argument("out")
    ap.add_argument("--mode", choices=["smoke", "full"], default="smoke")
    ap.add_argument("--n-jobs", type=int, default=1)
    a = ap.parse_args()
    if a.mode == "full" and not (Path(a.fl_out) / "APPROVED_FULL_RUN").exists():
        raise SystemExit("full mode (training seeds 1-4) needs the same owner approval file as fixed_learner.py")
    main(*(Path(p).expanduser() for p in (a.restored, a.fl_out, a.out)), a.mode, a.n_jobs)
