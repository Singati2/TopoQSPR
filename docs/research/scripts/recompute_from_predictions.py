"""Independent recomputation of the frozen Experiment B statistics from the archived
out-of-fold predictions (mission §6). Uses no topoqspr code: only numpy/scipy/pandas.

For each archived run directory:
  * fold integrity: each molecule exactly once per (target, learner, set, repeat); identical
    fold membership across sets and learners within a target and repeat; distinct outer seeds;
  * fold metrics: R^2 = 1 - SSE/SST with SST around the *test-fold* mean (sklearn r2_score
    convention), MAE, RMSE, compared with folds.csv;
  * pooled R^2 per repeat, compared with summary.csv r2_pooled_mean (mean over repeats);
  * paired TI - size differences, Nadeau-Bengio corrected repeated-CV t (independent code),
    BH and BY over the 4-endpoint family, compared with confirmatory.csv.

    python -I docs/research/scripts/recompute_from_predictions.py \
        ~/topoqspr_runs/main_confirmatory ~/topoqspr_runs/main_confirmatory_protocol_svr
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st


def r2(y, p):
    sst = float(((y - y.mean()) ** 2).sum())
    return np.nan if sst == 0 else 1.0 - float(((y - p) ** 2).sum()) / sst


def corrected_t(d, n1, n2, alpha=0.05):
    d = np.asarray(d, float)
    J = d.size
    s2 = d.var(ddof=1)
    se = np.sqrt((1 / J + np.mean(n2) / np.mean(n1)) * s2)
    t = d.mean() / se
    c = st.t.ppf(1 - alpha / 2, J - 1)
    return d.mean(), d.mean() - c * se, d.mean() + c * se, 2 * st.t.sf(abs(t), J - 1)


def bh(p):
    p = np.asarray(p, float); m = p.size; o = np.argsort(p); q = np.empty(m); run = 1.0
    for j in range(m - 1, -1, -1):
        run = min(run, p[o[j]] * m / (j + 1)); q[o[j]] = run
    return q


def by(p):
    m = len(p)
    return np.minimum(1.0, bh(p) * sum(1 / i for i in range(1, m + 1)))


def check(run: Path) -> dict:
    P = pd.read_csv(run / "predictions.csv", float_precision="round_trip")
    F = pd.read_csv(run / "folds.csv", float_precision="round_trip")
    S = pd.read_csv(run / "summary.csv", float_precision="round_trip")
    C = pd.read_csv(run / "confirmatory.csv", float_precision="round_trip")
    P["set"] = P["set"].fillna("null"); F["set"] = F["set"].fillna("null"); S["set"] = S["set"].fillna("null")
    out = dict(run=run.name)
    # --- integrity
    dup = P.groupby(["target", "learner", "set", "repeat", "row"]).size()
    out["rows_seen_more_than_once_per_repeat"] = int((dup > 1).sum())
    memb = P.groupby(["target", "repeat", "row"]).fold.nunique()
    out["rows_with_inconsistent_fold_across_sets"] = int((memb > 1).sum())
    out["distinct_outer_seeds"] = sorted(F.seed.unique().tolist())
    # --- fold metrics
    g = P.groupby(["target", "learner", "set", "repeat", "fold"])
    fm = g.apply(lambda d: pd.Series(dict(
        r2=r2(d.y_true.to_numpy(), d.y_pred.to_numpy()),
        mae=float(np.abs(d.y_true - d.y_pred).mean()),
        rmse=float(np.sqrt(((d.y_true - d.y_pred) ** 2).mean())), n_test=len(d))))
    fm = fm.reset_index().merge(F, on=["target", "learner", "set", "repeat", "fold"], suffixes=("", "_arch"))
    for k in ("r2", "mae", "rmse"):
        out[f"max_abs_diff_fold_{k}"] = float(np.nanmax(np.abs(fm[k] - fm[f"{k}_arch"])))
    out["fold_ntest_mismatch"] = int((fm.n_test != fm.n_test_arch).sum())
    # --- pooled R2 (mean over repeats)
    pr = P.groupby(["target", "learner", "set", "repeat"]).apply(
        lambda d: r2(d.y_true.to_numpy(), d.y_pred.to_numpy())).groupby(level=[0, 1, 2]).mean()
    pr = pr.rename("mine").reset_index().merge(S[["target", "learner", "set", "r2_pooled_mean"]],
                                               on=["target", "learner", "set"])
    out["max_abs_diff_pooled_r2"] = float(np.nanmax(np.abs(pr.mine - pr.r2_pooled_mean)))
    # --- confirmatory TI - size, selected learner
    rows = []
    for t in C.target:
        a = fm[(fm.target == t) & (fm.learner == "selected") & (fm.set == "TI")].sort_values(["repeat", "fold"])
        b = fm[(fm.target == t) & (fm.learner == "selected") & (fm.set == "size")].sort_values(["repeat", "fold"])
        mean, lo, hi, p = corrected_t(a.r2.to_numpy() - b.r2.to_numpy(), a.n_train.to_numpy(), a.n_test.to_numpy())
        rows.append(dict(target=t, mean=mean, lo=lo, hi=hi, p=p))
    R = pd.DataFrame(rows); R["q_bh"] = bh(R.p); R["q_by"] = by(R.p)
    M = R.merge(C, on="target")
    out["max_abs_diff_confirmatory"] = float(np.max(np.abs(np.r_[
        M["mean"] - M.d_r2_mean, M.lo - M.d_r2_ci_low, M.hi - M.d_r2_ci_high,
        M.p - M.d_r2_p, M.q_bh - M.d_r2_q_bh, M.q_by - M.d_r2_q_by])))
    out["informative_bh"] = sorted(M.loc[(M.q_bh <= 0.05) & (M["mean"] > 0), "target"].tolist())
    out["practically_significant"] = sorted(M.loc[M["mean"] >= 0.05, "target"].tolist())
    out["table"] = M[["target", "mean", "lo", "hi", "p", "q_bh", "q_by"]].round(6).to_dict("records")
    return out


if __name__ == "__main__":
    for a in sys.argv[1:]:
        o = check(Path(a).expanduser())
        tab = o.pop("table")
        for k, v in o.items():
            print(f"{k}: {v}")
        print(pd.DataFrame(tab).to_string(index=False)); print()
