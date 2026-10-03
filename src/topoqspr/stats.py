"""Statistics for the evaluation protocol.

Contents
--------
corrected_ttest
    Corrected repeated k-fold CV t-test and CI for a paired metric difference.
benjamini_hochberg
    BH step-up adjusted p-values (q-values) and rejections.
spearman_ci
    Spearman rho with a Fisher-z confidence interval.
y_scrambling
    Permutation (y-randomisation) test of any scoring function.
permutation_pvalue
    Empirical p-value (+1 correction) and 95th-percentile rule for a given null.
percentile_interval
    Percentile bootstrap interval that never drops failed resamples silently.
group_bootstrap
    Cluster (group) bootstrap resample of row indices.
ApplicabilityDomain
    kNN-distance or leverage applicability-domain flag.

References
----------
Nadeau C, Bengio Y (2003) Inference for the generalization error.
  Machine Learning 52:239-281. doi:10.1023/A:1024068626366
Bouckaert RR, Frank E (2004) Evaluating the replicability of significance
  tests for comparing learning algorithms. PAKDD 2004, LNCS 3056:3-12.
  doi:10.1007/978-3-540-24775-3_3
Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate.
  J R Stat Soc B 57:289-300. doi:10.1111/j.2517-6161.1995.tb02031.x
Fieller EC, Hartley HO, Pearson ES (1957) Tests for rank correlation
  coefficients. I. Biometrika 44:470-481. doi:10.1093/biomet/44.3-4.470
Bonett DG, Wright TA (2000) Sample size requirements for estimating Pearson,
  Kendall and Spearman correlations. Psychometrika 65:23-28.
  doi:10.1007/BF02294183
Phipson B, Smyth GK (2010) Permutation p-values should never be zero.
  Stat Appl Genet Mol Biol 9:39. doi:10.2202/1544-6115.1585
Ruecker C, Ruecker G, Meringer M (2007) y-Randomization and its variants in
  QSPR/QSAR. J Chem Inf Model 47:2345-2357. doi:10.1021/ci700157b
Tropsha A, Gramatica P, Gombar VK (2003) The importance of being earnest.
  QSAR Comb Sci 22:69-77. doi:10.1002/qsar.200390007
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

import numpy as np
from scipy import stats as _st
from scipy.spatial.distance import cdist as _cdist

__all__ = [
    "corrected_ttest",
    "benjamini_hochberg",
    "benjamini_yekutieli",
    "spearman_ci",
    "y_scrambling",
    "permutation_pvalue",
    "percentile_interval",
    "group_bootstrap",
    "ApplicabilityDomain",
]


def corrected_ttest(diffs: Sequence[float], n_train, n_test, alpha: float = 0.05,
                    nan_policy: str = "propagate") -> dict:
    """Corrected repeated k-fold CV t-test for a paired difference.

    For J = r*k paired per-fold differences d_j (r repeats of k-fold CV,
    every model evaluated on the same folds):

        d_bar   = (1/J) sum_j d_j
        s^2     = (1/(J-1)) sum_j (d_j - d_bar)^2
        se_corr = sqrt((1/J + n2/n1) * s^2)
        t       = d_bar / se_corr,  df = J - 1
        CI      = d_bar +/- t_{1-alpha/2, J-1} * se_corr

    n1 / n2 are the training / test sizes of one split (Nadeau & Bengio 2003,
    "corrected resampled t-test"; Bouckaert & Frank 2004 apply it to r x k
    repeated CV with df = r*k - 1). The term n2/n1 inflates the naive variance
    s^2/J to account for overlapping training sets. With unequal folds
    (grouped CV) the ratio is mean(n_test) / mean(n_train).

    If s^2 == 0 the CI collapses to [d_bar, d_bar] and p = 0 if d_bar != 0
    else 1. ``nan_policy``: "propagate" (any NaN -> NaN result) or "omit"
    (drop NaN differences and their sizes).

    Returns dict with keys mean, sd, se, t, df, p, ci_low, ci_high, n, ratio.
    """
    d = np.asarray(diffs, dtype=float).ravel()
    n1 = np.broadcast_to(np.asarray(n_train, dtype=float), d.shape)
    n2 = np.broadcast_to(np.asarray(n_test, dtype=float), d.shape)
    if nan_policy == "omit":
        ok = np.isfinite(d)
        d, n1, n2 = d[ok], n1[ok], n2[ok]
    elif nan_policy != "propagate":
        raise ValueError("nan_policy must be 'propagate' or 'omit'")
    J = d.size
    out = dict(mean=math.nan, sd=math.nan, se=math.nan, t=math.nan, df=J - 1,
               p=math.nan, ci_low=math.nan, ci_high=math.nan, n=J, ratio=math.nan)
    if J < 2 or not np.all(np.isfinite(d)):
        return out
    ratio = float(n2.mean() / n1.mean())
    mean = float(d.mean())
    var = 0.0 if np.ptp(d) == 0 else float(d.var(ddof=1))
    se = math.sqrt((1.0 / J + ratio) * var)
    out.update(mean=mean, sd=math.sqrt(var), se=se, ratio=ratio)
    if se == 0.0:
        out.update(t=math.copysign(math.inf, mean) if mean != 0 else math.nan,
                   p=0.0 if mean != 0 else 1.0, ci_low=mean, ci_high=mean)
        return out
    t = mean / se
    crit = float(_st.t.ppf(1 - alpha / 2, J - 1))
    out.update(t=t, p=float(2 * _st.t.sf(abs(t), J - 1)),
               ci_low=mean - crit * se, ci_high=mean + crit * se)
    return out


def _family_size(p: np.ndarray, m: int | None) -> int:
    k = int(np.sum(~np.isnan(p)))
    if m is None:
        return k
    m = int(m)
    if m < k:
        raise ValueError(f"family size m={m} is smaller than the {k} finite p-values")
    return m


def benjamini_hochberg(pvals: Sequence[float], q: float = 0.05, m: int | None = None):
    """Benjamini-Hochberg step-up procedure.

    For the k non-NaN p-values sorted p_(1) <= ... <= p_(k) and a family of
    m >= k hypotheses, the adjusted value is
    q_(i) = min_{j >= i} min(1, m * p_(j) / j); hypothesis i is rejected iff
    q_(i) <= q (equivalently: reject H_(1..r) with r = max{i : p_(i) <= i q / m}).
    NaN p-values stay NaN and are never rejected. By default (``m=None``)
    m = k, i.e. a NaN does not count towards m. Passing the pre-specified
    family size ``m`` keeps the untestable members in m (equivalent to
    giving them p = 1), so a member that fails to evaluate cannot loosen
    the correction for the others.

    Returns (adjusted, reject) as numpy arrays in the input order.
    """
    p = np.asarray(pvals, dtype=float).ravel()
    if np.any((p[~np.isnan(p)] < 0) | (p[~np.isnan(p)] > 1)):
        raise ValueError("p-values must lie in [0, 1]")
    m = _family_size(p, m)
    adj = np.full(p.shape, np.nan)
    ok = np.flatnonzero(~np.isnan(p))
    k = ok.size
    if k:
        order = ok[np.argsort(p[ok], kind="mergesort")]
        ranked = p[order] * m / np.arange(1, k + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        adj[order] = np.minimum(ranked, 1.0)
    reject = np.zeros(p.shape, dtype=bool)
    reject[ok] = adj[ok] <= q
    return adj, reject


def benjamini_yekutieli(pvals: Sequence[float], q: float = 0.05, m: int | None = None):
    """Benjamini-Yekutieli step-up procedure (valid under arbitrary dependence).

    Same as :func:`benjamini_hochberg` with every p-value multiplied by
    c(m) = sum_{i=1..m} 1/i (Benjamini & Yekutieli 2001, Thm 1.3):
    q_(i) = min_{j >= i} min(1, c(m) * m * p_(j) / j). Used as the sensitivity
    analysis for BH because the confirmatory endpoints share drugs and folds.
    ``m`` as in :func:`benjamini_hochberg` (default: the number of finite p).

    Returns (adjusted, reject) in the input order; NaN p-values stay NaN.
    """
    p = np.asarray(pvals, dtype=float).ravel()
    m = _family_size(p, m)
    if m == 0 or np.all(np.isnan(p)):
        return np.full(p.shape, np.nan), np.zeros(p.shape, dtype=bool)
    c = float(np.sum(1.0 / np.arange(1, m + 1)))
    scaled = np.where(np.isnan(p), np.nan, np.minimum(1.0, p * c))
    adj, _ = benjamini_hochberg(scaled, q, m=m)
    # BH on c*p caps at 1 before the step-up; recompute the reject set on the adjusted values
    adj = np.minimum(1.0, adj)
    return adj, np.where(np.isnan(adj), False, adj <= q)


def spearman_ci(x: Sequence[float], y: Sequence[float], alpha: float = 0.05,
                method: str = "fieller") -> dict:
    """Spearman rho with a Fisher-z confidence interval.

    z = atanh(rho); CI = tanh(z -/+ z_{1-alpha/2} * se) with
      se = sqrt(1.06 / (n - 3))              (method="fieller"; Fieller et al. 1957)
      se = sqrt((1 + rho^2 / 2) / (n - 3))   (method="bonett-wright"; 2000)
    The p-value is scipy's (t approximation). Pairs with a NaN are dropped.
    Returns rho, ci_low, ci_high, p, n and ``reason`` (empty if defined).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = int(x.size)
    out = dict(rho=math.nan, ci_low=math.nan, ci_high=math.nan, p=math.nan, n=n, reason="")
    if n < 4:
        out["reason"] = f"n={n} < 4"
        return out
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        out["reason"] = "constant input"
        return out
    res = _st.spearmanr(x, y)
    rho, p = float(res.statistic), float(res.pvalue)
    if method == "fieller":
        se = math.sqrt(1.06 / (n - 3))
    elif method == "bonett-wright":
        se = math.sqrt((1 + rho ** 2 / 2) / (n - 3))
    else:
        raise ValueError("method must be 'fieller' or 'bonett-wright'")
    zc = float(_st.norm.ppf(1 - alpha / 2))
    if abs(rho) >= 1.0:
        lo = hi = rho
    else:
        z = math.atanh(rho)
        lo, hi = math.tanh(z - zc * se), math.tanh(z + zc * se)
    out.update(rho=rho, ci_low=lo, ci_high=hi, p=p)
    return out


def y_scrambling(score_fn: Callable[[np.ndarray, np.ndarray], float], X, y,
                 n_permutations: int = 100, seed: int = 0) -> dict:
    """Y-randomisation (permutation) test.

    The observed statistic s_obs = score_fn(X, y) is compared with
    s_b = score_fn(X, pi_b(y)), b = 1..B, where pi_b are permutations drawn
    from ``numpy.random.default_rng(seed)``. score_fn must run the *whole*
    pipeline (including model selection) so that the null reflects every
    data-dependent step (Ruecker et al. 2007). Larger is better.

        p = (1 + #{b : s_b >= s_obs}) / (1 + B)       (Phipson & Smyth 2010)

    NaN null statistics are kept in ``null`` but excluded from the count and
    from B (``n_valid``).
    """
    y = np.asarray(y, dtype=float)
    rng = np.random.default_rng(seed)
    obs = float(score_fn(X, y))
    null = np.array([float(score_fn(X, y[rng.permutation(y.size)]))
                     for _ in range(n_permutations)])
    valid = null[np.isfinite(null)]
    p = (1 + int(np.sum(valid >= obs))) / (1 + valid.size) if np.isfinite(obs) else math.nan
    return dict(observed=obs, null=null, p=p, n_valid=int(valid.size),
                null_mean=float(valid.mean()) if valid.size else math.nan,
                null_max=float(valid.max()) if valid.size else math.nan)


def permutation_pvalue(observed: float, null: Sequence[float], quantile: float = 0.95) -> dict:
    """Empirical permutation p-value and the "uninformative" rule.

    With s_obs the observed statistic and s_1..s_K the statistics of K
    permuted targets (larger is better),

        p = (1 + #{b : s_b >= s_obs}) / (1 + K_valid)     (Phipson & Smyth 2010)
        q95 = numpy.quantile(valid null, quantile)        (linear interpolation)
        uninformative = s_obs <= q95                       (EXPERIMENT_PLAN 7)

    NaN null statistics (a permutation whose nested CV failed) are reported as
    ``n_failed`` and excluded from K_valid; a NaN observed statistic or an
    empty valid null makes p, q95 and the flag NaN/None with ``reason``.
    Returns observed, p, n, n_valid, n_failed, null_mean, null_sd, null_q95,
    null_max, uninformative (bool or None), reason.
    """
    null = np.asarray(null, dtype=float).ravel()
    valid = null[np.isfinite(null)]
    out = dict(observed=float(observed), p=math.nan, n=int(null.size), n_valid=int(valid.size),
               n_failed=int(null.size - valid.size), null_mean=math.nan, null_sd=math.nan,
               null_q95=math.nan, null_max=math.nan, uninformative=None, reason="")
    why = []
    if out["n_failed"]:
        why.append(f"{out['n_failed']} of {null.size} permutations failed (excluded)")
    if valid.size:
        out.update(null_mean=float(valid.mean()),
                   null_sd=float(valid.std(ddof=1)) if valid.size > 1 else math.nan,
                   null_q95=float(np.quantile(valid, quantile)), null_max=float(valid.max()))
    if not math.isfinite(out["observed"]):
        why.append("observed statistic undefined")
    elif valid.size == 0:
        why.append("no valid permutation")
    else:
        out["p"] = (1 + int(np.sum(valid >= out["observed"]))) / (1 + valid.size)
        out["uninformative"] = bool(out["observed"] <= out["null_q95"])
    out["reason"] = "; ".join(why)
    return out


def percentile_interval(values: Sequence[float], alpha: float = 0.05) -> dict:
    """Percentile bootstrap interval [Q(alpha/2), Q(1 - alpha/2)] of ``values``.

    Quantiles use ``numpy.quantile`` (linear interpolation). A NaN value (a
    resample whose procedure failed or whose metric is undefined) is **not**
    dropped: the interval is then NaN with a reason, because keeping only the
    resamples that succeeded would condition on success. Returns ci_low,
    ci_high, mean, sd, n, n_valid, reason.
    """
    v = np.asarray(values, dtype=float).ravel()
    ok = np.isfinite(v)
    out = dict(ci_low=math.nan, ci_high=math.nan, mean=math.nan, sd=math.nan,
               n=int(v.size), n_valid=int(ok.sum()), reason="")
    if v.size == 0:
        out["reason"] = "no resamples"
        return out
    if not ok.all():
        out["reason"] = f"{int((~ok).sum())} of {v.size} resamples undefined; interval not computed"
        return out
    lo, hi = np.quantile(v, [alpha / 2, 1 - alpha / 2])
    out.update(ci_low=float(lo), ci_high=float(hi), mean=float(v.mean()),
               sd=float(v.std(ddof=1)) if v.size > 1 else math.nan)
    return out


def group_bootstrap(groups: Sequence, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """One cluster-bootstrap resample: draw G groups with replacement from the
    G distinct ``groups`` and take every row of each drawn group.

    Returns ``(rows, draw)``: ``rows`` are indices into the original rows (a
    row appears once per draw of its group) and ``draw`` labels every resampled
    row with the position of its group's draw (0..G-1). Rows from different
    draws of the same group are copies of one another; callers that run CV on
    the resample must keep all copies of an original group together (pass the
    original group label as the CV group), otherwise a copy lands in the
    training and the test fold at once. Groups are taken in sorted order of
    their labels, so the result depends only on ``rng`` and the labels.
    """
    g = np.asarray(groups, dtype=object).astype(str)
    uniq, inv = np.unique(g, return_inverse=True)
    members = [np.flatnonzero(inv == i) for i in range(uniq.size)]
    pick = rng.integers(0, uniq.size, size=uniq.size)
    rows = np.concatenate([members[i] for i in pick])
    draw = np.concatenate([np.full(members[i].size, j) for j, i in enumerate(pick)])
    return rows, draw


class ApplicabilityDomain:
    """Applicability-domain (AD) flag fitted on training rows only.

    Training features are median-imputed (training medians), columns constant
    in training are dropped, and the rest standardised with training mean/SD.

    method="knn" (default; usable when p >= n):
        d(x) = mean Euclidean distance from x to its k nearest training rows.
        For training rows the row itself is excluded (leave-one-out).
        Threshold = the ``quantile`` (default 0.95) of training d values;
        x is in the domain iff d(x) <= threshold.
    method="leverage" (Tropsha et al. 2003 "Williams plot" convention):
        h(x) = x~ (X~' X~)^+ x~', X~ = [1, Z] (intercept + standardised
        training matrix), threshold h* = 3 (p' + 1) / n with p' the number of
        retained columns. Requires n > p' + 1, else ``ValueError`` because all
        training leverages would equal 1 (use knn instead).
    """

    def __init__(self, method: str = "knn", k: int = 5, quantile: float = 0.95):
        if method not in ("knn", "leverage"):
            raise ValueError("method must be 'knn' or 'leverage'")
        self.method, self.k, self.quantile = method, k, quantile

    def _z(self, X):
        X = np.asarray(X, dtype=float)
        X = np.where(np.isnan(X), self.median_, X)
        return (X[:, self.keep_] - self.mean_) / self.sd_

    def fit(self, X):
        X = np.asarray(X, dtype=float)
        if X.ndim != 2 or X.shape[1] == 0:
            raise ValueError("AD needs at least one feature column")
        with np.errstate(all="ignore"):
            med = np.nanmedian(X, axis=0)
        self.median_ = np.where(np.isnan(med), 0.0, med)
        Xi = np.where(np.isnan(X), self.median_, X)
        sd = Xi.std(axis=0)
        self.keep_ = np.flatnonzero(sd > 0)
        if self.keep_.size == 0:
            raise ValueError("all features constant in training rows")
        self.mean_ = Xi[:, self.keep_].mean(axis=0)
        self.sd_ = sd[self.keep_]
        Z = self._z(X)
        n, p = Z.shape
        if self.method == "knn":
            k = min(self.k, n - 1)
            if k < 1:
                raise ValueError("knn AD needs at least 2 training rows")
            self.k_ = k
            self.train_ = Z
            D = _cdist(Z, Z)
            np.fill_diagonal(D, np.inf)
            d = np.sort(D, axis=1)[:, :k].mean(axis=1)
            self.threshold_ = float(np.quantile(d, self.quantile))
        else:
            if n <= p + 1:
                raise ValueError(f"leverage AD undefined for n={n} <= p'+1={p + 1}")
            Xt = np.hstack([np.ones((n, 1)), Z])
            self.xtx_pinv_ = np.linalg.pinv(Xt.T @ Xt)
            self.threshold_ = 3.0 * (p + 1) / n
        return self

    def score(self, X) -> np.ndarray:
        """kNN distance or leverage of each row of X (larger = further out)."""
        Z = self._z(X)
        if self.method == "knn":
            D = _cdist(Z, self.train_)
            return np.sort(D, axis=1)[:, : self.k_].mean(axis=1)
        Xt = np.hstack([np.ones((Z.shape[0], 1)), Z])
        return np.einsum("ij,jk,ik->i", Xt, self.xtx_pinv_, Xt)

    def in_domain(self, X) -> np.ndarray:
        return self.score(X) <= self.threshold_
