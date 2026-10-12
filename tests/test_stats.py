"""Tests for topoqspr.stats against hand-computed and published reference values."""

from __future__ import annotations

import math

import numpy as np
import pytest

from topoqspr import stats


# ------------------------------------------------------- corrected t-test
def test_corrected_ttest_hand_example():
    # J = 5 (k = 5, r = 1); d = 1..5; n_train = 80, n_test = 20 -> n2/n1 = 1/4
    # mean 3, s^2 = 2.5, var_corr = (1/5 + 1/4) * 2.5 = 1.125, se = 1.0606602
    # t = 3 / 1.0606602 = 2.8284271, df = 4, t_.975,4 = 2.7764451
    # CI = 3 -/+ 2.7764451 * 1.0606602 = [0.0551353, 5.9448647]
    out = stats.corrected_ttest([1, 2, 3, 4, 5], n_train=80, n_test=20)
    assert out["mean"] == pytest.approx(3.0)
    assert out["se"] == pytest.approx(math.sqrt(1.125))
    assert out["t"] == pytest.approx(2.8284271, rel=1e-7)
    assert out["df"] == 4
    assert out["ci_low"] == pytest.approx(0.0551353, abs=1e-6)
    assert out["ci_high"] == pytest.approx(5.9448647, abs=1e-6)
    # two-sided p for t = 2.8284271 with df = 4 (just beyond t_.975,4 = 2.776)
    assert out["p"] == pytest.approx(0.04742, abs=2e-4)


def test_corrected_ttest_is_wider_than_naive():
    rng = np.random.default_rng(0)
    d = rng.normal(0.02, 0.05, size=25)
    out = stats.corrected_ttest(d, 80, 20)
    naive_se = d.std(ddof=1) / math.sqrt(25)
    assert out["se"] / naive_se == pytest.approx(math.sqrt(1 + 25 / 4))
    assert out["ratio"] == pytest.approx(0.25)


def test_corrected_ttest_unequal_folds_and_nan():
    out = stats.corrected_ttest([0.1, 0.2, 0.3], n_train=[8, 9, 10], n_test=[2, 1, 2])
    assert out["ratio"] == pytest.approx((5 / 3) / 9)
    assert math.isnan(stats.corrected_ttest([0.1, np.nan, 0.3], 8, 2)["mean"])
    om = stats.corrected_ttest([0.1, np.nan, 0.3], 8, 2, nan_policy="omit")
    assert om["n"] == 2 and om["mean"] == pytest.approx(0.2)


def test_corrected_ttest_zero_variance():
    out = stats.corrected_ttest([0.0, 0.0, 0.0], 8, 2)
    assert out["p"] == 1.0 and out["ci_low"] == out["ci_high"] == 0.0


def test_corrected_ttest_identical_nonzero_differences_are_degenerate():
    # audit 2026-10-12: 25 identical non-zero fold differences (5 x 5 CV) used to
    # give p = 0 and the zero-width CI [0.1, 0.1]. The variance is not estimable,
    # so the test must not claim significance or a precise interval.
    out = stats.corrected_ttest([0.1] * 25, n_train=80, n_test=20)
    assert out["mean"] == pytest.approx(0.1) and out["se"] == 0.0
    assert math.isnan(out["p"]) and math.isnan(out["t"])
    assert math.isnan(out["ci_low"]) and math.isnan(out["ci_high"])
    assert "not estimable" in out["reason"]
    # any spread at all makes the variance estimable again
    out = stats.corrected_ttest([0.1] * 24 + [0.1 + 1e-9], 80, 20)
    assert math.isfinite(out["p"]) and out["reason"] == ""


# --------------------------------------------------- Benjamini-Hochberg
BH1995 = [0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344, 0.0459,
          0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.0000]


def test_bh_published_example():
    # Benjamini & Hochberg (1995), section 4: 15 p-values, q = 0.05 -> the
    # four smallest are rejected (Bonferroni rejects three).
    adj, rej = stats.benjamini_hochberg(BH1995, q=0.05)
    assert rej.sum() == 4 and rej[:4].all()
    # hand values q_(i) = min_{j>=i} 15 p_(j) / j:
    # 15p/j for j = 4..9: .035625 .0603 .0695 .063857 .0645 .0765
    assert adj[3] == pytest.approx(0.035625)
    assert adj[4] == pytest.approx(0.0603)
    assert adj[5] == pytest.approx(15 * 0.0298 / 7)
    assert adj[8] == pytest.approx(0.0765)
    assert adj[2] == pytest.approx(0.0095)
    assert adj[0] == pytest.approx(0.0015)
    assert adj[-1] == 1.0


def test_bh_small_hand_example_order_and_nan():
    adj, rej = stats.benjamini_hochberg([0.01, 0.04, np.nan, 0.03, 0.005], q=0.03)
    # sorted 0.005, 0.01, 0.03, 0.04 (m = 4): 4p/j = 0.02, 0.02, 0.04, 0.04
    np.testing.assert_allclose(adj[[0, 1, 3, 4]], [0.02, 0.04, 0.04, 0.02])
    assert np.isnan(adj[2]) and not rej[2]
    assert rej.tolist() == [True, False, False, False, True]


def test_bh_matches_step_up_definition():
    rng = np.random.default_rng(3)
    for _ in range(50):
        p = np.concatenate([rng.uniform(size=20), rng.uniform(0, 0.01, size=5)])
        _, rej = stats.benjamini_hochberg(p, q=0.1)
        ps = np.sort(p)
        ok = np.flatnonzero(ps <= np.arange(1, p.size + 1) * 0.1 / p.size)
        k = ok.max() + 1 if ok.size else 0
        assert rej.sum() == k
        if k:
            assert set(np.flatnonzero(rej)) == set(np.argsort(p, kind="mergesort")[:k])


def test_bh_rejects_bad_input():
    with pytest.raises(ValueError):
        stats.benjamini_hochberg([0.5, 1.2])


# -------------------------------------------------------------- Spearman
def test_spearman_ci_hand_example():
    x = np.arange(20, dtype=float)
    y = x.copy()
    y[[3, 7]] = y[[7, 3]]
    y[[10, 15]] = y[[15, 10]]
    out = stats.spearman_ci(x, y)
    # rho = 1 - 6 sum d^2 / (n (n^2 - 1)) = 1 - 6*(16+16+25+25)/(20*399)
    rho = 1 - 6 * 82 / (20 * 399)
    assert out["rho"] == pytest.approx(rho)
    se = math.sqrt(1.06 / 17)
    assert out["ci_low"] == pytest.approx(math.tanh(math.atanh(rho) - 1.959964 * se), abs=1e-6)
    assert out["ci_high"] == pytest.approx(math.tanh(math.atanh(rho) + 1.959964 * se), abs=1e-6)
    bw = stats.spearman_ci(x, y, method="bonett-wright")
    assert bw["ci_low"] < out["ci_low"]  # (1 + rho^2/2) > 1.06 for rho ~ 0.95


def test_spearman_edge_cases():
    assert stats.spearman_ci([1, 2, 3], [1, 2, 3])["reason"].startswith("n=3")
    assert stats.spearman_ci([1, 1, 1, 1, 1], [1, 2, 3, 4, 5])["reason"] == "constant input"
    out = stats.spearman_ci([1, 2, np.nan, 4, 5, 6], [2, 1, 3, 4, 6, 5])
    assert out["n"] == 5


# ------------------------------------------------------------ y-scrambling
def _abs_corr(X, y):
    return abs(np.corrcoef(np.asarray(X)[:, 0], y)[0, 1])


def test_y_scrambling_signal_and_noise():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(50, 1))
    y = X[:, 0] + rng.normal(scale=0.5, size=50)
    sig = stats.y_scrambling(_abs_corr, X, y, n_permutations=99, seed=1)
    assert sig["p"] == pytest.approx(1 / 100)
    noise = stats.y_scrambling(_abs_corr, X, rng.normal(size=50), n_permutations=99, seed=1)
    assert noise["p"] > 0.05
    again = stats.y_scrambling(_abs_corr, X, y, n_permutations=99, seed=1)
    np.testing.assert_array_equal(sig["null"], again["null"])


def test_y_scrambling_with_nested_pipeline():
    from topoqspr.evaluate import nested_r2_statistic

    rng = np.random.default_rng(2)
    X = rng.normal(size=(40, 3))
    y = 2 * X[:, 0] + rng.normal(scale=0.5, size=40)
    fn = nested_r2_statistic(families=("ridge",), seeds=(0,), outer_splits=3, ad=None)
    out = stats.y_scrambling(fn, X, y, n_permutations=9, seed=0)
    assert out["observed"] > 0.7
    assert out["null_max"] < 0.3
    assert out["p"] == pytest.approx(0.1)


# ------------------------------------------------- applicability domain
def test_ad_knn():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(100, 3))
    ad = stats.ApplicabilityDomain("knn", k=5).fit(X)
    assert ad.in_domain(np.zeros((1, 3)))[0]
    assert not ad.in_domain(np.full((1, 3), 10.0))[0]
    # ~95% of fresh in-distribution points are in the domain
    cov = ad.in_domain(rng.normal(size=(2000, 3))).mean()
    assert 0.85 < cov <= 1.0


def test_ad_leverage_matches_hat_matrix():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(30, 2))
    ad = stats.ApplicabilityDomain("leverage").fit(X)
    Z = (X - X.mean(0)) / X.std(0)
    Xt = np.hstack([np.ones((30, 1)), Z])
    H = Xt @ np.linalg.inv(Xt.T @ Xt) @ Xt.T
    np.testing.assert_allclose(ad.score(X), np.diag(H), atol=1e-10)
    assert ad.threshold_ == pytest.approx(3 * 3 / 30)
    assert np.diag(H).sum() == pytest.approx(3)
    with pytest.raises(ValueError):
        stats.ApplicabilityDomain("leverage").fit(rng.normal(size=(5, 4)))


def test_ad_handles_nan_and_constant_columns():
    X = np.column_stack([np.arange(10.0), np.ones(10), np.r_[np.nan, np.arange(9.0)]])
    ad = stats.ApplicabilityDomain("knn", k=3).fit(X)
    assert list(ad.keep_) == [0, 2]
    assert ad.score(X).shape == (10,)


def test_benjamini_yekutieli_matches_definition():
    """Brute-force check of q_(i) = min_{j>=i} min(1, c(m) m p_(j) / j), c(m) = sum 1/k."""
    from topoqspr.stats import benjamini_yekutieli

    rng = np.random.default_rng(1)
    for _ in range(200):
        m = int(rng.integers(1, 40))
        p = rng.uniform(0, 1, m) ** rng.uniform(0.5, 4)
        adj, rej = benjamini_yekutieli(p, 0.05)
        c = sum(1.0 / k for k in range(1, m + 1))
        o = np.argsort(p)
        ps = p[o]
        ref = np.empty(m)
        for i in range(m):
            ref[o[i]] = min(1.0, min(c * m * ps[j] / (j + 1) for j in range(i, m)))
        assert np.allclose(adj, ref, atol=1e-12)
        assert (rej == (ref <= 0.05)).all()


def test_benjamini_yekutieli_is_never_more_liberal_than_bh():
    from topoqspr.stats import benjamini_hochberg, benjamini_yekutieli

    rng = np.random.default_rng(2)
    p = rng.uniform(0, 0.2, 50)
    assert (benjamini_yekutieli(p)[0] >= benjamini_hochberg(p)[0] - 1e-15).all()


# ------------------------------------------- regressions (audit:math findings)
def test_benjamini_yekutieli_rejects_p_outside_unit_interval():
    # regression: BY clipped c*p to 1 before BH's range check, so p > 1 passed silently
    for bad in ([1.5, 0.01], [-0.1, 0.2]):
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            stats.benjamini_yekutieli(bad)
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            stats.benjamini_hochberg(bad)
    adj, rej = stats.benjamini_yekutieli([np.nan, 1.0, 0.0])  # NaN and the bounds are fine
    assert math.isnan(adj[0]) and not rej[0] and rej[2]


@pytest.mark.parametrize("y, sign", [([2, 3, 5, 9], 1.0), ([9, 5, 3, 2], -1.0)])
def test_spearman_ci_undefined_at_perfect_rank_correlation(y, sign):
    # regression: |rho| = 1 gave a zero-width CI [rho, rho] with an empty reason
    out = stats.spearman_ci([1, 2, 3, 4], y)
    assert out["rho"] == sign and out["n"] == 4
    assert math.isnan(out["ci_low"]) and math.isnan(out["ci_high"])
    assert out["reason"] == "rho = +/-1: Fisher-z interval undefined"
    bw = stats.spearman_ci([1, 2, 3, 4], y, method="bonett-wright")
    assert math.isnan(bw["ci_low"]) and bw["reason"] == out["reason"]


def test_ad_knn_groups_exclude_same_group_neighbours():
    # regression: duplicate training rows (connectivity groups: identical 2D features)
    # contributed distance 0 to each other's kNN distance and biased the threshold low
    rng = np.random.default_rng(5)
    base = rng.normal(size=(60, 3))
    X = np.vstack([base, base[:30]])                       # 30 molecules duplicated
    groups = np.r_[np.arange(60), np.arange(30)].astype(str)
    plain = stats.ApplicabilityDomain("knn", k=5).fit(X)
    grp = stats.ApplicabilityDomain("knn", k=5).fit(X, groups=groups)
    assert grp.threshold_ > plain.threshold_
    # threshold = quantile of the leave-one-group-out kNN distances
    Z = grp.train_
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    D[groups[:, None] == groups[None, :]] = np.inf
    d = np.sort(D, axis=1)[:, :5].mean(axis=1)
    assert np.all(np.isfinite(d))
    assert grp.threshold_ == pytest.approx(float(np.quantile(d, 0.95)))
    # singleton groups reproduce the leave-one-out rule exactly
    solo = stats.ApplicabilityDomain("knn", k=5).fit(base, groups=np.arange(60).astype(str))
    assert solo.threshold_ == stats.ApplicabilityDomain("knn", k=5).fit(base).threshold_
    # coverage of fresh in-distribution rows moves towards the nominal 95% with groups
    rng2 = np.random.default_rng(6)
    cov_plain, cov_grp = [], []
    for _ in range(40):
        b = rng2.normal(size=(60, 3))
        Xr, test = np.vstack([b, b[:30]]), rng2.normal(size=(400, 3))
        cov_plain.append(stats.ApplicabilityDomain("knn", k=5).fit(Xr).in_domain(test).mean())
        cov_grp.append(stats.ApplicabilityDomain("knn", k=5).fit(Xr, groups=groups)
                       .in_domain(test).mean())
    assert np.mean(cov_grp) > np.mean(cov_plain)
    with pytest.raises(ValueError, match="2 groups"):
        stats.ApplicabilityDomain("knn").fit(base[:4], groups=["a"] * 4)
    with pytest.raises(ValueError, match="labels"):
        stats.ApplicabilityDomain("knn").fit(base, groups=["a", "b"])


def test_nested_cv_passes_groups_to_knn_ad():
    # the outer-fold AD is fitted with the training rows' CV groups
    from unittest import mock

    from topoqspr import evaluate as ev

    rng = np.random.default_rng(1)
    X = rng.normal(size=(40, 2))
    y = X[:, 0] + rng.normal(scale=0.3, size=40)
    groups = np.array([f"g{i // 2}" for i in range(40)], dtype=object)
    seen = []
    real_fit = stats.ApplicabilityDomain.fit

    def spy(self, Xtr, groups=None):
        seen.append(None if groups is None else np.asarray(groups).copy())
        return real_fit(self, Xtr, groups=groups)

    with mock.patch.object(stats.ApplicabilityDomain, "fit", spy):
        ev.nested_cv(X, y, ("ridge",), outer_splits=3, seeds=(0,), inner_splits=2,
                     groups=groups, ad="knn", min_n=10, min_n_test=3)
    assert len(seen) == 3 and all(g is not None and g.size > 0 for g in seen)
