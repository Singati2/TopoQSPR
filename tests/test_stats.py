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
    out = stats.corrected_ttest([0.1, 0.1, 0.1], 8, 2)
    assert out["p"] == 0.0 and out["ci_low"] == pytest.approx(0.1)


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
