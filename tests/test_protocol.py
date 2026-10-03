"""Tests for the confirmatory protocol stage (pipeline.permutation_control,
pipeline.paired_bootstrap, pipeline.ad_report, pipeline.confirmatory_table and
the protocol section of run_study) and the stats helpers it uses. Offline and
synthetic only: random matrices, and the synthetic SMILES file
tests/fixtures/pipeline/tiny.csv with computed (RDKit) targets. Learners are cut
to OLS / ridge, one repeat, three outer folds.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from topoqspr import evaluate as ev
from topoqspr import pipeline as pl
from topoqspr import stats

REPO = Path(__file__).resolve().parents[1]
FIX = REPO / "tests" / "fixtures" / "pipeline"
TS = "2026-10-03T00:00:00Z"
KW = dict(outer_splits=3, seeds=[0], inner_splits=2, random_state=0, min_n=10, min_n_test=3)


# ------------------------------------------------------------------ stats
def test_permutation_pvalue_plus_one_and_rule():
    r = stats.permutation_pvalue(0.5, [0.1, 0.2, 0.6, math.nan])
    assert r["p"] == pytest.approx((1 + 1) / (1 + 3))
    assert r["n_valid"] == 3 and r["n_failed"] == 1 and "1 of 4" in r["reason"]
    assert r["null_q95"] == pytest.approx(np.quantile([0.1, 0.2, 0.6], 0.95))
    assert r["uninformative"] is True  # 0.5 <= q95 (0.56)
    strong = stats.permutation_pvalue(0.9, [0.0] * 19)
    assert strong["p"] == pytest.approx(1 / 20) and strong["uninformative"] is False
    assert stats.permutation_pvalue(0.9, [0.0, 0.9])["p"] == pytest.approx(2 / 3)  # ties count
    undef = stats.permutation_pvalue(math.nan, [0.1, 0.2])
    assert math.isnan(undef["p"]) and undef["uninformative"] is None and undef["reason"]
    none = stats.permutation_pvalue(0.3, [math.nan, math.nan])
    assert math.isnan(none["p"]) and "no valid permutation" in none["reason"]


def test_permutation_pvalue_is_calibrated_for_exchangeable_statistics():
    # Under H0 the observed and the K null statistics are exchangeable, so
    # P(p <= a) <= a for every a, with equality at a = j / (K + 1).
    rng = np.random.default_rng(7)
    K, reps = 19, 4000
    draws = rng.normal(size=(reps, K + 1))
    p = np.array([stats.permutation_pvalue(d[0], d[1:])["p"] for d in draws])
    for a in (0.05, 0.1, 0.25, 0.5):
        assert abs(np.mean(p <= a) - a) < 0.025
    assert abs(p.mean() - 0.525) < 0.02  # E[p] = (K + 2) / (2 (K + 1))


def test_percentile_interval_never_drops_failures():
    v = np.arange(1, 101, dtype=float)
    r = stats.percentile_interval(v, 0.05)
    assert r["ci_low"] == pytest.approx(np.quantile(v, 0.025))
    assert r["ci_high"] == pytest.approx(np.quantile(v, 0.975))
    bad = stats.percentile_interval([1.0, math.nan, 2.0])
    assert math.isnan(bad["ci_low"]) and "1 of 3" in bad["reason"] and bad["n_valid"] == 2
    assert "no resamples" in stats.percentile_interval([])["reason"]


def test_group_bootstrap_keeps_groups_whole_and_is_seeded():
    groups = ["a", "a", "b", "c", "c", "c", "d"]
    rows, draw = stats.group_bootstrap(groups, np.random.default_rng(3))
    g = np.asarray(groups)
    for j in np.unique(draw):  # every draw is one complete group
        members = rows[draw == j]
        assert len(set(g[members])) == 1
        assert sorted(members) == sorted(np.flatnonzero(g == g[members[0]]))
    assert len(np.unique(draw)) == 4  # G draws for G groups
    rows2, _ = stats.group_bootstrap(groups, np.random.default_rng(3))
    assert np.array_equal(rows, rows2)


# ---------------------------------------------------- permutation control
def _null_data(seed, n=36, p=3):
    rng = np.random.default_rng(seed)
    X = pd.DataFrame(rng.normal(size=(n, p)), columns=[f"x{i}" for i in range(p)])
    return X, rng.normal(size=n)


def _observed(X, y, families):
    return pl._nested_stat(ev.nested_cv(X, y, families, ad=None, **KW))["r2_pooled_mean"]


def test_permutation_control_null_and_signal():
    K = 9
    fams = ("linear",)
    # null: p roughly uniform on {0.1, ..., 1}
    ps = []
    for seed in range(6):
        X, y = _null_data(seed)
        perm = pl.permutation_control({"X": X}, y, families=fams, n_permutations=K,
                                      seed=seed, **KW)
        ps.append(stats.permutation_pvalue(_observed(X, y, fams), perm["r2_pooled_mean"])["p"])
    ps = np.array(ps)
    assert np.all((ps >= 1 / (K + 1)) & (ps <= 1))
    assert 0.25 <= ps.mean() <= 0.85, ps
    assert np.sum(ps <= 0.1) <= 2, ps
    # strong signal: every permutation is beaten, p = 1 / (K + 1), informative
    X, _ = _null_data(99)
    y = X.to_numpy() @ np.array([1.0, -2.0, 0.5]) + 0.05 * np.random.default_rng(1).normal(size=len(X))
    perm = pl.permutation_control({"X": X}, y, families=fams, n_permutations=K, seed=0, **KW)
    r = stats.permutation_pvalue(_observed(X, y, fams), perm["r2_pooled_mean"])
    assert r["p"] == pytest.approx(1 / (K + 1)) and r["uninformative"] is False
    assert r["observed"] > 0.9 and r["null_max"] < 0.5


def test_permutation_control_deterministic_shared_and_fold_preserving():
    X, y = _null_data(5)
    sets = {"A": X[["x0"]], "B": X}
    a = pl.permutation_control(sets, y, families=("linear",), n_permutations=3, seed=11, **KW)
    b = pl.permutation_control(sets, y, families=("linear",), n_permutations=3, seed=11, **KW)
    pd.testing.assert_frame_equal(a, b)
    c = pl.permutation_control(sets, y, families=("linear",), n_permutations=3, seed=12, **KW)
    assert not np.allclose(a["r2_pooled_mean"], c["r2_pooled_mean"])
    assert a["rng"].tolist()[:3] == [f"default_rng([11, {b}])" for b in range(3)]
    assert a["set"].tolist() == ["A"] * 3 + ["B"] * 3
    # permutation b is exactly default_rng([seed, b]).permutation(n), the same for every set
    yp = y[np.random.default_rng([11, 1]).permutation(y.size)]
    for nm in ("A", "B"):
        want = _observed(sets[nm], yp, ("linear",))
        got = a.loc[(a["set"] == nm) & (a["permutation"] == 1), "r2_pooled_mean"].iloc[0]
        assert got == pytest.approx(want, abs=0)
    # parallel permutations give identical rows
    d = pl.permutation_control(sets, y, families=("linear",), n_permutations=3, seed=11,
                               n_jobs=2, **KW)
    pd.testing.assert_frame_equal(a, d)


# -------------------------------------------------------------- bootstrap
def _effect_data(n=120, seed=0):
    rng = np.random.default_rng(seed)
    s, t, z = rng.normal(size=(3, n))
    y = s + t + rng.normal(size=n)
    idx = pd.Index([f"m{i:03d}" for i in range(n)], name="id")
    return {"ref": pd.DataFrame({"s": s}, index=idx),
            "A": pd.DataFrame({"s": s, "t": t}, index=idx),
            "N": pd.DataFrame({"s": s, "z": z}, index=idx)}, y


def test_bootstrap_ci_covers_known_effect_and_zero():
    sets, y = _effect_data()
    # population Delta R2 of A over ref = 1/3 - ... : R2_A = 2/3, R2_ref = 1/3
    rows, sm = pl.paired_bootstrap(sets, y, set_name="A", reference="ref", families=("linear",),
                                   n_boot=15, seed=0, **KW)
    assert sm["reason"] == "" and sm["n_valid"] == 15
    assert sm["d_r2_ci_low"] < 1 / 3 < sm["d_r2_ci_high"]
    assert sm["d_r2_ci_low"] > 0.1 and sm["d_mae_ci_low"] > 0
    # the original-data row is exactly the plain paired comparison
    cmp = ev.compare_feature_sets({"ref": sets["ref"], "A": sets["A"]}, pd.Series(y, index=sets["A"].index),
                                  reference="ref", families=("linear",), ad=None, **KW)
    assert sm["d_r2_original"] == pytest.approx(cmp.summary.set_index("set").loc["A", "d_r2_mean"])
    # a pure-noise column: the interval covers 0
    _, sn = pl.paired_bootstrap(sets, y, set_name="N", reference="ref", families=("linear",),
                                n_boot=15, seed=0, **KW)
    assert sn["d_r2_ci_low"] < 0 < sn["d_r2_ci_high"]
    assert list(rows["resample"]) == list(range(-1, 15))


def test_bootstrap_respects_units_and_is_deterministic():
    sets, y = _effect_data(n=45, seed=3)
    # pairs of rows form units (e.g. stereoisomers): drawn and folded together
    units = np.repeat([f"u{i}" for i in range(23)], 2)[:45]
    a, sa = pl.paired_bootstrap(sets, y, set_name="A", reference="ref", families=("linear",),
                                n_boot=4, seed=5, units=units, **KW)
    b, sb = pl.paired_bootstrap(sets, y, set_name="A", reference="ref", families=("linear",),
                                n_boot=4, seed=5, units=units, **KW)
    pd.testing.assert_frame_equal(a, b)
    assert sa == sb
    rows, _ = stats.group_bootstrap(units, np.random.default_rng([5, 0]))
    assert a.loc[a["resample"] == 0, "n_rows"].iloc[0] == rows.size
    # the resampled CV groups keep every copy of a unit in one outer fold
    g = units[rows]
    for _, _, _, tr, te in ev.make_outer_splits(rows.size, 3, [0], groups=g):
        assert not set(g[tr]) & set(g[te])


# ---------------------------------------------------------------- AD report
def test_ad_report_matches_predictions():
    preds = pd.DataFrame(dict(
        target="t", learner="ridge", set=["S"] * 6 + ["null"] * 2,
        repeat=0, fold=[0, 0, 0, 1, 1, 1, 0, 1], row=list("abcdefab"),
        pos=[0, 1, 2, 3, 4, 5, 0, 3],
        y_true=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 4.0],
        y_pred=[1.1, 2.1, 2.0, 4.0, 5.5, 6.5, 3.5, 3.5],
        in_domain=[1, 1, 0, 1, 1, 1, math.nan, math.nan]))
    folds = pd.DataFrame(dict(target="t", learner="ridge", set=["S", "S", "null", "null"],
                              repeat=0, fold=[0, 1, 0, 1],
                              ad_reason=["", "", "no feature columns", "no feature columns"]))
    ad, summ = pl.ad_report(preds, folds)
    s0 = ad[(ad["set"] == "S") & (ad["fold"] == 0)].iloc[0]
    assert s0["n_outside"] == 1 and s0["share_outside"] == pytest.approx(1 / 3)
    assert s0["r2_in_domain"] == pytest.approx(1 - (0.01 + 0.01) / 0.5)
    assert s0["mae_outside"] == pytest.approx(1.0)
    assert math.isnan(s0["r2_outside"])  # one outside row: R2 undefined, NaN (not 0)
    ss = summ.set_index("set")
    assert ss.loc["S", "share_outside_fold_mean"] == pytest.approx((1 / 3 + 0) / 2)
    assert ss.loc["S", "share_outside_pooled"] == pytest.approx(1 / 6)
    assert math.isnan(ss.loc["null", "share_outside_fold_mean"])
    assert "no feature columns" in ss.loc["null", "reason"]
    yt = preds.loc[preds["set"] == "S", "y_true"].to_numpy()
    yp = preds.loc[preds["set"] == "S", "y_pred"].to_numpy()
    dom = preds.loc[preds["set"] == "S", "in_domain"].to_numpy() == 1
    assert ss.loc["S", "r2_in_domain_pooled_mean"] == pytest.approx(ev._metrics(yt[dom], yp[dom])["r2"])


# ---------------------------------------------------------- confirmatory table
def test_confirmatory_table_bh_by_and_flags():
    summary = pd.DataFrame(dict(
        target=["a", "b", "c", "d"], learner="L", set="TI", n=[100, 100, 100, 30],
        underpowered=[False, False, False, True], underpowered_reason=["", "", "", "n=30 < min_n=50"],
        reason="", d_r2_mean=[0.08, 0.01, -0.02, 0.2], d_r2_ci_low=[0.02, -0.01, -0.05, -0.1],
        d_r2_ci_high=[0.14, 0.03, 0.01, 0.5], d_r2_p=[0.001, 0.3, 0.04, 0.2],
        d_mae_mean=0.1, d_mae_ci_low=0.0, d_mae_ci_high=0.2, d_mae_p=0.05))
    fam = [dict(target=t, set="TI", reference="size") for t in "abcd"] + \
          [dict(target="e", set="TI", reference="size")]
    tab = pl.confirmatory_table(fam, summary, None, learner="L", practical_threshold=0.05,
                                q=0.05, skipped={"e": "too few rows"})
    # d is underpowered: out of the family; e was not evaluated (n unknown): in the
    # family, untestable, counted in m
    assert tab["in_family"].tolist() == [True, True, True, False, True]
    assert (tab["family_m"] == 4).all()
    p = np.where(tab["in_family"], tab["d_r2_p"].to_numpy(dtype=float), np.nan)
    bh, _ = stats.benjamini_hochberg(p, 0.05, m=4)
    by, _ = stats.benjamini_yekutieli(p, 0.05, m=4)
    np.testing.assert_allclose(tab["d_r2_q_bh"].to_numpy(dtype=float), bh, equal_nan=True)
    np.testing.assert_allclose(tab["d_r2_q_by"].to_numpy(dtype=float), by, equal_nan=True)
    # c: 0.04 * 4 / 2 = 0.08 > 0.05 -> not rejected (m = 4: e counts, d does not); a rejected
    assert tab.loc[2, "d_r2_q_bh"] == pytest.approx(0.08)
    assert math.isnan(tab.loc[4, "d_r2_q_bh"]) and "too few rows" in tab.loc[4, "reason"]
    assert math.isnan(tab.loc[3, "d_r2_q_bh"])
    assert tab["informative_bh"].tolist() == [True, False, False, None, None]
    assert tab["practically_significant"].tolist() == [True, False, False, None, None]
    assert "underpowered" in tab.loc[3, "reason"]
    assert "not in confirmatory family" in tab.loc[3, "reason"]
    assert list(tab.columns) == pl.CONFIRMATORY_COLUMNS


def _summ(rows):
    d = pd.DataFrame(rows)
    d["learner"], d["set"], d["reason"] = "L", "TI", ""
    for c in ("d_r2_ci_low", "d_r2_ci_high", "d_mae_mean", "d_mae_ci_low", "d_mae_ci_high",
              "d_mae_p"):
        d[c] = 0.0
    d["underpowered_reason"] = np.where(d["underpowered"], "n < min_n=50", "")
    return d


def _fam(targets):
    return [dict(target=t, set="TI", reference="size") for t in targets]


def test_underpowered_rows_are_not_in_the_bh_family():
    # regression: an underpowered row (n < min_n) with a tiny p used to be BH-tested,
    # flagged informative, and to lower the other members' q
    s = _summ([dict(target="a", n=12, underpowered=True, d_r2_p=0.001, d_r2_mean=0.3),
               dict(target="b", n=200, underpowered=False, d_r2_p=0.04, d_r2_mean=0.06)])
    tab = pl.confirmatory_table(_fam("ab"), s, None, learner="L").set_index("target")
    assert not tab.loc["a", "in_family"] and tab.loc["b", "in_family"]
    assert (tab["family_m"] == 1).all()
    assert math.isnan(tab.loc["a", "d_r2_q_bh"]) and math.isnan(tab.loc["a", "d_r2_q_by"])
    assert tab.loc["a", "informative_bh"] is None and tab.loc["a", "informative_by"] is None
    assert tab.loc["a", "practically_significant"] is None
    assert "underpowered: not in confirmatory family" in tab.loc["a", "reason"]
    assert tab.loc["a", "d_r2_mean"] == pytest.approx(0.3)  # estimate kept (descriptive)
    assert tab.loc["b", "d_r2_q_bh"] == pytest.approx(0.04)
    s = _summ([dict(target="a", n=200, underpowered=False, d_r2_p=0.04, d_r2_mean=0.06),
               dict(target="b", n=200, underpowered=False, d_r2_p=0.5, d_r2_mean=0.0),
               dict(target="c", n=20, underpowered=True, d_r2_p=0.001, d_r2_mean=0.3)])
    tab = pl.confirmatory_table(_fam("abc"), s, None, learner="L").set_index("target")
    assert tab.loc["a", "d_r2_q_bh"] == pytest.approx(0.08)  # m = 2 (c excluded)
    assert tab.loc["c", "informative_bh"] is None
    # a target not evaluated with n < min_n is excluded; with n >= min_n it stays in m
    tab = pl.confirmatory_table(_fam("abxy"), s, None, learner="L", min_n=50,
                                skipped={"x": "too few CV units", "y": "too few CV units"},
                                n_values={"x": 4, "y": 60}).set_index("target")
    assert tab["in_family"].to_dict() == {"a": True, "b": True, "x": False, "y": True}
    assert tab.loc["x", "n"] == 4 and tab.loc["x", "underpowered"] is True
    assert (tab["family_m"] == 3).all()
    assert tab.loc["a", "d_r2_q_bh"] == pytest.approx(0.04 * 3)


def test_untestable_member_counts_in_m():
    # regression: a NaN p used to shrink m (anti-conservative)
    def run(last):
        s = _summ([dict(target=t, n=200, underpowered=False, d_r2_p=p, d_r2_mean=0.1)
                   for t, p in zip("abcd", [0.02, 0.03, 0.04, last])])
        return pl.confirmatory_table(_fam("abcd"), s, None, learner="L")
    nan, one = run(math.nan), run(1.0)
    np.testing.assert_allclose(nan["d_r2_q_bh"][:3], one["d_r2_q_bh"][:3])
    np.testing.assert_allclose(nan["d_r2_q_by"][:3], one["d_r2_q_by"][:3])
    np.testing.assert_allclose(nan["d_r2_q_bh"][:3], [0.04 * 4 / 3] * 3)
    assert nan["informative_bh"].tolist() == [False, False, False, None]
    assert (nan["family_m"] == 4).all() and nan["in_family"].all()
    assert math.isnan(nan.loc[3, "d_r2_q_bh"])
    assert "d_r2_p undefined" in nan.loc[3, "reason"]


def test_bh_by_fixed_family_size():
    p = np.array([0.01, np.nan, 0.03, 0.2])
    for f in (stats.benjamini_hochberg, stats.benjamini_yekutieli):
        a, ra = f(p, 0.05, m=4)
        b, rb = f(np.where(np.isnan(p), 1.0, p), 0.05)
        assert math.isnan(a[1]) and not ra[1]
        np.testing.assert_allclose(np.delete(a, 1), np.delete(b, 1))
        assert (np.delete(ra, 1) == np.delete(rb, 1)).all()
        # default m = number of finite p (unchanged behaviour)
        np.testing.assert_allclose(f(p, 0.05)[0], f(p, 0.05, m=3)[0], equal_nan=True)
        with pytest.raises(ValueError, match="smaller than"):
            f(p, 0.05, m=2)


def test_confirmatory_table_explains_missing_components():
    s = _summ([dict(target="a", n=200, underpowered=False, d_r2_p=0.04, d_r2_mean=0.06)])
    ps = pd.DataFrame(dict(target="a", learner="L", set=["TI"], n_permutations=5, p=0.5,
                           uninformative=True, reason="2 of 5 permutations failed (excluded)"))
    tab = pl.confirmatory_table(_fam("a"), s, None, learner="L", perm_summary=ps,
                                permutation_note="set not in protocol.permutation.feature_sets",
                                bootstrap_note="bootstrap disabled", ad_note="AD not fitted")
    why = tab.loc[0, "reason"]
    assert "permutation set: 2 of 5 permutations failed" in why
    assert "permutation reference (size): set not in protocol.permutation.feature_sets" in why
    assert "bootstrap disabled" in why and "AD set: AD not fitted" in why


# ------------------------------------------------- restricted permutation
def test_permutation_indices_row_wise_without_groups():
    a = pl.permutation_indices(20, None, np.random.default_rng([3, 1]))
    np.testing.assert_array_equal(a, np.random.default_rng([3, 1]).permutation(20))


def test_permutation_indices_moves_whole_blocks_among_equal_sizes():
    # pairs, triples, singletons and one group whose size (4) no other group shares
    groups = np.array(["p0", "p0", "p1", "p1", "p2", "p2", "t0", "t0", "t0", "t1", "t1", "t1",
                       "s0", "s1", "s2", "q", "q", "q", "q"], dtype=object)
    groups = groups[np.random.default_rng(0).permutation(groups.size)]  # non-contiguous groups
    codes = {g: float(i) for i, g in enumerate(sorted(set(groups)))}
    y = np.array([codes[g] for g in groups])  # y shared within a group
    rest = np.flatnonzero(np.isin(groups, ["q", "s0", "s1", "s2"]))
    seen = set()
    for b in range(30):
        perm = pl.permutation_indices(groups.size, groups, np.random.default_rng([7, b]))
        assert sorted(perm.tolist()) == list(range(groups.size))
        yp = y[perm]
        for g in set(groups):
            rows = np.flatnonzero(groups == g)
            if len(rows) in (2, 3):
                # a pair / triple receives the whole y block of one group of its own size
                src = set(groups[perm[rows]])
                assert len(src) == 1 and (groups == src.pop()).sum() == len(rows)
                assert len(set(yp[rows])) == 1
        # singletons and the unique-size group are permuted row-wise among themselves
        assert set(perm[rest]) == set(rest)
        seen.add(tuple(perm))
        np.testing.assert_array_equal(
            perm, pl.permutation_indices(groups.size, groups, np.random.default_rng([7, b])))
    assert len(seen) > 10  # it does permute


def test_permutation_control_uses_block_permutation_with_groups():
    X, y = _null_data(5)
    n = len(y)
    groups = np.array([f"g{i // 2}" for i in range(n)], dtype=object)
    perm = pl.permutation_control({"X": X}, y, families=("linear",), n_permutations=2, seed=4,
                                  groups=groups, **KW)
    idx = pl.permutation_indices(n, groups, np.random.default_rng([4, 1]))
    assert not np.array_equal(idx, np.random.default_rng([4, 1]).permutation(n))
    want = pl._nested_stat(ev.nested_cv(X, y[idx], ("linear",), groups=groups, ad=None, **KW))
    got = perm.loc[perm["permutation"] == 1, "r2_pooled_mean"].iloc[0]
    assert got == pytest.approx(want["r2_pooled_mean"], abs=0)


# ----------------------------------------------------- run_study integration
def _protocol_config(**proto):
    p = {"family": ["MolLogP:TI-size", "MolWt:TI-size", "MolLogP:TI_degree-size"],
         "learner": "selected", "practical_threshold": 0.05,
         "permutation": {"n_permutations": 3, "seed": 3},
         "bootstrap": {"n_boot": 2, "seed": 4}}
    p.update(proto)
    return {"name": "tiny-protocol", "cohort": {"source": "smiles_csv", "path": "tiny.csv"},
            "feature_sets": ["null", "size", "TI", "TI_family"],
            "targets": [{"provider": "rdkit", "names": ["MolLogP", "MolWt"]}],
            "indices": {"names": ["W", "M1", "M2", "R", "J"]},
            "evaluation": {"families": ["ridge"], "learners": ["selected"], "seeds": [0],
                           "outer_splits": 3, "inner_splits": 2, "min_n": 10, "min_n_test": 3},
            "protocol": p}


@pytest.fixture(scope="module")
def proto_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("proto") / "res"
    return pl.run_study(_protocol_config(), out, timestamp_utc=TS, base_dir=FIX)


def test_protocol_outputs_and_record(proto_run):
    out = proto_run.out_dir
    for f in pl.PROTOCOL_CSVS + pl.AD_CSVS:
        assert (out / f).exists(), f
    rec = proto_run.run_record
    assert set(pl.PROTOCOL_CSVS + pl.AD_CSVS) <= set(rec["result_csvs"])
    pr = rec["protocol"]
    assert pr["practical_threshold"] == 0.05 and pr["q"] == 0.05
    assert [(f["target"], f["set"], f["reference"]) for f in pr["family"]] == [
        ("MolLogP", "TI", "size"), ("MolWt", "TI", "size"), ("MolLogP", "TI_degree", "size")]
    assert pr["permutation"]["n_permutations"] == 3 and pr["permutation"]["seed"] == 3
    assert pr["bootstrap"]["n_boot"] == 2 and pr["bootstrap"]["families"] == ["ridge"]
    assert rec["config"]["protocol"]["family"][0] == {"target": "MolLogP", "set": "TI",
                                                      "reference": "size"}
    assert any(k.startswith("protocol:permutation") for k in rec["timings_s"])
    assert {t["name"]: t["target_class"] for t in rec["targets"]} == {"MolLogP": 2, "MolWt": 1}
    assert pr["family_membership"] == pl.FAMILY_RULE and pr["min_n"] == 10
    assert pr["permutation"]["enabled"] is True and pr["bootstrap"]["enabled"] is True
    assert pr["permutation"]["scheme"] == pl.PERMUTATION_SCHEME
    assert set(pr["permutation"]["scheme_by_target"]) == {"MolLogP", "MolWt"}
    assert set(pr["permutation"]["scheme_by_target"].values()) <= {"block", "row-wise"}


def test_protocol_disabled_components_are_explained(tmp_path):
    # regression: perm/boot columns were NaN with an empty reason, and the record had
    # no 'bootstrap' key when the bootstrap was off
    cfg = _protocol_config(family=["MolLogP:TI-size"],
                           permutation={"n_permutations": 2, "seed": 3, "feature_sets": ["TI"]},
                           bootstrap={"enabled": False})
    run = pl.run_study(cfg, tmp_path / "res", timestamp_utc=TS, base_dir=FIX)
    conf = pl.read_results(run.out_dir, "confirmatory")
    r = conf.iloc[0]
    assert math.isnan(r["perm_p_reference"]) and math.isnan(r["boot_d_r2_ci_low"])
    assert "permutation reference (size): set not in protocol.permutation.feature_sets" \
        in r["reason"]
    assert "bootstrap disabled (protocol.bootstrap.enabled = false)" in r["reason"]
    pr = run.run_record["protocol"]
    assert pr["bootstrap"] == {"enabled": False, "reason": "protocol.bootstrap.enabled = false"}
    assert pr["permutation"]["enabled"] is True


def test_protocol_confirmatory_consistent_with_components(proto_run):
    out = proto_run.out_dir
    conf = pl.read_results(out, "confirmatory")
    summ = pl.read_results(out, "summary").set_index(["target", "learner", "set"])
    assert len(conf) == 3 and (conf["comparison"] == ["TI-size", "TI-size", "TI_degree-size"]).all()
    for _, r in conf.iterrows():
        s = summ.loc[(r["target"], "selected", r["set"])]
        for c in ("d_r2_mean", "d_r2_ci_low", "d_r2_ci_high", "d_r2_p", "d_mae_mean", "d_mae_p"):
            assert r[c] == pytest.approx(s[c], abs=1e-12)
    bh, _ = stats.benjamini_hochberg(conf["d_r2_p"].to_numpy(float), 0.05)
    by, _ = stats.benjamini_yekutieli(conf["d_r2_p"].to_numpy(float), 0.05)
    np.testing.assert_allclose(conf["d_r2_q_bh"], bh)
    np.testing.assert_allclose(conf["d_r2_q_by"], by)
    assert (conf["practically_significant"] == (conf["d_r2_mean"] >= 0.05)).all()
    # permutation p recomputed from permutation.csv and the observed statistic
    perm = pl.read_results(out, "permutation")
    ps = pl.read_results(out, "permutation_summary").set_index(["target", "set"])
    assert len(perm) == 5 * 3  # MolLogP: TI, size, TI_degree; MolWt: TI, size
    for _, r in conf.iterrows():
        null = perm.loc[(perm["target"] == r["target"]) & (perm["set"] == r["set"]), "r2_pooled_mean"]
        assert len(null) == 3
        obs = ps.loc[(r["target"], r["set"]), "observed_r2_pooled_mean"]
        assert obs == pytest.approx(summ.loc[(r["target"], "selected", r["set"]), "r2_pooled_mean"])
        assert r["perm_p_set"] == pytest.approx((1 + (null >= obs).sum()) / 4)
        assert r["perm_k"] == 3
    # AD share = share of outer-test predictions flagged outside, per fold, averaged
    preds = pl.read_results(out, "predictions")
    for _, r in conf.iterrows():
        g = preds[(preds["target"] == r["target"]) & (preds["learner"] == "selected")
                  & (preds["set"] == r["set"])]
        want = g.groupby(["repeat", "fold"])["in_domain"].apply(lambda x: np.mean(x == 0)).mean()
        assert r["ad_share_outside_set"] == pytest.approx(want)
    boot = pl.read_results(out, "bootstrap")
    assert len(boot) == 3 * (2 + 1)
    for _, r in conf.iterrows():
        b = boot[(boot["target"] == r["target"]) & (boot["set"] == r["set"])]
        assert r["boot_d_r2_original"] == pytest.approx(r["d_r2_mean"])
        assert r["boot_d_r2_ci_low"] == pytest.approx(
            np.quantile(b.loc[b["resample"] >= 0, "d_r2"], 0.025))


def test_protocol_is_deterministic(proto_run, tmp_path):
    rep = pl.rerun(proto_run.out_dir / "run_record.json", tmp_path / "again", timestamp_utc=TS,
                   base_dir=FIX)
    assert rep["ok"], rep
    assert set(pl.PROTOCOL_CSVS) <= set(rep["identical"])


def test_protocol_auto_family_and_off_by_default(tmp_path):
    cfg = pl.resolve_config({"cohort": {"source": "smiles_csv", "path": "x.csv"}})
    assert "protocol" not in cfg  # old configs keep their hash
    cfg = pl.resolve_config({"cohort": {"source": "smiles_csv", "path": "x.csv"},
                             "protocol": {"family": "auto"}})
    assert cfg["protocol"]["permutation"]["n_permutations"] == 200  # EXPERIMENT_PLAN §7
    assert cfg["protocol"]["practical_threshold"] == 0.05
    from topoqspr.target_providers import TargetSpec
    specs = [TargetSpec("a", 3, None, "", "x"), TargetSpec("b", 2, None, "", "x")]
    assert pl._expand_family(cfg["protocol"], specs, "size") == [
        {"target": "a", "set": "TI", "reference": "size"}]


@pytest.mark.parametrize("proto, match", [
    ({"family": ["MolLogP:fingerprint-size"]}, "not an evaluated feature set"),
    ({"family": ["MolLogP:TI-null"]}, "differs from evaluation.reference"),
    ({"family": ["MolLogP"]}, "target:SET-REFERENCE"),
    ({"family": ["MolLogP:TI", "MolLogP:TI-size"]}, "duplicate"),
    ({"family": "everything"}, "null, 'auto'"),
    ({"learner": "rf"}, "not in evaluation.learners"),
    ({"permutation": {"n_permutations": 0}}, "n_permutations"),
    ({"bootstrap": {"families": ["lasso"]}}, "bootstrap.families"),
    ({"bootstrap": {"n_boot": 1}}, "n_boot"),
    ({"perm": {}}, "unknown key"),
])
def test_protocol_config_validation(proto, match):
    cfg = _protocol_config()
    cfg["protocol"] = {"family": ["MolLogP:TI"], **proto}
    with pytest.raises((ValueError, KeyError), match=match):
        pl.resolve_config(cfg)


def test_family_target_must_exist(tmp_path):
    cfg = _protocol_config(family=["logS:TI-size"])
    with pytest.raises(ValueError, match="no provider supplies"):
        pl.run_study(cfg, tmp_path / "x", timestamp_utc=TS, base_dir=FIX)
    assert not (tmp_path / "x").exists()
