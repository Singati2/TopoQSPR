"""Tests for features.py, models.py and evaluate.py (synthetic data, no network)."""

from __future__ import annotations

import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.ensemble import RandomForestRegressor
from sklearn.exceptions import NotFittedError

from topoqspr import evaluate as ev
from topoqspr import features as ft
from topoqspr import models as md

REPO = Path(__file__).resolve().parents[1]
FAST = ("ridge", "knn")


# ----------------------------------------------------------------- features
@pytest.mark.parametrize("smiles,n,m,rings", [
    ("c1ccccc1", 6, 6, 1),
    ("CCO", 3, 2, 0),
    ("C1CC2CCC1C2", 7, 8, 2),            # norbornane
    ("C12C3C4C1C5C2C3C45", 8, 12, 5),    # cubane: cyclomatic number 5
    ("[Na+].[Cl-]", 2, 0, 0),            # two components, no bonds
    ("CC(=O)[O-].[Na+]", 5, 3, 0),
])
def test_size_and_rings(smiles, n, m, rings):
    mol = Chem.MolFromSmiles(smiles)
    df = ft.size_rings_features([mol])
    assert df.loc[0, "n_heavy"] == n
    assert df.loc[0, "n_bonds"] == m
    assert df.loc[0, "n_rings"] == rings
    # fresh copies: GetSSSR overwrites the molecule's ring perception
    assert df.loc[0, "n_rings"] == len(Chem.GetSSSR(Chem.MolFromSmiles(smiles)))
    # RDKit's CalcNumRings counts the *symmetrised* SSSR: cubane 6, not 5
    if smiles.startswith("C12C3"):
        assert rdMolDescriptors.CalcNumRings(Chem.MolFromSmiles(smiles)) == 6


def test_explicit_hydrogens_do_not_change_graph():
    mol = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    a = ft.size_rings_features([mol])
    b = ft.size_rings_features([Chem.AddHs(mol)])
    pd.testing.assert_frame_equal(a, b)
    ca = ft.composition_features([mol])
    cb = ft.composition_features([Chem.AddHs(mol)])
    pd.testing.assert_frame_equal(ca, cb)


def test_composition_counts():
    mol = Chem.MolFromSmiles("OC(=O)c1ccc(Cl)cc1[Na]")  # odd but valid
    df = ft.composition_features([mol]).iloc[0]
    assert df["n_C"] == 7 and df["n_O"] == 2 and df["n_Cl"] == 1
    assert df["n_other"] == 1  # Na
    assert df["n_N"] == 0 and df["n_H"] == 4
    eth = ft.composition_features([Chem.MolFromSmiles("CCO")]).iloc[0]
    assert (eth["n_C"], eth["n_O"], eth["n_H"]) == (2, 1, 6)


def test_none_mol_is_nan_with_reason():
    df = ft.size_features([None, Chem.MolFromSmiles("CC")], index=["a", "b"])
    assert df.loc["a"].isna().all()
    assert df.loc["b", "n_heavy"] == 2
    assert {f["reason"] for f in df.attrs["failures"]} == {"mol is None"}
    assert {f["row"] for f in df.attrs["failures"]} == {"a"}


def test_exclusion_names_exist_in_rdkit():
    """Typos in the exclusion map would silently keep a leaking descriptor."""
    names = {n for n, _ in Descriptors.descList}
    for target, excl in ft.TARGET_EXCLUSIONS.items():
        missing = set(excl) - names
        assert not missing, (target, missing)


def test_rdkit2d_exclusions():
    keep = set(ft.rdkit2d_names("LogP"))
    for bad in ("MolLogP", "SlogP_VSA1", "SlogP_VSA12", "BCUT2D_LOGPHI", "MolMR",
                "SMR_VSA10", "qed"):
        assert bad not in keep
    assert "TPSA" in keep and "MolWt" in keep
    assert {"TPSA", "qed"}.isdisjoint(ft.rdkit2d_names("TPSA"))
    assert "LabuteASA" in ft.rdkit2d_names("TPSA")
    assert {"MolWt", "ExactMolWt", "HeavyAtomMolWt", "BCUT2D_MWHI"}.isdisjoint(
        ft.rdkit2d_names("MolWt"))
    assert "MolLogP" in ft.rdkit2d_names("XLogP3")
    assert len(ft.rdkit2d_names(None)) == len(Descriptors.descList)
    with pytest.raises(KeyError):
        ft.rdkit2d_names("SomethingNew")


def test_rdkit2d_values_and_build_sets():
    mols = [Chem.MolFromSmiles(s) for s in ("CCO", "c1ccccc1O", "CC(=O)O")]
    sets = ft.build_feature_sets(mols, "MolWt", index=["x", "y", "z"],
                                 ti=pd.DataFrame({"T1": [1.0, np.inf]}, index=["x", "y"]))
    assert list(sets) == ["null", "size", "size_rings", "composition", "rdkit2d", "TI"]
    assert sets["null"].shape == (3, 0)
    r = sets["rdkit2d"]
    assert "MolWt" not in r.columns
    assert r.loc["y", "TPSA"] == pytest.approx(20.23)
    ti = sets["TI"]
    assert np.isnan(ti.loc["y", "T1"]) and np.isnan(ti.loc["z", "T1"])
    reasons = {(f["row"], f["reason"]) for f in ti.attrs["failures"]}
    assert ("z", "row absent from TI frame") in reasons
    assert ("y", "non-finite value") in reasons


# ------------------------------------------------------------------- models
def test_drop_constant_and_duplicate():
    X = np.array([[1, 5, 1, 2], [2, 5, 2, 3], [3, 5, 3, 1.0]])
    t = md.DropConstantAndDuplicate().fit(X)
    assert list(t.keep_) == [0, 3]
    with pytest.raises(ValueError):
        md.DropConstantAndDuplicate().fit(np.ones((4, 2)))


def test_linear_family_feasibility():
    fam = md.get_family("linear")
    assert not fam.feasible(10, 9)[0]
    assert fam.feasible(11, 9)[0]
    pls = md.get_family("pls")
    assert [g["model__n_components"] for g in pls.param_grid(50, 3)] == [1, 2, 3]
    assert not md.get_family("ridge").feasible(50, 0)[0]
    assert md.get_family("mean").feasible(50, 0)[0]


def test_legacy_models_match_v0x_source():
    """LEGACY_MODELS reproduces ml.MODELS_AND_GRIDS (grids, classes, kwargs)."""
    src = (REPO / "ml.py").read_text()
    tree = ast.parse(src)
    node = next(n for n in tree.body if isinstance(n, ast.Assign)
                and getattr(n.targets[0], "id", "") == "MODELS_AND_GRIDS")
    legacy = {}
    for k, v in zip(node.value.keys, node.value.values):
        call, grid = v.elts
        kwargs = {kw.arg: ast.unparse(kw.value) for kw in call.keywords}
        legacy[ast.literal_eval(k)] = (call.func.id, kwargs, ast.literal_eval(grid))
    assert list(legacy) == list(md.LEGACY_MODELS)
    for name, (cls, kwargs, grid) in legacy.items():
        factory, mygrid = md.LEGACY_MODELS[name]
        est = factory()
        assert type(est).__name__ == cls, name
        assert mygrid == grid, name
        for kw, val in kwargs.items():
            if kw == "n_jobs":
                continue  # v0.x -1, here 1: identical seeded results
            expect = md.LEGACY_RANDOM_STATE if val == "RANDOM_STATE" else ast.literal_eval(val)
            assert est.get_params()[kw] == expect, (name, kw)


# ------------------------------------------------------------ nested CV core
def _noise(seed, n=60, p=10):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.normal(size=(n, p))), rng.normal(size=n)


def test_noise_nested_not_optimistic_but_legacy_is():
    light_legacy = {k: md.LEGACY_MODELS[k] for k in
                    ("LinearReg", "Ridge", "Lasso", "ElasticNet", "SVR", "RandomForest")}
    gaps, btrain, nested = [], [], []
    for seed in range(4):
        X, y = _noise(seed)
        res = ev.nested_cv(X, y, ("ridge", "knn", "svr"), seeds=(0, 1), ad=None)
        leg = ev.legacy_protocol(X, y, families=light_legacy)
        o = ev.optimism(leg, res)
        nested.append(o["r2_nested_pooled"])
        gaps.append(o["B_CV"])
        btrain.append(o["B_train"])
    assert max(nested) <= 0.05
    assert np.mean(nested) <= 0.0
    assert np.mean(gaps) > 0          # tuning + best-of-family on the test folds
    assert min(btrain) > 0 and np.mean(btrain) > 0.2   # in-sample R2 of the winner


def _size_like(seed, n=120):
    rng = np.random.default_rng(seed)
    nh = rng.integers(8, 45, size=n).astype(float)
    rings = rng.integers(0, 5, size=n).astype(float)
    nb = nh - 1 + rings
    size = pd.DataFrame({"n_heavy": nh, "n_bonds": nb})
    branch = rng.normal(size=n)  # a structural quantity independent of size
    ti = pd.DataFrame({"T_degsum": 2 * nb, "T_nm": nh + nb, "T_geo": np.sqrt(nh * nb),
                       "T_ratio": nb / nh, "T_branch": branch})
    return size, ti, branch, rng


def test_size_only_target_delta_ci_covers_zero():
    size, ti, _, rng = _size_like(0)
    y = 0.5 * size["n_heavy"] + 0.3 * size["n_bonds"] + rng.normal(scale=2.0, size=len(size))
    comp = ev.compare_feature_sets({"size": size, "TI": ti}, y, families=("ridge",),
                                   seeds=(0, 1, 2), ad=None)
    row = comp.summary.set_index("set").loc["TI"]
    assert row["d_r2_ci_low"] <= 0 <= row["d_r2_ci_high"]
    # paired: identical outer folds
    a = comp.results["size"].predictions[["repeat", "fold", "pos"]]
    b = comp.results["TI"].predictions[["repeat", "fold", "pos"]]
    pd.testing.assert_frame_equal(a, b)


def test_extra_information_is_detected():
    """Positive control: when y depends on a non-size TI the CI excludes 0."""
    size, ti, branch, rng = _size_like(1)
    y = 0.5 * size["n_heavy"] + 3.0 * branch + rng.normal(scale=1.0, size=len(size))
    comp = ev.compare_feature_sets({"size": size, "TI": ti, "null": pd.DataFrame(index=size.index)},
                                   y, families=("ridge",), seeds=(0, 1), ad=None)
    s = comp.summary.set_index("set")
    assert s.loc["TI", "d_r2_ci_low"] > 0
    assert s.loc["TI", "d_r2_q"] < 0.05
    assert s.loc["null", "d_r2_ci_high"] < 0
    assert set(comp.results["null"].folds["family"]) == {"mean"}


def test_null_model_is_slightly_negative():
    X, y = _noise(3, n=50)
    res = ev.nested_cv(pd.DataFrame(index=X.index), y, seeds=(0, 1))
    s = res.summary()
    assert s["family_counts"] == {"mean": 10}
    assert -0.2 < s["r2_pooled_mean"] < 0.0


# -------------------------------------------------------------- safety
def test_clone_safety_stored_model_unchanged():
    rng = np.random.default_rng(5)
    X = rng.normal(size=(60, 4))
    ya = X @ np.array([1.0, 2.0, 0.0, 0.0]) + rng.normal(scale=0.1, size=60)
    yb = -3 * X[:, 3] + rng.normal(scale=0.1, size=60)
    fams = ("linear", "ridge", "knn")
    sel_a = ev.select_and_refit(X, ya, fams, inner=3)
    pred_a = sel_a.estimator.predict(X).copy()
    sel_b = ev.select_and_refit(X, yb, fams, inner=3)
    assert sel_a.estimator is not sel_b.estimator
    np.testing.assert_array_equal(sel_a.estimator.predict(X), pred_a)
    assert np.corrcoef(pred_a, ya)[0, 1] > 0.99
    # registry prototypes are never fitted
    for name in fams:
        proto = md.get_family(name).build(0)
        with pytest.raises(NotFittedError):
            proto.predict(X)
    res = ev.nested_cv(X, ya, fams, seeds=(0,), keep_models=True, ad=None)
    assert len({id(m) for m in res.models}) == len(res.models) == 5


def test_determinism_and_parallel_identity():
    X, y = _noise(7, n=50)
    X.iloc[3, 2] = np.nan  # imputation path
    rf_small = md.ModelFamily(  # stochastic estimator, kept small for speed
        "rf_small", lambda rs: md._pipe(RandomForestRegressor(10, random_state=rs, n_jobs=1),
                                        scale=False),
        {"model__max_features": [1.0, 0.5]})
    kw = dict(families=("ridge", "knn", rf_small), seeds=(0, 1))
    a = ev.nested_cv(X, y, **kw)
    b = ev.nested_cv(X, y, **kw)
    c = ev.nested_cv(X, y, n_jobs=2, **kw)
    for r1, r2 in ((a, b), (a, c)):
        pd.testing.assert_frame_equal(r1.predictions, r2.predictions)
        pd.testing.assert_frame_equal(r1.folds, r2.folds)
        pd.testing.assert_frame_equal(r1.repeats, r2.repeats)
        pd.testing.assert_frame_equal(r1.inner, r2.inner)
    assert a.meta["n_missing_cells"] == 1


def test_outputs_shape_and_pooled_metrics():
    X, y = _noise(8, n=40)
    res = ev.nested_cv(X, y, FAST, seeds=(0, 1))
    assert len(res.folds) == 10 and len(res.repeats) == 2
    assert len(res.predictions) == 80
    for r, g in res.predictions.groupby("repeat"):
        assert sorted(g["pos"]) == list(range(40))  # each row once per repeat
        yt, yp = g["y_true"].to_numpy(), g["y_pred"].to_numpy()
        r2 = 1 - ((yt - yp) ** 2).sum() / ((yt - yt.mean()) ** 2).sum()
        assert res.repeats.loc[r, "r2"] == pytest.approx(r2)
    assert set(res.folds["status"]) == {"ok"}
    assert res.predictions["in_domain"].isin([0.0, 1.0]).all()
    assert all(json.loads(p) is not None for p in res.folds["params"])


def test_ols_skipped_when_p_too_large():
    X, y = _noise(9, n=30, p=25)
    res = ev.nested_cv(X, y, ("linear", "ridge"), seeds=(0,), ad=None)
    lin = res.inner[res.inner["family"] == "linear"]
    assert (lin["status"] == "skipped").all()
    assert set(res.folds["family"]) == {"ridge"}


def test_bad_inputs_raise():
    X, y = _noise(1, n=20)
    y2 = y.copy()
    y2[0] = np.nan
    with pytest.raises(ValueError, match="y contains"):
        ev.nested_cv(X, y2, FAST)
    X2 = X.copy()
    X2.iloc[0, 0] = np.inf
    with pytest.raises(ValueError, match="inf"):
        ev.nested_cv(X2, y, FAST)


# ------------------------------------------------------------ fold leakage
_SPY_LOG: list[tuple[str, frozenset]] = []


class _Spy(RegressorMixin, BaseEstimator):
    """Records which row ids (column 0) every fit/predict call sees."""

    def __init__(self, shrink=1.0):
        self.shrink = shrink

    def fit(self, X, y):
        _SPY_LOG.append(("fit", frozenset(np.asarray(X)[:, 0].astype(int).tolist())))
        self.mean_ = float(np.mean(y)) * self.shrink
        return self

    def predict(self, X):
        _SPY_LOG.append(("predict", frozenset(np.asarray(X)[:, 0].astype(int).tolist())))
        return np.full(len(X), self.mean_)


SPY = md.ModelFamily("spy", lambda rs: _Spy(), {"shrink": [1.0, 0.5]})


def _check_spy_log(splits, groups=None):
    log = list(_SPY_LOG)
    pos = 0
    for (_, _, _, tr, te) in splits:
        test, train = frozenset(te.tolist()), frozenset(tr.tolist())
        end = next(i for i in range(pos, len(log)) if log[i] == ("predict", test))
        seg = log[pos:end]
        assert seg, "no selection happened"
        for kind, ids in seg:
            assert not (ids & test), f"outer-test rows seen during {kind}"
            if groups is not None:
                assert not ({groups[i] for i in ids} & {groups[i] for i in test})
        inner_fits = [ids for kind, ids in seg[:-1] if kind == "fit"]
        assert inner_fits and all(ids < train for ids in inner_fits)
        if groups is not None:  # inner splits are grouped as well
            for (k1, fit_ids), (k2, pred_ids) in zip(seg[:-1], seg[1:-1]):
                if (k1, k2) == ("fit", "predict"):
                    assert not ({groups[i] for i in fit_ids} & {groups[i] for i in pred_ids})
        assert seg[-1] == ("fit", train)  # refit on the full outer-training part
        pos = end + 1
    assert pos == len(log)


def test_no_fold_leakage_instrumented():
    n = 40
    X = np.column_stack([np.arange(n), np.random.default_rng(0).normal(size=n)])
    y = np.random.default_rng(1).normal(size=n)
    _SPY_LOG.clear()
    ev.nested_cv(X, y, (SPY,), seeds=(0, 1), ad=None)
    _check_spy_log(ev.make_outer_splits(n, 5, (0, 1)))


def test_grouped_folds_keep_groups_together():
    n = 48
    groups = np.repeat(np.arange(16), 3)
    X = np.column_stack([np.arange(n), np.random.default_rng(0).normal(size=n)])
    y = np.random.default_rng(1).normal(size=n)
    _SPY_LOG.clear()
    res = ev.nested_cv(X, y, (SPY,), seeds=(3,), groups=groups, ad=None)
    splits = ev.make_outer_splits(n, 5, (3,), groups)
    for _, _, _, tr, te in splits:
        assert not set(groups[tr]) & set(groups[te])
    _check_spy_log(splits, groups)
    assert res.meta["grouped"]


def test_scaffold_groups():
    mols = [Chem.MolFromSmiles(s) for s in
            ("c1ccccc1C", "c1ccccc1CC", "c1ccncc1", "CCO", "CCCO")]
    g = ev.scaffold_groups(mols)
    assert g[0] == g[1] == "c1ccccc1"
    assert g[2] != g[0]
    assert g[3] != g[4]
    g2 = ev.scaffold_groups(mols, acyclic="shared")
    assert g2[3] == g2[4]


def test_legacy_protocol_reports_training_r2_and_drops_nonfinite():
    X, y = _noise(2, n=40)
    X.iloc[0, 0] = np.nan
    leg = ev.legacy_protocol(X, y, families={k: md.LEGACY_MODELS[k] for k in ("Ridge", "SVR")})
    assert leg.n_used == 39 and leg.dropped == [0]
    assert leg.best_family in ("Ridge", "SVR")
    assert leg.cv_r2 == pytest.approx(leg.table["R2_mean"].max())
    assert leg.train_r2 > leg.cv_r2


@pytest.mark.parametrize("name", md.DEFAULT_FAMILIES + ("mean",))
def test_every_family_fits_and_predicts(name):
    rng = np.random.default_rng(4)
    X = rng.normal(size=(40, 4))
    X[5, 1] = np.nan
    y = 100 + 50 * X[:, 0] + rng.normal(size=40)  # large-scale y (MolWt-like)
    sel = ev.select_and_refit(X, y, (name,), inner=3)
    row = sel.table.iloc[0]
    assert row["status"] == "ok", row["reason"]
    assert row["n_failed"] == 0
    pred = sel.estimator.predict(X)
    assert pred.shape == (40,) and np.all(np.isfinite(pred))
    if name not in ("mean", "knn"):
        assert np.corrcoef(pred, y)[0, 1] > 0.8


class _AlwaysFails(RegressorMixin, BaseEstimator):
    def __init__(self, a=1):
        self.a = a

    def fit(self, X, y):
        raise RuntimeError("boom")


def test_failing_family_is_reported_not_imputed():
    X, y = _noise(11, n=30)
    bad = md.ModelFamily("bad", lambda rs: _AlwaysFails(), {"a": [1, 2]})
    res = ev.nested_cv(X, y, (bad, "ridge"), seeds=(0,), ad=None)
    assert set(res.folds["family"]) == {"ridge"}
    b = res.inner[res.inner["family"] == "bad"]
    assert (b["status"] == "failed").all() and b["reason"].str.contains("failed").all()
    only = ev.nested_cv(X, y, (bad,), seeds=(0,), ad=None)
    assert (only.folds["status"] == "failed").all()
    assert only.predictions["y_pred"].isna().all()
    assert only.folds["r2"].isna().all()
