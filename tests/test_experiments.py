"""Tests for topoqspr.experiments and tools/run_experiment.py (offline).

Synthetic data where the answer is known; the only real data is the quarantined pilot
cohort (Migraine Disorders, D008881) from the recorded cache in
tests/fixtures/pipeline/chembl. Learners, seeds and folds are cut down to keep the
module quick on a loaded machine; no performance number of the pilot is asserted
beyond "defined".
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from topoqspr import experiments as ex
from topoqspr.cohort import MANIFEST_COLUMNS, Cohort
from topoqspr.models import LEGACY_MODELS

REPO = Path(__file__).resolve().parents[1]
FIX = REPO / "tests" / "fixtures"
SPEC = FIX / "experiments" / "synthetic_study.json"
TINY = FIX / "pipeline" / "tiny.csv"
PILOT_CACHE = FIX / "pipeline" / "chembl"
MIGRAINE = "D008881"
FAST = dict(seeds=[0], outer_splits=3, inner_splits=2, n_jobs=1)


def _tool():
    spec = importlib.util.spec_from_file_location("run_experiment", REPO / "tools" / "run_experiment.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _tiny():
    df = pd.read_csv(TINY, dtype=str, keep_default_na=False)
    return df[df["id"].str.startswith("m")].reset_index(drop=True)


# ------------------------------------------------------------------ guard
def test_quarantine_guard():
    ex.check_modelling_allowed([MIGRAINE])
    ex.check_modelling_allowed([" d008881 "])
    ex.check_modelling_allowed(None)
    with pytest.raises(ex.QuarantineError, match="quarantined"):
        ex.check_modelling_allowed(["D003924"])
    with pytest.raises(ex.QuarantineError):
        ex.check_modelling_allowed([MIGRAINE, "D006973"])
    ex.check_modelling_allowed(["D003924"], plan_frozen=True)


# ------------------------------------------------------------- structures
def test_prepare_structures_dedup_and_connectivity():
    st = ex.prepare_structures(
        ["a", "b", "c", "d", "e"],
        ["C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "CC(N)C(=O)O.Cl", "not_smiles", "CCCCO"])
    t = st.table.set_index("id")
    assert st.ids == ["a", "b", "c", "e"]  # c is the racemate: different full key
    assert t.loc["d", "status"] == "rejected"
    assert t.loc["a", "connectivity_group"] == "a" == t.loc["b", "connectivity_group"]
    assert t.loc["c", "connectivity_group"] == "a"
    assert t.loc["e", "connectivity_group"] == ""
    assert list(st.connectivity_labels()) == ["a", "a", "a", "e"]
    assert st.n_connectivity_groups == 1
    st2 = ex.prepare_structures(["x", "y"], ["CCO", "OCC"])
    assert st2.table.set_index("id").loc["y", "reason"] == "duplicate_inchikey:x"
    assert ex.cv_groups(st2, "random") is None
    g = ex.cv_groups(st, "scaffold")
    assert g[0] == g[1] == g[2]  # acyclic members share the first member's group


def test_feature_sets():
    df = _tiny()
    st = ex.prepare_structures(df["id"], df["smiles"])
    assert list(ex.feature_set(st, "size").columns) == ["n_heavy", "n_bonds"]
    fam = ex.feature_set(st, "TI_distance")
    assert "W" in fam.columns and "M1" not in fam.columns
    with pytest.raises(ValueError, match="exclusion"):
        ex.feature_set(st, "rdkit2d")
    rd = ex.feature_set(st, "rdkit2d", rdkit2d_exclusion="MolLogP")
    assert "MolLogP" not in rd.columns and "TPSA" in rd.columns
    with pytest.raises(ValueError):
        ex.feature_set(st, "nope")


# ---------------------------------------------------------- small helpers
def test_ccc_and_q2f3():
    y = np.array([1.0, 2.0, 3.0, 4.0])
    assert ex.ccc(y, y) == pytest.approx(1.0)
    assert ex.ccc(y, y + 1) < 1.0
    # Lin 1989: CCC = 2 s_xy / (s_x^2 + s_y^2 + (mx - my)^2)
    p = np.array([1.5, 1.5, 3.5, 3.0])
    num = 2 * np.mean((y - y.mean()) * (p - p.mean()))
    assert ex.ccc(y, p) == pytest.approx(num / (y.var() + p.var() + (y.mean() - p.mean()) ** 2))
    train = np.array([0.0, 2.0, 4.0])  # variance 8/3
    assert ex.q2_f3([1.0, 3.0], [2.0, 3.0], [train]) == pytest.approx(1 - 0.5 / (8 / 3))
    assert math.isnan(ex.q2_f3([1.0], [1.0], [np.array([2.0, 2.0])]))


def test_single_index_polynomial_exact():
    rng = np.random.default_rng(0)
    X = rng.uniform(-3, 3, size=(30, 3))
    y = 1.0 - 2.0 * X[:, 1] + 0.5 * X[:, 1] ** 3
    est = ex.SingleIndexPolynomial(column=1, degree=3).fit(X, y)
    assert np.allclose(est.predict(X), y, atol=1e-9)
    with pytest.raises(ValueError):
        ex.SingleIndexPolynomial(column=0, degree=3).fit(X[:4], y[:4])
    fam = ex.single_index_family(3, (1, 3))
    assert len(fam.param_grid(20, 3)) == 6
    assert len(fam.param_grid(4, 3)) == 3  # cubic needs >= 5 training rows


def test_map_index_name():
    assert ex.map_index_name("M1") == ("M1", "registry name")
    assert ex.map_index_name("first Zagreb index")[0] == "M1"
    assert ex.map_index_name("Randić")[0] == "R"
    assert ex.map_index_name("ABC_4") == ("S_ABC", "literature name in registry")
    assert ex.map_index_name("ND_3")[0] == "S_ReZ3"
    assert ex.map_index_name("wiener")[0] == "W"
    assert ex.map_index_name("sz_e")[0] == "Sz_e"
    name, why = ex.map_index_name("XYZ index")
    assert name is None and "no registry" in why
    assert ex.map_index_name("first Gourava")[0] is None


# ------------------------------------------------------------ Experiment A
def _manifest(rows, kind):
    out = []
    for mid, parent, smi in rows:
        key = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
        out.append(dict(disease_label="synthetic", mesh_ids=MIGRAINE, chembl_release="ChEMBL_37",
                        molecule_chembl_id=mid, parent_chembl_id=parent, pref_name=mid,
                        canonical_smiles=smi, standard_inchi_key=key,
                        molecule_type="Small molecule", phase_for_indication=4.0,
                        indication_ref_types="", efo_ids="", decision="include",
                        reason="included"))
    man = pd.DataFrame(out, columns=list(MANIFEST_COLUMNS))
    return Cohort(kind, "synthetic", "ChEMBL_37", "", {"mesh_ids": [MIGRAINE]}, man)


NAIVE = [("N1", "P1", "CCCCCC"), ("N2", "P2", "Cc1ccccc1"), ("N3", "P3", "CCN"),
         ("N4", "P3", "CCN.Cl"), ("N5", "P4", "O=C(O)c1ccccc1"),
         ("N6", "P4", "O=C([O-])c1ccccc1.[Na+]"), ("N7", "P5", "CCCCCCCC"),
         ("N8", "P6", "Oc1ccccc1")]
CONTROLLED = [("P2", "P2", "Cc1ccccc1"), ("P3", "P3", "CCN"), ("P4", "P4", "O=C(O)c1ccccc1"),
              ("P7", "P7", "c1ccc2ccccc2c1"), ("P8", "P8", "C1CCCCC1")]


def test_experiment_a_overlap_and_tables():
    res = ex.experiment_a(_manifest(NAIVE, "naive"), _manifest(CONTROLLED, "controlled"),
                          target_names=["MolWt", "MolLogP"], index_names=["M1", "M2", "W", "R"],
                          top_k=2)
    s = res.summary
    # standardised parents: naive 6 (2 salts collapse), controlled 5, 3 shared
    assert s["counts"]["naive"]["standardised_ok"] == 6
    assert s["counts"]["naive"]["rejected"] == {"duplicate_inchikey": 2}
    assert s["counts"]["controlled"]["standardised_ok"] == 5
    assert s["jaccard"]["standardised_parent_inchikey"] == pytest.approx(3 / 8)
    assert s["jaccard"]["structures_used"] == pytest.approx(3 / 10)  # salts are other keys
    assert s["jaccard"]["parent_chembl_id"] == pytest.approx(3 / 8)
    assert s["intersection_standardised"] == 3
    assert s["controlled_in_naive_share"] == pytest.approx(3 / 5)
    dist = res.tables["distributions"]
    assert list(dist["variable"]) == ["n_heavy", "n_bonds", "n_rings", "MolWt", "MolLogP"]
    assert dist["ks_q"].notna().all() and (dist["ks_q"] >= dist["ks_p"] - 1e-12).all()
    cor = res.tables["correlations"]
    assert len(cor) == 2 * 2 * 4 and set(cor["cohort"]) == {"naive", "controlled"}
    ov = res.tables["top_overlap"]
    assert set(ov["variant"]) == {"all", "bh_significant"}
    assert (ov.loc[ov["variant"] == "all", "n_naive"] == 2).all()
    assert (ov["n_overlap"] <= ov[["n_naive", "n_controlled"]].min(axis=1)).all()
    mem = res.tables["membership"]
    assert int((mem["in_naive"] & mem["in_controlled"]).sum()) == 3
    assert "model" not in res.tables


def test_experiment_a_top_k_order_is_deterministic():
    cor = pd.DataFrame(dict(index=list("abcd"), rho=[0.5, -0.9, 0.5, np.nan],
                            q=[0.01, 0.2, 0.01, np.nan]))
    assert ex._top_k(cor, 3) == ["b", "a", "c"]
    assert ex._top_k(cor, 3, alpha=0.05) == ["a", "c"]


def test_experiment_a_model_guard_and_pilot_model():
    naive, ctl = _manifest(NAIVE, "naive"), _manifest(CONTROLLED, "controlled")
    with pytest.raises(ex.QuarantineError):
        ex.experiment_a(naive, ctl, model=True, mesh_ids=["D003924"])
    df = _tiny()
    rows = [(i, i, s) for i, s in zip(df["id"], df["smiles"])]
    res = ex.experiment_a(_manifest(rows, "naive"), _manifest(rows[:20], "controlled"),
                          target_names=["MolLogP"], index_names=["M1", "W", "R", "J"],
                          model=True, families=("ridge",), min_n=10, min_n_test=3, **FAST)
    m = res.tables["model"]
    assert list(m["cohort"]) == ["naive", "controlled"]
    assert list(m["n"]) == [30, 20]
    assert m["d_r2_mean"].notna().all() and not m["underpowered"].any()


# ------------------------------------------------------------ Experiment C
def test_experiment_c_pure_noise():
    rng = np.random.default_rng(3)
    n = 60
    X = pd.DataFrame(rng.standard_normal((n, 15)), index=[f"r{i}" for i in range(n)])
    y = rng.standard_normal(n)
    legacy = {k: LEGACY_MODELS[k] for k in ("LinearReg", "RandomForest")}
    res = ex.experiment_c(X, y, scaffolds=[f"s{i % 15}" for i in range(n)],
                          families=("ridge", "knn"), legacy_families=legacy, seeds=[0],
                          outer_splits=5, inner_splits=3)
    s = res.summary
    assert s["r2_train"] > 0.6  # in-sample fit of the best v0.x family on pure noise
    assert s["r2_nested"] < 0.15 and s["r2_scaffold"] < 0.15
    assert s["B_train"] > 0.5
    assert math.isfinite(s["B_scaffold"]) and math.isfinite(s["B_CV"])
    assert s["scaffold"] == dict(n_groups=15, largest_group=4, largest_group_share=4 / 60,
                                 n_singletons=0)
    opt = res.tables["optimism"].set_index("term")
    assert opt.loc["B_train", "ci_low"] <= s["B_train"] <= opt.loc["B_train", "ci_high"]
    assert opt.loc["B_scaffold", "ci_low"] <= s["B_scaffold"] <= opt.loc["B_scaffold", "ci_high"]
    assert list(res.tables["chain"]["step"]) == ["train", "train_selected", "cv_legacy",
                                                 "nested_random", "nested_scaffold"]
    assert s["random_split"] == "KFold(shuffle=True)"


def test_experiment_c_signal_and_missing_and_groups():
    rng = np.random.default_rng(1)
    n = 50
    X = pd.DataFrame(rng.standard_normal((n, 4)))
    y = 2.0 * X[0].to_numpy() + 0.1 * rng.standard_normal(n)
    y[[3, 7]] = np.nan
    conn = [f"g{i}" if i >= 10 else "pair" if i < 2 else f"g{i}" for i in range(n)]
    res = ex.experiment_c(X, y, connectivity=conn, scaffolds=None, families=("ridge",),
                          legacy_families={"Ridge": LEGACY_MODELS["Ridge"]}, **FAST)
    s = res.summary
    assert s["n"] == 48 and s["n_target_missing"] == 2
    assert s["r2_nested"] > 0.9 and abs(s["B_CV"]) < 0.1
    assert s["connectivity_groups"] == 1
    assert s["random_split"].startswith("balanced_group_kfold")
    assert math.isnan(s["r2_scaffold"]) and "no scaffolds" in s["scaffold_reason"]


# ------------------------------------------------------------ Experiment D
@pytest.fixture(scope="module")
def d_result():
    spec = ex.load_study_spec(SPEC)
    return ex.experiment_d(spec, families=("ridge",), n_random=100, n_permutations=100, **FAST)


def test_experiment_d_reproduces_exact_polynomial(d_result):
    fits = d_result.tables["fits"].set_index(["property", "index", "form"])
    bp = fits.loc[("BP", "M1", "quadratic")]
    assert bp["r2"] == pytest.approx(1.0, abs=1e-12)
    assert np.allclose(json.loads(bp["coefficients"]), [12.5, 3.2, -0.015], atol=1e-8)
    assert bp["r2_minus_reported"] == pytest.approx(0.0, abs=1e-12)
    assert bp["loo_q2"] == pytest.approx(1.0, abs=1e-10)
    assert bp["loo_ccc"] == pytest.approx(1.0, abs=1e-10)
    mv = fits.loc[("MV", "R", "linear")]
    assert mv["r2"] == pytest.approx(1.0, abs=1e-12) and mv["reported_r2"] == 1.0
    assert fits.loc[("BP", "M1", "linear"), "r2"] < 1.0
    assert len(d_result.tables["fits"]) == 2 * 4 * 2  # properties x mapped indices x forms


def test_experiment_d_mapping_compounds_and_index_check(d_result):
    s = d_result.summary
    assert s["n_compounds_spec"] == 21 and s["n_compounds_used"] == 20
    assert s["unmappable"] == ["XYZ index"] and s["n_mapped"] == 4
    comp = d_result.tables["compounds"].set_index("name")
    assert comp.loc["unresolved drug", "status"] == "unresolved"
    chk = d_result.tables["index_check"].set_index(["compound", "study_label"])
    assert chk.loc[("hexane", "first Zagreb"), "agree"]
    assert chk.loc[("hexane", "Randic"), "agree"]  # printed to 4 decimals
    assert not chk.loc[("octane", "first Zagreb"), "agree"]  # wrong printed value
    assert s["index_values_agree"] == 3


def test_experiment_d_nested_and_controls(d_result):
    nest = d_result.tables["nested"]
    g = nest[(nest["learner"] == "genre") & (nest["split"] == "random")].set_index(["property", "set"])
    assert g.loc[("BP", "TI_study"), "r2_fold_mean"] == pytest.approx(1.0, abs=1e-8)
    assert g.loc[("BP", "TI_study"), "q2_f3"] == pytest.approx(1.0, abs=1e-8)
    assert g.loc[("BP", "TI_study"), "family_counts"] == '{"single_index": 3}'
    assert g.loc[("BP", "null"), "family_counts"] == '{"mean": 3}'
    assert math.isfinite(g.loc[("BP", "TI_study"), "d_r2_vs_size"])
    assert set(nest["learner"]) == {"genre", "families"}
    assert ((nest["split"] == "scaffold").sum()) == 2
    ctl = d_result.tables["controls"].set_index("property")
    assert ctl.loc["BP", "observed_best_r2"] == pytest.approx(1.0, abs=1e-12)
    assert ctl.loc["BP", "observed_best"] == "M1:quadratic"
    # 100 draws: the exact fit beats every pseudo-descriptor and every permutation
    assert ctl.loc["BP", "pseudo_p"] == pytest.approx(1 / 101)
    assert ctl.loc["BP", "perm_p"] == pytest.approx(1 / 101)
    assert ctl.loc["BP", "pseudo_null_p95"] < 0.8
    assert ctl.loc["BP", "k"] == 4


def test_experiment_d_spec_validation():
    spec = ex.load_study_spec(SPEC)
    with pytest.raises(ValueError, match="unknown study-spec keys"):
        ex.experiment_d(dict(spec, bogus=1))
    with pytest.raises(ValueError, match="model forms"):
        ex.experiment_d(dict(spec, model_forms=["exponential"]))


def test_experiment_d_deterministic(tmp_path):
    spec = ex.load_study_spec(SPEC)
    spec = dict(spec, properties=["MV"])
    kw = dict(families=("ridge",), n_random=50, n_permutations=50, scaffold=False, **FAST)
    a = ex.experiment_d(spec, **kw)
    b = ex.experiment_d(spec, **kw)
    a.write(tmp_path / "a")
    b.write(tmp_path / "b")
    for p in sorted((tmp_path / "a").iterdir()):
        assert p.read_bytes() == (tmp_path / "b" / p.name).read_bytes(), p.name


# --------------------------------------------------------- size dependence
def test_size_dependence_known_answers():
    df = _tiny()
    st = ex.prepare_structures(df["id"], df["smiles"])
    size = ex.size_table(st)
    ti, _ = ex.index_table(st, ["M1", "W"])
    rng = np.random.default_rng(0)
    desc = ti.assign(n_heavy_copy=size["n_heavy"].astype(float),
                     noise=rng.standard_normal(len(ti)), constant=1.0)
    res = ex.size_dependence(desc, size, families=("linear",),
                             family_of={"n_heavy_copy": "size", "noise": "other",
                                        "constant": "other"},
                             min_n=10, min_n_test=3, **FAST)
    per = res.tables["per_descriptor"].set_index(["descriptor", "set"])
    for s in ("size", "size_rings"):
        assert per.loc[("n_heavy_copy", s), "r2_fold_mean"] == pytest.approx(1.0, abs=1e-12)
        assert per.loc[("noise", s), "r2_fold_mean"] < 0.3
        assert per.loc[("constant", s), "reason"] == "constant descriptor"
    assert per.loc[("M1", "size"), "r2_fold_mean"] > 0.8
    fam = res.tables["per_family"]
    assert {"degree", "distance", "size", "other"} <= set(fam["family"])
    row = fam[(fam["family"] == "size") & (fam["set"] == "size")].iloc[0]
    assert row["n_r2_ge_0_99"] == 1


def test_poly2_family_fits_quadratic_size():
    rng = np.random.default_rng(0)
    n = 40
    size = pd.DataFrame({"n_heavy": rng.integers(5, 40, n).astype(float)})
    size["n_bonds"] = size["n_heavy"] + rng.integers(0, 4, n)
    size["n_rings"] = size["n_bonds"] - size["n_heavy"] + 1
    desc = pd.DataFrame({"quad": size["n_heavy"] ** 2}, index=size.index)
    res = ex.size_dependence(desc, size, families=(ex.POLY2_FAMILY,), min_n=10, min_n_test=3,
                             **FAST)
    assert (res.tables["per_descriptor"]["r2_fold_mean"] > 0.999).all()


# ------------------------------------------------------------ pilot smoke
def _pilot_structures():
    from topoqspr.chembl import ChemblClient
    from topoqspr.cohort import build_cohort

    coh = build_cohort([MIGRAINE], ChemblClient(PILOT_CACHE, offline=True))
    inc = coh.included
    return ex.prepare_structures(inc["parent_chembl_id"], inc["canonical_smiles"])


def test_pilot_smoke_c_and_size_dependence():
    st = _pilot_structures()
    assert len(st.ids) == 33
    y, _, _ = ex.compute_targets(st, "rdkit", ["MolLogP"])
    X = ex.feature_set(st, "TI")
    res = ex.experiment_c(X, y["MolLogP"], connectivity=st.connectivity_labels(), mols=st.mols,
                          families=("ridge",), legacy_families={"Ridge": LEGACY_MODELS["Ridge"]},
                          min_n=10, min_n_test=3, **FAST)
    s = res.summary
    for k in ("r2_train", "r2_cv", "r2_nested", "r2_scaffold", "B_train", "B_CV", "B_scaffold"):
        assert math.isfinite(s[k]), k
    assert s["scaffold"]["n_groups"] >= 3
    ti, _ = ex.index_table(st, ["M1", "W", "J"])
    sd = ex.size_dependence(ti, ex.size_table(st), families=("linear",),
                            groups=ex.cv_groups(st, "random"), min_n=10, min_n_test=3, **FAST)
    assert sd.tables["per_descriptor"]["r2_fold_mean"].notna().all()


# --------------------------------------------------------------------- CLI
def test_cli_exp_d_and_size_dependence(tmp_path):
    tool = _tool()
    out = tmp_path / "d"
    rc = tool.main(["exp-d", "--spec", str(SPEC), "--out", str(out), "--families", "ridge",
                    "--seeds", "0", "--outer-splits", "3", "--inner-splits", "2",
                    "--n-random", "20", "--n-permutations", "20", "--no-scaffold",
                    "--timestamp", "2026-10-03T00:00:00Z"])
    assert rc == 0
    names = sorted(p.name for p in out.iterdir())
    assert names == ["compounds.csv", "controls.csv", "fits.csv", "index_check.csv",
                     "index_map.csv", "index_values.csv", "nested.csv", "run_info.json",
                     "summary.json"]
    info = json.loads((out / "run_info.json").read_text())
    assert info["command"] == "exp-d" and len(info["spec"]["sha256"]) == 64
    assert not (tmp_path / ".d.partial").exists()
    # non-empty output directory is refused
    assert tool.main(["exp-d", "--spec", str(SPEC), "--out", str(out)]) == 2
    out2 = tmp_path / "sd"
    rc = tool.main(["size-dependence", "--smiles-csv", str(TINY), "--out", str(out2),
                    "--index", "M1", "W", "--families", "linear", ex.POLY2_FAMILY.name,
                    "--seeds", "0", "--outer-splits", "3", "--inner-splits", "2"])
    assert rc == 0
    per = pd.read_csv(out2 / "per_descriptor.csv", keep_default_na=False, na_values=[""])
    assert len(per) == 4
    mol = pd.read_csv(out2 / "molecules.csv", keep_default_na=False, na_values=[""])
    assert (mol["status"] == "rejected").sum() == 3  # bad SMILES, salt, duplicate


def test_cli_quarantine_refuses_before_network(tmp_path, capsys):
    tool = _tool()
    empty_cache = tmp_path / "cache"
    rc = tool.main(["exp-c", "--mesh", "D003924", "--target", "MolLogP", "--out",
                    str(tmp_path / "c"), "--offline", "--cache-dir", str(empty_cache)])
    assert rc == 2 and "quarantined" in capsys.readouterr().err
    rc = tool.main(["exp-a", "--mesh", "D006973", "--naive-text", "hypertension", "--model",
                    "--out", str(tmp_path / "a"), "--offline", "--cache-dir", str(empty_cache)])
    assert rc == 2 and "quarantined" in capsys.readouterr().err
    rc = tool.main(["size-dependence", "--mesh", "D001249", "--out", str(tmp_path / "s"),
                    "--offline", "--cache-dir", str(empty_cache)])
    assert rc == 2
    assert not empty_cache.exists()
    assert not (tmp_path / "c").exists() or not any((tmp_path / "c").iterdir())


def test_cli_exp_c_pilot_offline(tmp_path):
    tool = _tool()
    out = tmp_path / "c"
    rc = tool.main(["exp-c", "--mesh", MIGRAINE, "--target", "TPSA", "--feature-set", "size",
                    "--out", str(out), "--offline", "--cache-dir", str(PILOT_CACHE),
                    "--families", "ridge", "--legacy-families", "Ridge", "LinearReg",
                    "--seeds", "0", "--outer-splits", "3", "--inner-splits", "2"])
    assert rc == 0
    s = json.loads((out / "summary.json").read_text())
    assert s["n"] == 33 and s["experiment"] == "C"
    info = json.loads((out / "run_info.json").read_text())
    assert info["cohort"]["release"] == "ChEMBL_37" and "CC BY-SA" in info["data_licence"]
