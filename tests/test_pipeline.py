"""Tests for pipeline.py, provenance.py, target_providers.py and cli.py (offline).

Real data: only the quarantined pilot cohort (Migraine Disorders, D008881) from the
recorded ChEMBL_37 cache in tests/fixtures/pipeline/chembl; everything else uses the
synthetic SMILES file tests/fixtures/pipeline/tiny.csv. Learners and folds are cut
down (ridge, one repeat, three outer folds) so the module runs in about a minute.
"""

from __future__ import annotations

import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from topoqspr import cli
from topoqspr import pipeline as pl
from topoqspr import provenance as pv
from topoqspr import target_providers as tp
from topoqspr.chembl import ChemblClient

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
FIX = REPO / "tests" / "fixtures" / "pipeline"
CHEMBL = FIX / "chembl"
TINY = FIX / "tiny.csv"
TS = "2026-10-03T00:00:00Z"
MIGRAINE = "D008881"
FAST_EVAL = {"families": ["ridge"], "learners": ["selected"], "seeds": [0], "outer_splits": 3,
             "inner_splits": 2, "min_n": 10, "min_n_test": 3}


def _tiny_config(**over):
    cfg = {"name": "tiny", "cohort": {"source": "smiles_csv", "path": "tiny.csv"},
           "feature_sets": ["null", "size", "TI", "TI_family", "rdkit2d"],
           "targets": [{"provider": "rdkit", "names": ["MolLogP", "MolWt"]}],
           "evaluation": dict(FAST_EVAL)}
    cfg.update(over)
    return cfg


def _pilot_config():
    return {"name": "pilot-fast",
            "cohort": {"source": "chembl", "mesh_ids": [MIGRAINE], "chembl_release": "ChEMBL_37"},
            "feature_sets": ["null", "size", "composition", "TI", "rdkit2d"],
            "targets": [{"provider": "rdkit", "names": ["TPSA"]}],
            "evaluation": dict(FAST_EVAL),
            "runtime": {"chembl_cache_dir": str(CHEMBL), "offline": True}}


@pytest.fixture(scope="module")
def tiny_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("tiny") / "res"
    return pl.run_study(_tiny_config(), out, timestamp_utc=TS, base_dir=FIX)


@pytest.fixture(scope="module")
def pilot_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("pilot") / "res"
    return pl.run_study(_pilot_config(), out, timestamp_utc=TS)


# ------------------------------------------------------------------ config
def test_resolve_fills_defaults_and_is_idempotent():
    cfg = pl.resolve_config({"cohort": {"source": "chembl", "mesh_ids": "d008881"}})
    assert cfg["cohort"]["mesh_ids"] == [MIGRAINE]
    assert cfg["cohort"]["min_phase"] == 4.0
    assert cfg["evaluation"]["seeds"] == [0, 1, 2, 3, 4]
    assert cfg["evaluation"]["learners"][0] == "selected"
    assert cfg["evaluation"]["learners"][1:] == cfg["evaluation"]["families"]
    assert [t["names"] for t in cfg["targets"]] == [["MolWt", "MolLogP", "TPSA", "MolMR"]]
    assert pl.resolve_config(cfg) == cfg


@pytest.mark.parametrize("bad, match", [
    ({"cohort": {"source": "chembl", "mesh_ids": ["D008881"]}, "seed": 1}, "unknown key"),
    ({"cohort": {"source": "chembl", "mesh_id": ["D008881"]}}, "unknown key"),
    ({"cohort": {"source": "pubchem"}}, "cohort.source"),
    ({"cohort": {"source": "smiles_csv"}}, "cohort.path"),
    ({"cohort": {"source": "chembl", "mesh_ids": ["D008881"]},
      "feature_sets": ["null", "TI"]}, "reference"),
    ({"cohort": {"source": "chembl", "mesh_ids": ["D008881"]},
      "evaluation": {"learners": ["selected", "lasso"]}}, "unknown learners"),
    ({"cohort": {"source": "chembl", "mesh_ids": ["D008881"]},
      "feature_sets": ["size", "fingerprints"]}, "unknown feature sets"),
])
def test_resolve_rejects_bad_configs(bad, match):
    with pytest.raises((ValueError, KeyError), match=match):
        pl.resolve_config(bad)


def test_config_hash_stability():
    a = pl.resolve_config({"cohort": {"source": "chembl", "mesh_ids": [MIGRAINE]}})
    # same study, keys in another order, defaults written out, other runtime
    b_raw = json.loads(json.dumps(a))
    b_raw = {k: b_raw[k] for k in reversed(list(b_raw))}
    b_raw["runtime"] = {"n_jobs": 8, "chembl_cache_dir": "/elsewhere", "offline": True}
    b = pl.resolve_config(b_raw)
    assert pv.config_sha256(a) == pv.config_sha256(b)
    c = pl.resolve_config({"cohort": {"source": "chembl", "mesh_ids": [MIGRAINE]},
                           "evaluation": {"seeds": [0, 1, 2, 3, 5]}})
    assert pv.config_sha256(c) != pv.config_sha256(a)
    # golden value for the shipped pilot example: changes only with a deliberate
    # change of the config schema or its defaults
    ex = pl.resolve_config(pl.load_config(REPO / "examples" / "pilot_migraine.json"))
    assert pv.config_sha256(ex) == PILOT_CONFIG_SHA256


PILOT_CONFIG_SHA256 = "f9478a2c5b12d846aa0486bac345c1d05904ee3c705730a1760c5a5971e2cdea"


def test_config_hash_independent_of_hash_seed():
    code = ("import sys; sys.path.insert(0, %r); from topoqspr import pipeline as p, "
            "provenance as v; print(v.config_sha256(p.resolve_config(p.load_config(%r))))"
            % (str(SRC), str(REPO / "examples" / "pilot_migraine.json")))
    outs = set()
    for seed in ("0", "12345"):
        env = {**os.environ, "PYTHONHASHSEED": seed}
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                           env=env, check=True)
        outs.add(r.stdout.strip())
    assert len(outs) == 1


def test_canonical_json():
    assert pv.canonical_json({"b": 1, "a": [1.5, None]}) == '{"a":[1.5,null],"b":1}'
    with pytest.raises(ValueError):
        pv.canonical_json({"x": math.nan})


def test_yaml_config_needs_pyyaml(tmp_path):
    p = tmp_path / "s.yaml"
    p.write_text("name: x\n")
    try:
        import yaml  # noqa: F401
    except ImportError:
        with pytest.raises(RuntimeError, match="PyYAML"):
            pl.load_config(p)
    else:
        assert pl.load_config(p) == {"name": "x"}


# ------------------------------------------------------------- run_study
def test_tiny_outputs(tiny_run):
    out = tiny_run.out_dir
    for name in pl.RESULT_CSVS + ("manifest.json", "study.json", "run_record.json"):
        assert (out / name).exists(), name
    mol = pd.read_csv(out / "molecules.csv", dtype=str, keep_default_na=False)
    rej = mol[mol["status"] == "rejected"].set_index("id")["reason"].to_dict()
    assert rej == {"bad1": "unparsable", "bad2": "disallowed_element",
                   "dup1": "duplicate_inchikey:m05"}
    feats = pd.read_csv(out / "features.csv")
    assert len(feats) == 30 and feats["id"].tolist()[0] == "m01"
    # an undefined index is NaN with a reason, never 0
    fail = pd.read_csv(out / "failures.csv")
    r3 = fail[(fail["row"] == "m01") & (fail["feature"] == "R3_M1")]
    assert len(r3) == 1 and "k=3" in r3["reason"].iloc[0]
    assert math.isnan(feats.loc[feats["id"] == "m01", "ti.R3_M1"].iloc[0])


def test_tiny_summary_ladder(tiny_run):
    s = tiny_run.summary
    sets = ["null", "size", "TI", "TI_degree", "TI_reverse_degree",
            "TI_neighbourhood_degree_sum", "TI_distance", "TI_matching", "TI_entropy", "rdkit2d"]
    for target in ("MolLogP", "MolWt"):
        sub = s[s["target"] == target]
        assert sub["set"].tolist() == sets
        assert (sub["n"] == 30).all() and (sub["learner"] == "selected").all()
        ref = sub[sub["set"] == "size"].iloc[0]
        assert np.isnan(ref["d_r2_mean"])
        others = sub[sub["set"] != "size"]
        assert np.isfinite(others["d_r2_ci_low"]).all()
        assert (others["d_r2_ci_low"] <= others["d_r2_mean"]).all()
    cls = s.drop_duplicates("target").set_index("target")["target_class"].to_dict()
    assert cls == {"MolLogP": 2, "MolWt": 1}
    # rdkit2d exclusions follow the target
    study = json.loads((tiny_run.out_dir / "study.json").read_text())
    rd = study["feature_sets"]["rdkit2d"]
    assert "rdkit2d.MolLogP" not in rd["MolLogP"] and "rdkit2d.MolMR" not in rd["MolLogP"]
    assert "rdkit2d.MolWt" not in rd["MolWt"] and "rdkit2d.MolLogP" in rd["MolWt"]
    assert s.loc[(s["target"] == "MolLogP") & (s["set"] == "rdkit2d"), "p"].iloc[0] == len(rd["MolLogP"])


def test_folds_identical_across_sets(tiny_run):
    f = pd.read_csv(tiny_run.out_dir / "predictions.csv")
    f = f[f["target"] == "MolWt"]
    by_set = {k: g[["repeat", "fold", "row"]].reset_index(drop=True) for k, g in f.groupby("set")}
    first = next(iter(by_set.values()))
    assert all(v.equals(first) for v in by_set.values())


def test_refuses_non_empty_out_dir(tiny_run):
    with pytest.raises(FileExistsError):
        pl.run_study(_tiny_config(), tiny_run.out_dir, timestamp_utc=TS, base_dir=FIX)


# ------------------------------------------------------------- provenance
def test_run_record_complete(tiny_run):
    rec = json.loads((tiny_run.out_dir / "run_record.json").read_text())
    for key in ("schema", "created_utc", "topoqspr_version", "code", "environment", "chembl",
                "config", "config_sha256", "seeds", "inputs", "manifest_sha256", "outputs",
                "timings_s", "counts", "result_csvs", "base_dir"):
        assert key in rec, key
    assert rec["created_utc"] == TS
    code = rec["code"]
    assert set(code) >= {"commit", "dirty", "untracked"}
    if code["commit"] is not None:
        assert re.fullmatch(r"[0-9a-f]{40}", code["commit"]) and isinstance(code["dirty"], bool)
    else:
        assert code["reason"]
    env = rec["environment"]
    for pkg in ("numpy", "pandas", "scikit-learn", "scipy", "rdkit"):
        assert pkg in env["packages"], pkg
    assert env["python_version"] and env["platform"] and "PYTHONHASHSEED" in env["env"]
    assert rec["config_sha256"] == pv.config_sha256(rec["config"])
    assert rec["seeds"]["outer_seeds"] == [0]
    assert rec["chembl"] is None  # SMILES input
    assert rec["inputs"]["cohort"]["sha256"] == pv.sha256_file(TINY)
    assert rec["manifest_sha256"] == pv.sha256_file(tiny_run.out_dir / "manifest.csv")
    for name, info in rec["outputs"].items():
        assert pv.sha256_file(tiny_run.out_dir / name) == info["sha256"], name
    assert set(rec["result_csvs"]) == set(pl.RESULT_CSVS + pl.AD_CSVS)  # no protocol section
    assert {"cohort", "standardize", "features", "targets", "evaluate"} <= set(rec["timings_s"])
    assert pv.read_run_record(tiny_run.out_dir / "run_record.json")["config_sha256"] == rec["config_sha256"]


def test_tampered_run_record_rejected(tiny_run, tmp_path):
    rec = json.loads((tiny_run.out_dir / "run_record.json").read_text())
    rec["config"]["evaluation"]["seeds"] = [7]
    p = tmp_path / "run_record.json"
    p.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match="hash mismatch"):
        pv.read_run_record(p)


def test_rerun_detects_changed_input(tiny_run, tmp_path):
    (tmp_path / "tiny.csv").write_text(TINY.read_text().replace("m01,CCCC", "m01,CCCCC"))
    with pytest.raises(ValueError, match="differs"):
        pl.rerun(tiny_run.out_dir / "run_record.json", tmp_path / "out", timestamp_utc=TS,
                 base_dir=tmp_path)


# ---------------------------------------------------- pilot fixture + rerun
def test_pilot_cohort_from_recorded_cache(pilot_run):
    rec = pilot_run.run_record
    assert rec["chembl"]["release"] == "ChEMBL_37"
    assert rec["chembl"]["release_date"] == "2026-05-01"
    man = pd.read_csv(pilot_run.out_dir / "manifest.csv", dtype=str, keep_default_na=False)
    assert (man["decision"] == "include").sum() == rec["counts"]["candidates"] == 33
    side = json.loads((pilot_run.out_dir / "manifest.json").read_text())
    assert side["manifest_sha256"] == rec["manifest_sha256"]
    assert side["created_utc"] == TS


def test_rerun_byte_identical(pilot_run, tmp_path):
    rep = pl.rerun(pilot_run.out_dir / "run_record.json", tmp_path / "again",
                   timestamp_utc="2027-01-01T00:00:00Z", cache_dir=CHEMBL)
    assert rep["ok"], rep
    assert sorted(rep["identical"]) == sorted(pl.RESULT_CSVS + pl.AD_CSVS)
    for name in pl.RESULT_CSVS + pl.AD_CSVS:
        assert (tmp_path / "again" / name).read_bytes() == (pilot_run.out_dir / name).read_bytes()
    assert json.loads((tmp_path / "again" / "rerun_report.json").read_text())["ok"]


def test_rerun_needs_the_release_cache(pilot_run, tmp_path):
    from topoqspr.chembl import CacheMiss

    (tmp_path / "empty").mkdir()
    with pytest.raises(CacheMiss):
        pl.rerun(pilot_run.out_dir / "run_record.json", tmp_path / "x", timestamp_utc=TS,
                 cache_dir=tmp_path / "empty")


# ----------------------------------------------------------- target providers
class _Experimental:
    """A class-3 provider that has values for some molecules only."""

    name = "toy_experimental"

    def __init__(self, names=None, missing=("m02", "m03")):
        self.missing = set(missing)

    def specs(self):
        return [tp.TargetSpec("logS_toy", 3, None, "log10(mol/L)", "toy",
                              "toy experimental endpoint")]

    def compute(self, molecules, mols):
        vals, fails = [], []
        for rid, smi in zip(molecules.index, molecules["smiles"]):
            if rid in self.missing:
                vals.append(math.nan)
                fails.append({"row": rid, "feature": "logS_toy", "reason": "no measurement"})
            else:
                vals.append(-0.1 * len(smi) + 0.05 * smi.count("O"))
        return tp.TargetValues(pd.DataFrame({"logS_toy": vals}, index=molecules.index), fails)

    def provenance(self):
        return {"provider": self.name, "matching": "full InChIKey"}


def make_toy(names=None, **opts):
    return _Experimental(names, **opts)


def test_external_provider_plugs_in(tmp_path):
    cfg = _tiny_config(
        feature_sets=["null", "size"],
        targets=[{"provider": f"{__name__}:make_toy", "names": None, "options": {}}])
    res = pl.run_study(cfg, tmp_path / "res", timestamp_utc=TS, base_dir=FIX)
    s = res.summary
    assert (s["target"] == "logS_toy").all() and (s["target_class"] == 3).all()
    assert (s["n"] == 28).all() and (s["n_target_missing"] == 2).all()
    fail = pd.read_csv(tmp_path / "res" / "failures.csv")
    assert set(fail.loc[fail["feature"] == "logS_toy", "row"]) == {"m02", "m03"}
    tv = pd.read_csv(tmp_path / "res" / "targets.csv")
    assert tv.loc[tv["id"] == "m02", "logS_toy"].isna().all()  # NaN, not 0
    study = json.loads((tmp_path / "res" / "study.json").read_text())
    assert study["target_providers"] == [{"provider": "toy_experimental", "matching": "full InChIKey"}]


def test_target_spec_and_registry_checks():
    with pytest.raises(KeyError):
        tp.TargetSpec("x", 3, "NoSuchRule", "", "src")
    with pytest.raises(ValueError):
        tp.TargetSpec("x", 4, None, "", "src")
    with pytest.raises(KeyError):
        tp.get_target_provider("nope")
    with pytest.raises(KeyError):
        tp.RDKitTargetProvider(["MolWt", "XLogP3"])
    specs = {s.name: s for s in tp.RDKitTargetProvider().specs()}
    assert {k: v.target_class for k, v in specs.items()} == {
        "MolWt": 1, "MolLogP": 2, "TPSA": 2, "MolMR": 2}
    tp.register_target_provider("toy_registered", make_toy)
    try:
        assert isinstance(tp.get_target_provider("toy_registered"), tp.TargetProvider)
        with pytest.raises(ValueError):
            tp.register_target_provider("toy_registered", make_toy)
    finally:
        tp.PROVIDERS.pop("toy_registered")


def test_rdkit_targets_match_rdkit():
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors

    mols = [Chem.MolFromSmiles(s) for s in ("CCO", "c1ccccc1O")]
    table = pd.DataFrame({"inchikey": ["a", "b"]}, index=["x", "y"])
    v = tp.RDKitTargetProvider().compute(table, mols).values
    assert v.loc["x", "MolLogP"] == Crippen.MolLogP(mols[0])
    assert v.loc["y", "TPSA"] == Descriptors.TPSA(mols[1])


# ------------------------------------------------------------------- CLI
def test_cli_indices(tmp_path, capsys):
    out = tmp_path / "idx.csv"
    assert cli.main(["indices", "--smiles-csv", str(TINY), "--out", str(out)]) == 0
    df = pd.read_csv(out, dtype={"Z": object})
    assert len(df) == 33
    assert df.set_index("id").loc["bad1", "status"] == "rejected"
    assert df.set_index("id").loc["m01", "W"] == 10  # butane: Wiener index 10
    assert df.set_index("id").loc["m09", "Z"] == "18"  # benzene: Hosoya Z 18
    fail = pd.read_csv(tmp_path / "idx.failures.csv")
    assert "rejected: unparsable" in set(fail["reason"])


def test_cli_cohort_offline_deterministic(tmp_path, capsys):
    args = ["cohort", "--mesh", MIGRAINE, "--cache-dir", str(CHEMBL), "--offline",
            "--timestamp", TS]
    assert cli.main(args + ["--out", str(tmp_path / "a")]) == 0
    assert cli.main(args + ["--out", str(tmp_path / "b")]) == 0
    a = (tmp_path / "a" / "manifest.csv").read_bytes()
    assert a == (tmp_path / "b" / "manifest.csv").read_bytes()
    assert (tmp_path / "a" / "manifest.json").read_bytes() == (tmp_path / "b" / "manifest.json").read_bytes()
    assert "33 included" in capsys.readouterr().out


def test_cli_headings_offline(capsys):
    assert cli.main(["headings", "migraine", "--cache-dir", str(CHEMBL), "--offline"]) == 0
    out = capsys.readouterr().out
    assert f"{MIGRAINE},Migraine Disorders" in out


def test_cli_run_and_rerun(tmp_path, capsys):
    cfg = _tiny_config(feature_sets=["null", "size", "TI"],
                       targets=[{"provider": "rdkit", "names": ["TPSA"]}])
    cfg["cohort"]["path"] = str(TINY)
    cfg_path = tmp_path / "study.json"
    cfg_path.write_text(json.dumps(cfg))
    out = tmp_path / "res"
    assert cli.main(["run", "--config", str(cfg_path), "--out", str(out), "--timestamp", TS]) == 0
    assert (out / "run_record.json").exists()
    assert cli.main(["run", "--config", str(cfg_path), "--out", str(out)]) == 2  # not empty
    assert cli.main(["rerun", str(out / "run_record.json"), "--out", str(tmp_path / "re")]) == 0
    assert '"ok": true' in capsys.readouterr().out


def test_cli_bad_config_exit_code(tmp_path, capsys):
    p = tmp_path / "bad.json"
    p.write_text(json.dumps({"cohort": {"source": "chembl", "mesh_ids": ["D008881"]},
                             "evaluation": {"folds": 5}}))
    assert cli.main(["run", "--config", str(p), "--out", str(tmp_path / "o")]) == 2
    assert "unknown key" in capsys.readouterr().err


def test_cli_offline_cache_miss_exit_code(tmp_path, capsys):
    (tmp_path / "empty").mkdir()
    rc = cli.main(["cohort", "--mesh", MIGRAINE, "--cache-dir", str(tmp_path / "empty"),
                   "--offline", "--out", str(tmp_path / "o")])
    assert rc == 2 and "no cached ChEMBL status" in capsys.readouterr().err


def test_module_entry_point():
    env = {**os.environ, "PYTHONPATH": str(SRC)}
    r = subprocess.run([sys.executable, "-m", "topoqspr", "--version"], capture_output=True,
                       text=True, env=env, check=True)
    assert r.stdout.startswith("topoqspr ")


# ------------------------------------------------- regression tests (phase-1 review)
@pytest.mark.parametrize("over, match", [
    ({"evaluation": {"seeds": [0, 0, 0]}}, "duplicate outer seeds"),
    ({"evaluation": {"seeds": []}}, "seeds must be a non-empty"),
    ({"evaluation": {"outer_splits": 1}}, "outer_splits must be >= 2"),
    ({"evaluation": {"inner_splits": 1}}, "inner_splits must be >= 2"),
    ({"evaluation": {"families": [], "learners": ["selected"]}}, "families is empty"),
    ({"evaluation": {"learners": []}}, "learners is empty"),
    ({"targets": []}, "targets must be a non-empty list"),
    ({"targets": [{"provider": "rdkit", "names": []}]}, "'names' is empty"),
    ({"targets": [{"provider": "rdkit", "names": ["TPSA", "TPSA"]}]}, "duplicate names"),
])
def test_resolve_rejects_degenerate_evaluation(over, match):
    cfg = {"cohort": {"source": "chembl", "mesh_ids": [MIGRAINE]}, **over}
    with pytest.raises(ValueError, match=match):
        pl.resolve_config(cfg)


def test_non_mapping_config_is_a_value_error(tmp_path):
    with pytest.raises(ValueError, match="mapping"):
        pl.resolve_config([1])
    p = tmp_path / "list.json"
    p.write_text("[1]")
    with pytest.raises(ValueError, match="mapping"):
        pl.load_config(p)


def test_parquet_needs_an_engine_before_anything_runs(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(pl, "_parquet_engine", lambda: None)
    cfg = _tiny_config(output={"features_format": "parquet"})
    with pytest.raises(ValueError, match="pyarrow or fastparquet"):
        pl.resolve_config(cfg)
    cfg["cohort"]["path"] = str(TINY)
    p = tmp_path / "pq.json"
    p.write_text(json.dumps(cfg))
    out = tmp_path / "res"
    assert cli.main(["run", "--config", str(p), "--out", str(out)]) == 2
    assert "pyarrow or fastparquet" in capsys.readouterr().err
    assert not out.exists()


class _Failing(_Experimental):
    def compute(self, molecules, mols):
        raise RuntimeError("provider exploded")


def make_failing(names=None, **opts):
    return _Failing(names, **opts)


def test_failed_run_leaves_no_partial_results(tmp_path):
    cfg = _tiny_config(feature_sets=["null", "size"],
                       targets=[{"provider": f"{__name__}:make_failing"}])
    out = tmp_path / "res"
    with pytest.raises(RuntimeError, match="exploded"):
        pl.run_study(cfg, out, timestamp_utc=TS, base_dir=FIX)
    assert not out.exists()
    assert list(tmp_path.iterdir()) == []  # no staging directory left either
    # a failed overwrite keeps the previous results
    out.mkdir()
    (out / "keep.txt").write_text("old")
    with pytest.raises(RuntimeError):
        pl.run_study(cfg, out, timestamp_utc=TS, base_dir=FIX, overwrite=True)
    assert [x.name for x in out.iterdir()] == ["keep.txt"]


def test_cli_exit_codes_for_config_and_internal_errors(tmp_path, monkeypatch, capsys):
    p = tmp_path / "list.json"
    p.write_text("[]")
    assert cli.main(["run", "--config", str(p), "--out", str(tmp_path / "a")]) == 2
    cfg = _tiny_config(targets=[{"provider": "nosuch.mod:f"}])
    cfg["cohort"]["path"] = str(TINY)
    p = tmp_path / "mod.json"
    p.write_text(json.dumps(cfg))
    assert cli.main(["run", "--config", str(p), "--out", str(tmp_path / "b")]) == 2
    assert "could not be loaded" in capsys.readouterr().err
    assert not (tmp_path / "b").exists()  # provider checked before the cohort stage

    def boom(*a, **k):
        raise ZeroDivisionError("bug")

    monkeypatch.setattr(pl, "run_study", boom)
    p = tmp_path / "ok.json"
    p.write_text(json.dumps({"cohort": {"source": "smiles_csv", "path": str(TINY)}}))
    rc = cli.main(["run", "--config", str(p), "--out", str(tmp_path / "c")])
    assert rc == cli.EXIT_INTERNAL == 3
    assert "ZeroDivisionError" in capsys.readouterr().err


def test_read_results_keeps_the_null_set(tiny_run):
    assert pd.read_csv(tiny_run.out_dir / "summary.csv")["set"].isna().any()  # the pitfall
    s = pl.read_results(tiny_run.out_dir, "summary")
    assert "null" in set(s["set"]) and not s["set"].isna().any()
    assert s["d_r2_mean"].dtype == float  # empty cells are still NaN
    f = pl.read_results(tiny_run.out_dir, "folds.csv")
    assert (f["set"] == "null").sum() == 3 * 2  # 3 folds x 2 targets


class _NoReason(_Experimental):
    def compute(self, molecules, mols):
        tv = super().compute(molecules, mols)
        return tp.TargetValues(tv.values, [])  # NaN for m02, m03 but no reasons


def make_no_reason(names=None, **opts):
    return _NoReason(names, **opts)


def test_nan_target_without_reason_is_rejected(tmp_path):
    cfg = _tiny_config(feature_sets=["null", "size"],
                       targets=[{"provider": f"{__name__}:make_no_reason"}])
    with pytest.raises(ValueError, match="2 NaN value.*without a failure reason"):
        pl.run_study(cfg, tmp_path / "res", timestamp_utc=TS, base_dir=FIX)
    v = pd.DataFrame({"t": [1.0, math.inf]}, index=["a", "b"])
    with pytest.raises(ValueError, match="infinite"):
        pl._check_target_contract("x", v, [])


STEREO = {"citalopram": "CN(C)CCCC1(c2ccc(F)cc2)OCc2cc(C#N)ccc21",
          "escitalopram": "CN(C)CCC[C@@]1(c2ccc(F)cc2)OCc2cc(C#N)ccc21",
          "cetirizine": "OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1",
          "levocetirizine": "OC(=O)COCCN1CCN(CC1)[C@@H](c1ccccc1)c1ccc(Cl)cc1"}


def test_stereoisomers_share_a_fold(tmp_path):
    tiny = pd.read_csv(TINY, dtype=str, keep_default_na=False)
    extra = pd.DataFrame({"id": list(STEREO), "smiles": list(STEREO.values())})
    pd.concat([tiny, extra]).to_csv(tmp_path / "s.csv", index=False)
    cfg = _tiny_config(cohort={"source": "smiles_csv", "path": "s.csv"},
                       feature_sets=["null", "size", "composition", "TI"],
                       targets=[{"provider": "rdkit", "names": ["MolLogP"]}])
    cfg["evaluation"] = dict(FAST_EVAL, seeds=[0, 1, 2])
    res = pl.run_study(cfg, tmp_path / "res", timestamp_utc=TS, base_dir=tmp_path)
    mol = pl.read_results(res.out_dir, "molecules", dtype=str).set_index("id")
    grp = mol["connectivity_group"].fillna("")
    assert (mol.loc[list(STEREO), "status"] == "ok").all()  # all four kept
    assert grp["citalopram"] == grp["escitalopram"] == "citalopram"
    assert grp["cetirizine"] == grp["levocetirizine"] == "cetirizine"
    assert (grp.drop(list(STEREO)) == "").all()
    rec = res.run_record
    assert rec["counts"]["connectivity_groups"]["groups"] == 2
    assert rec["counts"]["connectivity_groups"]["molecules"] == 4
    assert rec["seeds"]["outer_splitter"].startswith("balanced_group_kfold (connectivity")
    pred = pl.read_results(res.out_dir, "predictions")
    fold = pred[pred["set"] == "size"].set_index(["repeat", "row"])["fold"]
    for r in range(3):
        assert fold[(r, "citalopram")] == fold[(r, "escitalopram")]
        assert fold[(r, "cetirizine")] == fold[(r, "levocetirizine")]
    # identical 2D features within a pair: exactly why they must not straddle folds
    feats = pl.read_results(res.out_dir, "features").set_index("id")
    cols = [c for c in feats.columns if c.startswith(("ti.", "size.", "composition."))]
    assert feats.loc["citalopram", cols].equals(feats.loc["escitalopram", cols])


def test_cv_groups_rules():
    from rdkit import Chem

    mols = [Chem.MolFromSmiles(s) for s in ("C[C@H](O)CC", "C[C@@H](O)CC", "c1ccccc1CC", "CCCO")]
    ids = ["a", "b", "c", "d"]
    conn = ["a", "a", "", ""]
    assert pl._cv_groups("random", mols, ["", "", "", ""], ids) is None  # plain KFold
    assert list(pl._cv_groups("random", mols, conn, ids)) == ["a", "a", "c", "d"]
    g = pl._cv_groups("scaffold", mols, conn, ids)
    assert g[0] == g[1] and len({g[0], g[2], g[3]}) == 3  # acyclic pair merged


def _pilot_manifest_copy(pilot_run, dst: Path) -> Path:
    dst.mkdir()
    for n in ("manifest.csv", "manifest.json"):
        (dst / n).write_bytes((pilot_run.out_dir / n).read_bytes())
    return dst / "manifest.csv"


def _manifest_config(path):
    return {"name": "from-manifest", "cohort": {"source": "manifest", "path": str(path)},
            "feature_sets": ["null", "size"], "targets": [{"provider": "rdkit", "names": ["TPSA"]}],
            "evaluation": dict(FAST_EVAL)}


def test_edited_manifest_is_refused(pilot_run, tmp_path):
    man = _pilot_manifest_copy(pilot_run, tmp_path / "m")
    res = pl.run_study(_manifest_config(man), tmp_path / "ok", timestamp_utc=TS)
    assert res.run_record["chembl"]["manifest_verified"] is True
    assert res.run_record["chembl"]["release"] == "ChEMBL_37"
    df = pd.read_csv(man, dtype=str, keep_default_na=False)
    i = df.index[df["decision"] == "include"][0]
    df.loc[i, "canonical_smiles"] = "CCCCCCCCCC"
    man.write_bytes(df.to_csv(index=False, lineterminator="\n").encode())
    with pytest.raises(ValueError, match="edited after it was built"):
        pl.run_study(_manifest_config(man), tmp_path / "bad", timestamp_utc=TS)
    assert not (tmp_path / "bad").exists()
    (tmp_path / "m" / "manifest.json").unlink()  # explicit opt-out: runs, labelled unverified
    res = pl.run_study(_manifest_config(man), tmp_path / "unverified", timestamp_utc=TS)
    assert res.run_record["chembl"]["manifest_verified"] is False


def test_chembl_licence_recorded(pilot_run):
    from topoqspr.chembl import CHEMBL_DATA_LICENCE

    assert "CC BY-SA 3.0" in CHEMBL_DATA_LICENCE
    assert pilot_run.run_record["chembl"]["data_licence"] == CHEMBL_DATA_LICENCE
    side = json.loads((pilot_run.out_dir / "manifest.json").read_text())
    assert side["data_licence"] == CHEMBL_DATA_LICENCE
    notice = (CHEMBL / "LICENSE").read_text()
    assert "CC BY-SA 3.0" in notice and "ChEMBL_37" in notice and "EMBL-EBI" in notice


def _dockerfile_build_info(dst: Path, commit: str, dirty: str) -> None:
    """Write _build_info.py exactly as the Dockerfile's RUN line does."""
    text = (REPO / "Dockerfile").read_text()
    m = re.search(r"printf ('[^']*')", text)
    assert m, "Dockerfile no longer writes _build_info.py with printf"
    subprocess.run(["sh", "-c", f'printf {m.group(1)} "$GIT_COMMIT" "$GIT_DIRTY" > "$OUT"'],
                   env={**os.environ, "GIT_COMMIT": commit, "GIT_DIRTY": dirty,
                        "OUT": str(dst / "_build_info.py")}, check=True)


def test_installed_copy_outside_git_names_the_build_commit(tmp_path):
    """The Docker situation: the package lives outside any git work tree."""
    import shutil

    site = tmp_path / "site"
    shutil.copytree(SRC / "topoqspr", site / "topoqspr",
                    ignore=shutil.ignore_patterns("__pycache__", "_build_info.py"))
    probe = "import json, topoqspr.provenance as p; print(json.dumps(p.git_info()))"

    def probe_git_info():
        r = subprocess.run([sys.executable, "-c", probe], cwd=tmp_path, capture_output=True,
                           text=True, env={**os.environ, "PYTHONPATH": str(site)}, check=True)
        return json.loads(r.stdout)

    no_info = probe_git_info()
    assert no_info["commit"] is None and "_build_info" in no_info["reason"]
    commit = "0123456789abcdef0123456789abcdef01234567"
    _dockerfile_build_info(site / "topoqspr", commit, "1")
    info = probe_git_info()
    assert info["commit"] == commit and info["dirty"] is True
    assert info["source"] == "build_info"
    # inside the repository the live git state is used, never a build-info file
    assert pv.git_info()["source"] in ("git", None)


# ------------------------------------------------- expdata provider (class 3)
def _std_key(smi):
    from topoqspr.graph import standardize

    return standardize(smi).inchikey


def _fake_expdata_loader(name, cache_dir, *, offline, drop_disagreeing, keep_salt_forms,
                         client=None):
    """get_target stand-in: values for the tiny molecules m01-m24 (as get_target would
    return them: one row per full InChIKey, sorted), escitalopram only (not citalopram),
    and one parent (m02) whose replicates disagree."""
    tiny = pd.read_csv(TINY, dtype=str, keep_default_na=False).set_index("id")["smiles"]
    rows = []
    for i in range(1, 25):
        rows.append(dict(inchikey=_std_key(tiny[f"m{i:02d}"]),
                         value=float(i) + (0.5 if name == "mp_bradley" else 0.0),
                         n_measurements=1, source="fake", sd=math.nan,
                         replicates_disagree=False))
    rows[1].update(sd=2.0, replicates_disagree=True, n_measurements=2)  # m02
    rows.append(dict(inchikey=_std_key(STEREO["escitalopram"]), value=-1.0, n_measurements=1,
                     source="fake", sd=math.nan, replicates_disagree=False))
    df = pd.DataFrame(rows).sort_values("inchikey").reset_index(drop=True)
    df.attrs.update(raw_sha256="f" * 64, processing_fingerprint="fp-test")
    return df


def make_fake_expdata(names=None, **opts):
    return tp.ExpdataTargetProvider(names, loader=_fake_expdata_loader, **opts)


def _stereo_csv(tmp_path):
    tiny = pd.read_csv(TINY, dtype=str, keep_default_na=False)
    extra = pd.DataFrame({"id": list(STEREO), "smiles": list(STEREO.values())})
    pd.concat([tiny, extra]).to_csv(tmp_path / "s.csv", index=False)


def test_expdata_provider_full_key_join_and_reasons(tmp_path):
    _stereo_csv(tmp_path)
    cfg = _tiny_config(cohort={"source": "smiles_csv", "path": "s.csv"},
                       feature_sets=["null", "size"],
                       targets=[{"provider": f"{__name__}:make_fake_expdata",
                                 "names": ["logs_aqsoldb", "mp_bradley"]}])
    res = pl.run_study(cfg, tmp_path / "res", timestamp_utc=TS, base_dir=tmp_path)
    tv = pl.read_results(res.out_dir, "targets").set_index("id")
    assert tv.loc["m01", "logs_aqsoldb"] == 1.0 and tv.loc["m01", "mp_bradley"] == 1.5
    assert tv.loc["escitalopram", "logs_aqsoldb"] == -1.0
    for rid in ("m02", "citalopram", "m25", "cetirizine"):  # NaN, never 0
        assert math.isnan(tv.loc[rid, "logs_aqsoldb"])
    fail = pl.read_results(res.out_dir, "failures")
    why = fail[fail["feature"] == "logs_aqsoldb"].set_index("row")["reason"]
    assert why["m02"].startswith("replicates_disagree")
    assert why["citalopram"].startswith("no_measurement_stereo_mismatch")
    assert why["m25"].startswith("no_measurement:")
    assert "escitalopram" not in why.index and "m01" not in why.index
    s = res.summary
    assert (s["target_class"] == 3).all() and set(s["target"]) == {"logs_aqsoldb", "mp_bradley"}
    assert (s["n"] == 24).all()  # m01, m03-m24, escitalopram
    rec = res.run_record
    tg = {t["name"]: t for t in rec["targets"]}
    assert tg["logs_aqsoldb"]["target_class"] == 3
    assert tg["logs_aqsoldb"]["source"] == "expdata:aqsoldb"
    assert tg["logs_aqsoldb"]["rdkit2d_exclusion"] is None
    prov = rec["target_providers"][0]
    assert prov["provider"] == "expdata" and "full standard InChIKey" in prov["matching"]
    by_t = {t["target"]: t for t in prov["targets"]}
    t0 = by_t["logs_aqsoldb"]
    assert t0["licence"] == "CC0 1.0" and t0["doi"] == "10.1038/s41597-019-0151-1"
    assert t0["raw_sha256"] == "f" * 64 and t0["n_values"] == 24
    assert t0["n_dropped_disagreeing"] == 1 and t0["n_matched_full_key"] == 25
    assert by_t["mp_bradley"]["licence"] == "CC0"


def test_expdata_provider_options_and_client():
    p = tp.get_target_provider("expdata", ["logd74_az"], {"drop_disagreeing": False})
    assert [s.name for s in p.specs()] == ["logd74_az"] and p.specs()[0].target_class == 3
    with pytest.raises(KeyError, match="unknown expdata targets"):
        tp.get_target_provider("expdata", ["logS_esol"])
    seen = {}

    def loader(name, cache_dir, **kw):
        seen.update(kw, cache_dir=cache_dir)
        return pd.DataFrame(dict(inchikey=["AAAAAAAAAAAAAA-BBBBBBBBBB-N"], value=[9.0], sd=[3.0],
                                 replicates_disagree=[True]))

    class _Client:
        release = "ChEMBL_37"

    prov = tp.ExpdataTargetProvider(["logd74_az"], drop_disagreeing=False, loader=loader)
    prov.set_runtime(cache_dir="/c", offline=True, client=_Client())
    mt = pd.DataFrame({"inchikey": ["AAAAAAAAAAAAAA-BBBBBBBBBB-N", ""]},
                      index=pd.Index(["a", "b"], name="id"))
    out = prov.compute(mt, [None, None])
    assert out.values["logd74_az"].tolist()[0] == 9.0  # kept: drop_disagreeing=False
    assert out.failures == [{"row": "b", "feature": "logd74_az",
                             "reason": "no_inchikey: standardised parent has no InChIKey"}]
    assert seen["client"].release == "ChEMBL_37" and seen["offline"] is True
    assert seen["cache_dir"] == "/c" and seen["drop_disagreeing"] is False
    pinned = tp.ExpdataTargetProvider(["logd74_az"], chembl_release="ChEMBL_36", loader=loader)
    pinned.set_runtime(client=_Client())
    with pytest.raises(RuntimeError, match="pins ChEMBL_36"):
        pinned.compute(mt, [None, None])
    json.dumps(prov.provenance(), allow_nan=False)  # run-record safe


# ---------------------------------------------------------- pooled cohort
class _PooledClient:
    """ChEMBL stand-in for expdata.approved_drug_pool: the first 24 tiny molecules as
    approved parents (one with a salt form), plus a protein."""

    release = "ChEMBL_37"
    offline = True

    def __init__(self):
        from rdkit import Chem

        tiny = pd.read_csv(TINY, dtype=str, keep_default_na=False)
        self.db, rows = {}, []
        for i, smi in enumerate(tiny["smiles"][:24], start=1):
            cid = f"CHEMBL{i}"
            self.db[cid] = self._rec(cid, smi, cid, Chem.MolToInchiKey(Chem.MolFromSmiles(smi)))
            rows.append(dict(drugind_id=i, molecule_chembl_id=cid, parent_molecule_chembl_id=cid,
                             max_phase_for_ind="4.0", mesh_id=f"D{i % 3}"))
        self.db["CHEMBL900"] = self._rec("CHEMBL900", "CCN.Cl", "CHEMBL22", "SALT-KEY")
        rows.append(dict(drugind_id=90, molecule_chembl_id="CHEMBL900",
                         parent_molecule_chembl_id="CHEMBL22", max_phase_for_ind="4.0",
                         mesh_id="D9"))
        self.db["CHEMBL999"] = self._rec("CHEMBL999", "", "CHEMBL999", "", mtype="Protein")
        rows.append(dict(drugind_id=99, molecule_chembl_id="CHEMBL999",
                         parent_molecule_chembl_id="CHEMBL999", max_phase_for_ind="4.0",
                         mesh_id="D1"))
        self.rows = rows

    @staticmethod
    def _rec(cid, smi, parent, key, mtype="Small molecule"):
        return dict(molecule_chembl_id=cid, pref_name=cid, molecule_type=mtype,
                    structure_type="MOL", inorganic_flag=0, max_phase="4.0",
                    atc_classifications=[], molecule_hierarchy={"parent_chembl_id": parent},
                    molecule_structures={"canonical_smiles": smi, "standard_inchi_key": key})

    def fetch_all(self, endpoint, params):
        assert endpoint == "drug_indication" and params["max_phase_for_ind"] == "4"
        return list(self.rows)

    def molecules(self, ids):
        return {i: self.db.get(i) for i in ids}

    def verify_release(self):
        pass

    def status(self):
        return {"chembl_db_version": self.release, "chembl_release_date": "2026-05-01"}


def test_pooled_cohort_source(tmp_path):
    cfg = _tiny_config(cohort={"source": "pooled", "chembl_release": "ChEMBL_37"},
                       feature_sets=["null", "size"],
                       targets=[{"provider": "rdkit", "names": ["MolLogP"]}])
    res = pl.run_study(cfg, tmp_path / "res", timestamp_utc=TS, client=_PooledClient())
    side = json.loads((res.out_dir / "manifest.json").read_text())
    assert side["kind"] == "pooled" and side["chembl_release"] == "ChEMBL_37"
    assert side["indication_rows"] == 26 and side["counts"]["included"] == 24
    man = pl.read_results(res.out_dir, "manifest", dtype=str).set_index("molecule_chembl_id")
    assert man.loc["CHEMBL900", "reason"] == "collapsed_to_parent"
    assert man.loc["CHEMBL999", "reason"] == "not_small_molecule"
    assert res.run_record["chembl"]["release"] == "ChEMBL_37"
    assert res.run_record["config"]["cohort"] == {"source": "pooled", "chembl_release": "ChEMBL_37"}
    assert (res.summary["n"] == 24).all()
    cfg["cohort"]["chembl_release"] = "ChEMBL_36"
    with pytest.raises(RuntimeError, match="pins ChEMBL_36"):
        pl.run_study(cfg, tmp_path / "res2", timestamp_utc=TS, client=_PooledClient())


def test_target_with_too_few_values_is_skipped_not_crashed(tmp_path):
    cfg = _tiny_config(feature_sets=["null", "size"],
                       targets=[{"provider": f"{__name__}:make_toy", "names": None,
                                 "options": {"missing": [f"m{i:02d}" for i in range(1, 29)]}}],
                       protocol={"family": ["logS_toy:size-null"],
                                 "permutation": {"n_permutations": 2},
                                 "bootstrap": {"n_boot": 2}})
    cfg["evaluation"] = dict(FAST_EVAL, reference="null")
    res = pl.run_study(cfg, tmp_path / "res", timestamp_utc=TS, base_dir=FIX)
    assert "logS_toy" in res.run_record["counts"]["skipped_targets"]
    assert len(res.summary) == 0
    conf = pl.read_results(res.out_dir, "confirmatory")
    assert len(conf) == 1 and "target not evaluated" in conf["reason"].iloc[0]
    assert math.isnan(conf["d_r2_mean"].iloc[0])
