"""The replication folders must convert into specs the Experiment D harness accepts (the two
packages were built separately; this is the integration test that was missing)."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from replication_to_spec import convert, map_label, normalise_label  # noqa: E402

from topoqspr.experiments import _check_spec  # noqa: E402


@pytest.mark.parametrize("label,expected", [
    ("ND1(G)", "ND_1"), ("ND_4", "ND_4"), ("ABC^d", "ABC"), ("ReZ_{1}", "ReZ1"), ("NF*(G)", "NF"),
    ("mM2^d", "mM2"), ("SO^{nd}", "SO"),
])
def test_normalise_label(label, expected):
    assert normalise_label(label) == expected


def test_map_label_handles_neighbourhood_suffix():
    assert map_label("ABC^nd")[0] == "S_ABC" and map_label("ABC^d")[0] == "ABC"
    # no silent fall-back from X^nd to the degree index X (Godlin 2026 M1^nd was mapped to M1)
    assert map_label("M1^nd")[0] is None and map_label("M1^d")[0] == "M1"
    assert map_label("ND5(G)")[0] == "S_SDD" and map_label("ECI")[0] is None


@pytest.mark.parametrize("study", ["kabeer2026_antidepressants", "almulla2026_antifungal", "godlin2026_mdd"])
def test_in_scope_studies_convert_and_validate(study):
    spec, report = convert(ROOT / "data" / "replication" / study)
    assert not report["blocking_unmapped_predictors"]
    assert report["n_compounds"] >= 10 and report["n_comparable_fits"] > 0
    _check_spec(spec)  # the harness's own validator


def test_out_of_scope_study_is_refused_not_guessed(tmp_path):
    spec, report = convert(ROOT / "data" / "replication" / "vimala2026_anticancer")
    assert report["blocking_unmapped_predictors"]
    r = subprocess.run([sys.executable, str(ROOT / "tools" / "replication_to_spec.py"),
                        str(ROOT / "data" / "replication" / "vimala2026_anticancer"), str(tmp_path / "v.json")],
                       capture_output=True, text=True)
    assert r.returncode == 2 and "refused" in r.stderr and not (tmp_path / "v.json").exists()


def test_unsupported_reported_fit_form_is_listed_not_refused():
    """A study whose paper reports a logarithmic/ridge fit must still load; the harness lists
    such fits as not compared (regression for the Kabeer 2026 run that was refused)."""
    spec, _ = convert(ROOT / "data" / "replication" / "kabeer2026_antidepressants")
    assert any(f["form"] == "logarithmic" for f in spec["reported_fits"])
    checked = _check_spec(spec)
    assert any(f["form"] == "logarithmic" for f in checked["reported_fits"])
    assert "logarithmic" not in checked["model_forms"]


def test_statistics_keys_are_case_insensitive():
    """Papers print R2, r2 or R^2; Godlin 2026 was listed as 'no reported r2' because of case."""
    from replication_to_spec import _stat
    assert _stat({"r2": "0.963"}, "r2", "r^2") == 0.963
    assert _stat({"R2": "0.773"}, "r2") == 0.773
    assert _stat({"R": "x"}, "r") is None and _stat({}, "r2") is None


def test_loocv_models_go_to_the_report_not_the_in_sample_fits():
    spec, report = convert(ROOT / "data" / "replication" / "godlin2026_mdd")
    assert all(f.get("r2") is not None for f in spec["reported_fits"]) and len(spec["reported_fits"]) == 9
    assert len(report["reported_loocv"]) == 9 and all(r["q2_loo"] is not None for r in report["reported_loocv"])
    assert len(report["reported_other_validation"]) == 9  # external-set rows: listed only
    # Kabeer's Table 14 (LOO) rows are no longer duplicates of the in-sample fits
    spec, report = convert(ROOT / "data" / "replication" / "kabeer2026_antidepressants")
    assert len(spec["reported_fits"]) == 12 and len(report["reported_loocv"]) == 11


@pytest.mark.parametrize("study,n_values", [("almulla2026_antifungal", 180), ("kabeer2026_antidepressants", 225),
                                            ("godlin2026_mdd", 21 * 24)])
def test_printed_index_values_are_passed_on_by_registry_name(study, n_values):
    spec, report = convert(ROOT / "data" / "replication" / study)
    riv = spec["reported_index_values"]
    assert report["n_reported_index_values"] == n_values and len(riv) == report["n_compounds"]
    names = {c["name"] for c in spec["compounds"]}
    assert set(riv) <= names
    mapped = set(report["index_map"].values())
    keys = {k for v in riv.values() for k in v}
    assert mapped <= keys                      # every registry index the study uses has printed values
    assert all(isinstance(x, float) for v in riv.values() for x in v.values())
    _check_spec(spec)
