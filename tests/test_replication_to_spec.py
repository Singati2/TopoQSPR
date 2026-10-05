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
