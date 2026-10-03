"""Tests for the Experiment D replication data (data/replication/, tools/extract_replication.py).

Offline: schema and consistency validation of all ten study folders against the committed
PubChem cache, deterministic re-resolution, and unit tests of the JATS rendering, value parsing,
property extraction and resolution rules on synthetic inputs
(``tests/fixtures/replication/pubchem_cache.json`` is a hand-made cache, not PubChem data).
The live re-extraction from Europe PMC is marked ``network``.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import re
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE_CACHE = ROOT / "tests" / "fixtures" / "replication" / "pubchem_cache.json"


def _load_tool():
    spec = importlib.util.spec_from_file_location("extract_replication",
                                                  ROOT / "tools" / "extract_replication.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["extract_replication"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


er = _load_tool()
STUDY_DIRS = er.study_dirs()


def _shortlist_rows() -> dict[int, str]:
    """Row number -> DOI from docs/REPLICATION_SHORTLIST.md."""
    out = {}
    for line in (ROOT / "docs" / "REPLICATION_SHORTLIST.md").read_text().splitlines():
        m = re.match(r"\|\s*(\d+)\s*\|[^|]*\|\s*\[([^\]]+)\]", line)
        if m:
            out[int(m.group(1))] = m.group(2)
    return out


# --------------------------------------------------------------------------- the ten folders

def test_ten_studies_rows_1_to_10():
    assert len(STUDY_DIRS) == 10
    rows = _shortlist_rows()
    seen = {}
    for d in STUDY_DIRS:
        s = er.load_study(d)
        seen[s["shortlist_row"]] = s
        assert s["citation"]["doi"] == rows[s["shortlist_row"]]
        assert s["status"] == ("proposed" if s["shortlist_row"] <= 5 else "reserve")
    assert sorted(seen) == list(range(1, 11))


@pytest.mark.parametrize("d", STUDY_DIRS, ids=[d.name for d in STUDY_DIRS])
def test_study_folder_validates_offline(d):
    cache = er.PubChemCache(offline=True)
    assert er.validate_study(d, cache) == []


@pytest.mark.parametrize("d", STUDY_DIRS, ids=[d.name for d in STUDY_DIRS])
def test_resolution_is_deterministic_from_cache(d):
    cache = er.PubChemCache(offline=True)
    for row in er.read_csv(d / "drugs.csv"):
        out = er.resolve_drug(row, cache)
        assert all(out[c] == row[c] for c in er.DRUG_OUTPUT_COLUMNS), row["drug_id"]


@pytest.mark.parametrize("d", STUDY_DIRS, ids=[d.name for d in STUDY_DIRS])
def test_licence_and_files(d):
    s = er.load_study(d)
    assert (d / "SOURCE.md").exists()
    src = (d / "SOURCE.md").read_text()
    assert s["citation"]["doi"] in src and s["licence"]["statement_verbatim"][:40] in src
    if s["licence"]["spdx"] == "CC-BY-4.0":
        assert s["licence"]["data_redistribution"] == "stored"
        assert (d / "properties.csv").exists()
    else:
        assert s["licence"]["data_redistribution"] == "withheld"
        assert not (d / "properties.csv").exists()
        assert s["reported_tables"] == []


@pytest.mark.parametrize("d", STUDY_DIRS, ids=[d.name for d in STUDY_DIRS])
def test_missing_values_are_nan_with_reason_never_zero(d):
    p = d / "properties.csv"
    if not p.exists():
        pytest.skip("licence-withheld study")
    for r in er.read_csv(p):
        if r["value"] == "":
            assert r["missing_reason"] in {"empty_cell", "dash_in_table", "unparseable"}
        else:
            assert r["missing_reason"] == ""
            v, _ = er.parse_value(r["value_verbatim"])
            assert float(r["value"]) == v


def test_known_resolution_cases():
    by = {d.name: {r["drug_id"]: r for r in er.read_csv(d / "drugs.csv")} for d in STUDY_DIRS}
    risp = by["colakoglu2026_bipolar"]["16"]
    assert (risp["name_verbatim"], risp["rule"], risp["pubchem_title"]) == \
        ("riseridone", "typo_corrected", "Risperidone")
    valp = by["ali2026_stiff_person"]["09"]
    assert valp["rule"] == "parent_form" and valp["molecular_formula"] == "C8H16O2"
    assert by["tang2026_arthritis"]["34"]["status"] == "unresolved"       # chondroitin sulfate
    assert by["supriya2026_tuberculosis"]["06"]["status"] == "unresolved"  # imipenem-cilastatin
    acs = er.load_study(er.DATA_DIR / "ahmed2026_antiepileptic")["automated_checks"]
    assert acs["duplicate_structures"] == [["10", "18"]]                   # lacosamide twice


def test_shared_connectivity_flagged():
    s = er.load_study(er.DATA_DIR / "godlin2026_mdd")
    assert ["07", "18"] in s["automated_checks"]["shared_connectivity_groups"]  # (es)citalopram


# --------------------------------------------------------------------------- unit tests

def test_parse_value():
    assert er.parse_value("−1.10") == (-1.1, "")
    assert er.parse_value("− 3.259") == (-3.259, "")
    assert er.parse_value("–0.310101") == (-0.310101, "")
    v, why = er.parse_value("–")
    assert math.isnan(v) and why == "dash_in_table"
    v, why = er.parse_value("")
    assert math.isnan(v) and why == "empty_cell"
    v, why = er.parse_value("12a")
    assert math.isnan(v) and why == "unparseable"


def test_matches_printed():
    assert er.matches_printed("306.23", "306.2")    # PubChem prints 306.2
    assert er.matches_printed("437", "436.9")
    assert not er.matches_printed("159.12", "159.23")
    assert not er.matches_printed("296", "308")


JATS = b"""<article xmlns:mml="http://www.w3.org/1998/Math/MathML"
 xmlns:xlink="http://www.w3.org/1999/xlink"><body>
<table-wrap id="t1"><label>2</label><caption><p>Props</p></caption><table>
<thead><tr><th>Drug</th><th>BP</th><th><inline-formula><mml:math><mml:msub><mml:mi>M</mml:mi><mml:mn>1</mml:mn></mml:msub></mml:math></inline-formula></th></tr></thead>
<tbody>
<tr><td>A</td><td>101.5</td><td>&#8722;2</td></tr>
<tr><td>B</td><td>&#8211;</td><td><inline-formula><alternatives><tex-math>\\documentclass{minimal}\\begin{document}$$x^2$$\\end{document}</tex-math><mml:math><mml:mi>x</mml:mi></mml:math></alternatives></inline-formula></td></tr>
<tr><td>C</td><td></td><td><inline-graphic xlink:href="img.gif"/></td></tr>
</tbody></table></table-wrap>
<p><mml:math><mml:mfrac><mml:mn>1</mml:mn><mml:msqrt><mml:mi>d</mml:mi></mml:msqrt></mml:mfrac></mml:math></p>
</body></article>"""


def test_read_table_and_render():
    t = er.read_table(JATS, "t1")
    assert t["label"] == "2" and t["caption"] == "Props"
    assert t["rows"][0] == ["Drug", "BP", "$M_{1}$"]
    assert t["rows"][2][2] == "$x^2$"           # TeX preferred over MathML
    assert t["rows"][3][2] == "[IMG img.gif]"
    import xml.etree.ElementTree as ET
    p = ET.fromstring(JATS).find("body/p")
    assert er.render(p) == "$(1)/(sqrt(d))$"
    with pytest.raises(KeyError):
        er.read_table(JATS, "nope")


def _drugs(labels):
    out = []
    for i, lab in enumerate(labels, 1):
        r = {c: "" for c in er.DRUG_COLUMNS}
        r.update(drug_id=f"{i:02d}", row_label_verbatim=lab)
        out.append(r)
    return out


def test_extract_properties():
    spec = {"xml_table_id": "t1", "skip_rows": 1, "label_column": 0,
            "columns": [{"property": "BP", "column": 1}]}
    rows = er.extract_properties(JATS, spec, _drugs(["A", "B", "C"]))
    assert [(r["value"], r["missing_reason"]) for r in rows] == \
        [("101.5", ""), ("", "dash_in_table"), ("", "empty_cell")]
    with pytest.raises(ValueError):
        er.extract_properties(JATS, spec, _drugs(["A", "X", "C"]))
    with pytest.raises(ValueError):
        er.extract_properties(JATS, spec, _drugs(["A", "B"]))


def _row(name, rule="as_printed", query=None):
    r = {c: "" for c in er.DRUG_COLUMNS}
    r.update(drug_id="01", name_verbatim=name, query_name=query if query is not None else name,
             rule=rule)
    return r


def test_resolution_rules_with_fixture_cache():
    cache = er.PubChemCache(FIXTURE_CACHE, offline=True)
    out = er.resolve_drug(_row("ethanol"), cache)
    assert (out["status"], out["pubchem_cid"], out["n_components"]) == ("resolved", "1", "1")
    out = er.resolve_drug(_row("spirits"), cache)
    assert (out["status"], out["pubchem_cid"], out["n_cids"]) == ("ambiguous_first_cid", "1", "2")
    out = er.resolve_drug(_row("sodium acetate"), cache)
    assert out["status"] == "resolved_multicomponent" and out["n_components"] == "2"
    out = er.resolve_drug(_row("ethanolum"), cache)
    assert out["status"] == "unresolved" and out["n_cids"] == "0" and out["pubchem_cid"] == ""
    out = er.resolve_drug(_row("ethanolum", "typo_corrected", "ethanol"), cache)
    assert out["pubchem_cid"] == "1"
    out = er.resolve_drug(_row("gin and tonic", "unresolvable", ""), cache)
    assert out["status"] == "unresolved"
    with pytest.raises(er.CacheMissError):
        er.resolve_drug(_row("methanol"), cache)


def test_cache_fills_and_saves(tmp_path):
    calls = []

    def getter(url):
        calls.append(url)
        if "/name/" in url:
            return 200, {"IdentifierList": {"CID": [1]}}
        return 200, json.loads(FIXTURE_CACHE.read_text())[er.props_path(1)]["json"]

    p = tmp_path / "cache.json"
    cache = er.PubChemCache(p, offline=False, getter=getter)
    out = er.resolve_drug(_row("ethyl alcohol"), cache)
    assert out["pubchem_title"] == "Ethanol" and len(calls) == 2
    cache.save()
    again = er.PubChemCache(p, offline=True)
    assert er.resolve_drug(_row("ethyl alcohol"), again) == out
    assert er.name_path("ethyl alcohol") == "/compound/name/ethyl%20alcohol/cids/JSON"


def test_automated_checks_synthetic():
    cache = er.PubChemCache(FIXTURE_CACHE, offline=True)
    drugs = []
    for i, name in enumerate(["2-butanol", "R-2-butanol", "ethanol"], 1):
        r = _row(name)
        r.update(drug_id=f"{i:02d}", row_label_verbatim=name)
        r.update(er.resolve_drug(r, cache))
        drugs.append(r)
    drugs[2]["formula_verbatim"] = "$C_{2}H_{6}O$"
    drugs[0]["formula_verbatim"] = "C4H8O"
    study = {"n_reported": 3, "properties": [
        {"code": "MW", "pubchem_field": "MolecularWeight"}, {"code": "BP"}]}
    props = []
    for d, mw, bp in zip(drugs, ["74.12", "74.07", "46.07"], ["99", "99", "78"]):
        for code, v in (("MW", mw), ("BP", bp)):
            val, why = er.parse_value(v)
            props.append({"drug_id": d["drug_id"], "property": code, "value_verbatim": v,
                          "value": er.fmt_float(val), "missing_reason": why})
    a = er.automated_checks(study, drugs, props)
    assert a["shared_connectivity_groups"] == [["01", "02"]]
    assert a["duplicate_structures"] == []
    assert a["formula_checked"] == 2 and a["formula_mismatch"] == ["01"]
    mw = a["properties"]["MW"]
    assert mw["pubchem_match"] == 2 and mw["pubchem_mismatch"] == ["02"]
    assert mw["mismatch_equal_to_exact_mass"] == ["02"]    # 74.07 = monoisotopic mass
    assert a["identical_property_rows"] == []
    assert er.automated_checks(study, drugs, None)["properties"].startswith("withheld")


def test_validator_catches_corruption(tmp_path):
    src = er.DATA_DIR / "ali2026_stiff_person"
    d = tmp_path / src.name
    shutil.copytree(src, d)
    cache = er.PubChemCache(offline=True)
    assert er.validate_study(d, cache) == []
    # a value that no longer matches its printed string
    rows = er.read_csv(d / "properties.csv")
    rows[0]["value"] = "0.0"
    er.write_csv(d / "properties.csv", rows, er.PROPERTY_COLUMNS)
    errs = er.validate_study(d, cache)
    assert any("inconsistent with value_verbatim" in e for e in errs)
    # an unknown resolution status and a tampered CID
    shutil.copy(src / "properties.csv", d / "properties.csv")
    drugs = er.read_csv(d / "drugs.csv")
    drugs[0]["pubchem_cid"] = "999"
    drugs[1]["status"] = "guessed"
    er.write_csv(d / "drugs.csv", drugs, er.DRUG_COLUMNS)
    errs = er.validate_study(d, cache)
    assert any("status 'guessed'" in e for e in errs)
    assert any("pubchem_cid '999' != recomputed" in e for e in errs)
    # a withheld licence must not ship property values
    shutil.copy(src / "drugs.csv", d / "drugs.csv")
    s = er.load_study(d)
    s["licence"]["data_redistribution"] = "withheld"
    er.dump_study(d, s)
    assert any("committed although licence is withheld" in e for e in er.validate_study(d, cache))
    # a model with a non-verbatim statistic
    s["licence"]["data_redistribution"] = "stored"
    s["models"][0]["statistics"]["R2"] = 0.98
    er.dump_study(d, s)
    assert any("verbatim strings" in e for e in er.validate_study(d, cache))


def test_drugs_csv_header_is_exact():
    for d in STUDY_DIRS:
        with open(d / "drugs.csv", newline="", encoding="utf-8") as fh:
            assert next(csv.reader(fh)) == er.DRUG_COLUMNS


@pytest.mark.network
@pytest.mark.parametrize("key", ["colakoglu2026_bipolar", "almulla2026_antifungal"])
def test_live_xml_reextraction_matches(key):  # pragma: no cover - network
    import hashlib

    import requests
    d = er.DATA_DIR / key
    s = er.load_study(d)
    r = requests.get(er.EPMC.format(pmcid=s["citation"]["pmcid"]), timeout=600)
    r.raise_for_status()
    rows = er.extract_properties(r.content, s["property_table"], er.read_csv(d / "drugs.csv"))
    assert [{c: x[c] for c in er.PROPERTY_COLUMNS} for x in rows] == er.read_csv(d / "properties.csv")
    if hashlib.sha256(r.content).hexdigest() != s["source"]["xml_sha256"]:
        pytest.xfail("Europe PMC re-rendered the XML; table content still matches")
