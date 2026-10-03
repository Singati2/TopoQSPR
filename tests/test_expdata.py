"""Tests for topoqspr.expdata (measured property targets).

Offline by default.  ``tests/fixtures/expdata/`` holds small verbatim excerpts of the two
CC0 sources (AqSolDB, Bradley); OPERA, Biogen and ChEMBL records are synthetic.  The live
download test is marked ``network``.
"""

from __future__ import annotations

import dataclasses
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from topoqspr import expdata as ed
from topoqspr.chembl import CacheMiss
from topoqspr.expdata import (
    RECORD_COLUMNS,
    SOURCES,
    STATUS_CODES,
    TARGET_COLUMNS,
    TARGETS,
    ChecksumError,
    aggregate,
    convert_unit,
    coverage,
    fetch,
    get_target,
    parse_aqsoldb,
    parse_bradley,
    parse_biogen,
    parse_chembl_activities,
    parse_opera,
    process_records,
    read_sdf_properties,
    sha256_bytes,
)
from topoqspr.graph import standardize

FIX = Path(__file__).resolve().parent / "fixtures" / "expdata"
LOGS = TARGETS["logs_aqsoldb"]
LOGP = TARGETS["logp_opera"]
MP = TARGETS["mp_bradley"]


# ----------------------------------------------------------------------------- units
@pytest.mark.parametrize(
    "value, unit_from, unit_to, mw, expected",
    [
        (1000.0, "nM", "log10(mol/L)", None, -6.0),
        (81300.0, "nM", "log10(mol/L)", None, math.log10(81.3e-6)),
        (10.0, "uM", "log10(mol/L)", None, -5.0),
        (1.0, "mol/L", "log10(mol/L)", None, 0.0),
        (-2.0, "log10(mol/L)", "log10(mol/L)", None, -2.0),
        (100.0, "ug/mL", "log10(mol/L)", 100.0, -3.0),  # 0.1 g/L / 100 g/mol = 1e-3 M
        (2.0, "log10(ug/mL)", "log10(mol/L)", 100.0, -3.0),
        (180.16, "mg/L", "log10(mol/L)", 180.16, -3.0),
        (373.15, "K", "degC", None, 100.0),
        (212.0, "degF", "degC", None, 100.0),
        (25.0, "degC", "degC", None, 25.0),
    ],
)
def test_convert_unit(value, unit_from, unit_to, mw, expected):
    assert convert_unit(value, unit_from, unit_to, mw=mw) == pytest.approx(expected, abs=1e-12)


def test_convert_unit_errors_and_nan():
    assert math.isnan(convert_unit(float("nan"), "nM", "log10(mol/L)"))
    with pytest.raises(ValueError, match="molar mass"):
        convert_unit(1.0, "ug/mL", "log10(mol/L)")
    with pytest.raises(ValueError, match="molar mass"):
        convert_unit(1.0, "log10(ug/mL)", "log10(mol/L)", mw=float("nan"))
    with pytest.raises(ValueError, match="non-positive"):
        convert_unit(0.0, "nM", "log10(mol/L)")
    with pytest.raises(ValueError):
        convert_unit(1.0, "furlong", "log10(mol/L)")
    with pytest.raises(ValueError):
        convert_unit(1.0, "nM", "degC")


# ----------------------------------------------------------------------------- aggregation
def _recs(rows, source="aqsoldb"):
    base = dict(source=source, record_id="", smiles="", raw_value="", raw_unit="", n_reported=1,
                sd_reported=float("nan"), canonical_smiles="", status="ok", detail="")
    df = pd.DataFrame([{**base, **r} for r in rows])
    return df[list(RECORD_COLUMNS)]


def test_aggregate_median_counts_and_sd():
    recs = _recs(
        [
            dict(inchikey="BBB", value=1.0),
            dict(inchikey="AAA", value=-3.0, n_reported=3, sd_reported=0.2),
            dict(inchikey="BBB", value=2.0),
            dict(inchikey="BBB", value=10.0),
            dict(inchikey="CCC", value=float("nan"), status="unparsable"),
        ]
    )
    agg = aggregate(recs, LOGS)
    assert tuple(agg.columns) == TARGET_COLUMNS
    assert agg["inchikey"].tolist() == ["AAA", "BBB"]  # sorted; rejected record absent
    a, b = agg.iloc[0], agg.iloc[1]
    assert a["value"] == -3.0 and a["n_measurements"] == 3  # source aggregate counts its occurrences
    assert a["sd"] == pytest.approx(0.2) and not a["replicates_disagree"]
    assert b["value"] == 2.0  # median, not mean (4.33)
    assert b["n_measurements"] == 3
    assert b["sd"] == pytest.approx(np.std([1.0, 2.0, 10.0], ddof=1))
    assert b["replicates_disagree"]


def test_aggregate_even_count_median_and_threshold_is_strict():
    # two values 0.5*sqrt(2) apart -> sample SD exactly 0.5 -> not flagged (rule is sd > 0.5)
    d = 0.5 * math.sqrt(2.0)
    agg = aggregate(_recs([dict(inchikey="K", value=0.0), dict(inchikey="K", value=d)]), LOGS)
    assert agg.loc[0, "value"] == pytest.approx(d / 2)
    assert agg.loc[0, "sd"] == pytest.approx(0.5)
    assert not agg.loc[0, "replicates_disagree"]
    agg = aggregate(_recs([dict(inchikey="K", value=0.0), dict(inchikey="K", value=d + 1e-6)]), LOGS)
    assert agg.loc[0, "replicates_disagree"]


def test_aggregate_source_reported_sd_flags_single_record():
    # an AqSolDB row of group G2/G4 (its own SD > 0.5) is flagged even alone
    agg = aggregate(_recs([dict(inchikey="K", value=-1.0, n_reported=2, sd_reported=0.9)]), LOGS)
    assert agg.loc[0, "replicates_disagree"] and agg.loc[0, "n_measurements"] == 2


def test_aggregate_melting_point_threshold():
    rows = [dict(inchikey="K", value=100.0), dict(inchikey="K", value=112.0)]  # SD 8.49 degC
    assert not aggregate(_recs(rows, "bradley"), MP).loc[0, "replicates_disagree"]
    rows[1]["value"] = 116.0  # SD 11.3 degC > 10
    assert aggregate(_recs(rows, "bradley"), MP).loc[0, "replicates_disagree"]


def test_aggregate_empty():
    agg = aggregate(_recs([dict(inchikey="K", value=float("nan"), status="mixture")]), LOGS)
    assert len(agg) == 0 and tuple(agg.columns) == TARGET_COLUMNS


# ----------------------------------------------------------------------------- records
def _raw(rows):
    base = dict(record_id="r", raw_value="1", raw_unit="u", value=1.0, n_reported=1,
                sd_reported=float("nan"), status="", detail="")
    return pd.DataFrame([{**base, **r} for r in rows])[list(ed._RAW)]


def test_process_records_rejections_are_recorded_with_reasons():
    raw = _raw(
        [
            dict(record_id="ok", smiles="c1ccccc1O"),
            dict(record_id="bad", smiles="not_a_smiles"),
            dict(record_id="hg", smiles="C[Hg]C"),
            dict(record_id="tiny", smiles="CO"),
            dict(record_id="nocarbon", smiles="O=S(=O)(O)O"),
            dict(record_id="mix", smiles="c1ccccc1.CCCCCCO"),
            dict(record_id="empty", smiles=""),
            dict(record_id="flag", smiles="CCCCO", status="flagged_by_source", detail="donotuse: salt"),
        ]
    )
    df = process_records(raw, LOGP)
    st = dict(zip(df["record_id"], df["status"]))
    assert st == {
        "ok": "ok",
        "bad": "unparsable",
        "hg": "disallowed_element",
        "tiny": "too_few_heavy_atoms",
        "nocarbon": "no_carbon",
        "mix": "mixture",
        "empty": "missing_structure",
        "flag": "flagged_by_source",
    }
    assert set(df["status"]) <= set(STATUS_CODES)
    bad = df[df["status"] != "ok"]
    assert bad["value"].isna().all()  # never imputed
    assert (bad.loc[bad["status"] != "missing_structure", "detail"] != "").all()
    ok = df[df["status"] == "ok"].iloc[0]
    assert ok["inchikey"] == standardize("Oc1ccccc1").inchikey == "ISWSIDIOOBJBQZ-UHFFFAOYSA-N"


def test_salt_forms_excluded_for_solid_state_targets_only():
    raw = _raw(
        [
            dict(record_id="free", smiles="CN(C)C(=N)N=C(N)N", value=-0.5),
            dict(record_id="hcl", smiles="CN(C)C(=N)N=C(N)N.Cl", value=0.2),
            dict(record_id="hydrate", smiles="CC(=O)Nc1ccc(O)cc1.O", value=-1.0),
        ]
    )
    solid = process_records(raw, LOGS)
    assert solid["status"].tolist() == ["ok", "salt_or_solvate_form", "salt_or_solvate_form"]
    assert "2 components" in solid.loc[1, "detail"]
    # the parent of the HCl salt is the free base: same key, but its value is not used
    assert solid.loc[1, "inchikey"] == solid.loc[0, "inchikey"]
    agg = aggregate(solid, LOGS)
    assert agg["n_measurements"].tolist() == [1]
    # logP is a property of the neutral species: salts collapse onto the parent
    liquid = process_records(raw, LOGP)
    assert liquid["status"].tolist() == ["ok", "ok", "ok"]
    # sensitivity switch
    kept = process_records(raw, LOGS, keep_salt_forms=True)
    assert kept["status"].tolist() == ["ok", "ok", "ok"]
    agg = aggregate(kept, LOGS)
    met = agg[agg["inchikey"] == solid.loc[0, "inchikey"]].iloc[0]
    assert met["n_measurements"] == 2 and met["value"] == pytest.approx(-0.15)


def test_process_records_independent_of_n_jobs(monkeypatch):
    df = pd.read_csv(FIX / "bradley_excerpt.csv")
    raw = parse_bradley(df)
    serial = process_records(raw, MP, n_jobs=1)
    monkeypatch.setattr(ed, "_POOL_MIN", 0)
    monkeypatch.setattr(ed, "_STD_MEMO", {})  # force the pool to recompute
    pooled = process_records(raw, MP, n_jobs=2)
    pd.testing.assert_frame_equal(serial, pooled, check_exact=True)


# ----------------------------------------------------------------------------- parsers
def test_parse_aqsoldb_excerpt():
    df = pd.read_csv(FIX / "aqsoldb_excerpt.csv", dtype={"ID": str})
    raw = parse_aqsoldb(df)
    assert len(raw) == len(df)
    assert raw["raw_unit"].eq("log10(mol/L)").all()
    assert raw["n_reported"].tolist() == df["Ocurrences"].astype(int).tolist()
    assert np.allclose(raw["value"], df["Solubility"])
    recs = process_records(raw, LOGS)
    assert set(recs["status"]) <= set(STATUS_CODES)
    salts = df["SMILES"].str.contains(r"\.", regex=True)
    assert (recs.loc[salts.values, "status"] != "ok").all()


def test_parse_bradley_excerpt_marks_donotuse():
    df = pd.read_csv(FIX / "bradley_excerpt.csv")
    raw = parse_bradley(df)
    flagged = df["donotuse"].notna().values
    assert flagged.any()
    assert (raw.loc[flagged, "status"] == "flagged_by_source").all()
    assert raw.loc[flagged, "detail"].str.startswith("donotuse").all()
    assert (raw.loc[~flagged, "status"] == "").all()
    recs = process_records(raw, MP)
    agg = aggregate(recs, MP)
    # replicated compounds in the excerpt are aggregated by median
    assert (agg["n_measurements"] > 1).any()


SDF = """
  synthetic

  0  0  0  0  0  0  0  0  0  0999 V2000
M  END
> <ChemID>
1

> <Original_SMILES>
CCCCCCO

> <Canonical_QSARr>
CCCCCCO

> <LogP>
2.03

$$$$
x
M  END
> <ChemID>
2

> <Original_SMILES>


> <Canonical_QSARr>
OC(=O)c1ccccc1

> <LogP>
1.87

$$$$
x
M  END
> <ChemID>
3

> <Original_SMILES>
CCCCCCO

> <LogP>
2.2

$$$$
x
M  END
> <ChemID>
4

> <Original_SMILES>
CCCC

> <LogP>
n/a

$$$$
x
M  END
> <ChemID>
5

> <Original_SMILES>
C[C@H](N)CCCC|C[C@@H](N)CCCC

> <LogP>
1.5

$$$$
"""


def test_read_sdf_and_parse_opera():
    props = read_sdf_properties(SDF)
    assert props["ChemID"].tolist() == ["1", "2", "3", "4", "5"]
    raw = parse_opera(props, "LogP_QR.sdf")
    assert raw["smiles"].tolist()[:4] == ["CCCCCCO", "OC(=O)c1ccccc1", "CCCCCCO", "CCCC"]  # QSAR-ready fallback
    assert raw["status"].tolist() == ["", "", "", "missing_value", "ambiguous_structure"]
    assert raw.loc[4, "detail"].startswith("2 alternative")
    agg = aggregate(process_records(raw, LOGP), LOGP)
    hexanol = agg[agg["inchikey"] == standardize("CCCCCCO").inchikey].iloc[0]
    assert hexanol["n_measurements"] == 2 and hexanol["value"] == pytest.approx(2.115)
    assert len(agg) == 2


def test_parse_opera_counts_merged_physprop_ids():
    """An OPERA record merging several PHYSPROP entries (MPID ``a|b|c``) counts as that many
    measurements, like an AqSolDB row with Ocurrences > 1 (regression: was always 1)."""
    props = pd.DataFrame(
        {
            "ChemID": ["1", "2", "3", "4"],
            "Original_SMILES": ["CCCCCCO", "CCCCCCO", "OC(=O)c1ccccc1", "CCCCO"],
            "MPID": ["303|361286|362587", "17", "", "5| 5 |"],
            "LogP": ["2.0", "2.2", "1.9", "0.9"],
        }
    )
    raw = parse_opera(props, "LogP_QR.sdf")
    assert raw["n_reported"].tolist() == [3, 1, 1, 1]  # blank -> 1; duplicates/blanks ignored
    assert raw["sd_reported"].isna().all()  # OPERA publishes no spread
    agg = aggregate(process_records(raw, LOGP), LOGP)
    hexanol = agg[agg["inchikey"] == standardize("CCCCCCO").inchikey].iloc[0]
    assert hexanol["n_measurements"] == 4


def test_assay_limit_values_are_censored():
    """AstraZeneca values at (or beyond) the assay range limits are limit values, not
    measurements, although ChEMBL stores them with '=' (regression)."""
    logd = TARGETS["logd74_az"]
    assert logd.assay_range == (-1.5, 4.5)
    vals = ["-1.5", "-1.49", "2.0", "4.49", "4.5", "4.7"]
    acts = [dict(activity_id=i, canonical_smiles="CCCCCCO", standard_relation="=", standard_value=v)
            for i, v in enumerate(vals, 1)]
    recs = process_records(parse_chembl_activities(acts, "log10"), logd)
    assert recs["status"].tolist() == ["censored_value", "ok", "ok", "ok", "censored_value", "censored_value"]
    assert recs.loc[0, "detail"].startswith("-1.5 at or beyond the assay range")
    assert recs.loc[recs["status"] != "ok", "value"].isna().all()
    sol = TARGETS["logs74_az"]
    acts = [dict(activity_id=i, canonical_smiles="CCCCCCO", standard_relation="=", standard_value=v,
                 standard_units="nM") for i, v in enumerate(["100", "101", "1513600", "1500000"], 1)]
    recs = process_records(parse_chembl_activities(acts, sol.unit), sol)
    assert recs["status"].tolist() == ["censored_value", "ok", "censored_value", "ok"]
    # targets without an assay range are unaffected
    assert all(TARGETS[t].assay_range is None for t in TARGETS if TARGETS[t].source not in ("az_logd", "az_solubility"))


def test_parse_biogen_converts_with_recorded_mw():
    df = pd.DataFrame(
        {
            "Internal ID": ["Mol1", "Mol2", "Mol3"],
            "SMILES": ["c1ccccc1O", "c1ccccc1O", "xx"],
            "LOG SOLUBILITY PH 6.8 (ug/mL)": [2.0, float("nan"), 1.0],
        }
    )
    raw = parse_biogen(df, "LOG SOLUBILITY PH 6.8 (ug/mL)")
    mw = 94.113
    assert raw.loc[0, "value"] == pytest.approx(2.0 - 3.0 - math.log10(mw), abs=1e-4)
    assert raw["status"].tolist() == ["", "missing_value", "unit_conversion_failed"]
    assert math.isnan(raw.loc[2, "value"])


def test_parse_chembl_activities_censoring_flags_and_units():
    acts = [
        dict(activity_id=3, canonical_smiles="CCCCCCO", standard_relation="=", standard_value="1000.0",
             standard_units="nM", data_validity_comment=None),
        dict(activity_id=1, canonical_smiles="CCCCCCO", standard_relation=">", standard_value="1500000",
             standard_units="nM", data_validity_comment=None),
        dict(activity_id=2, canonical_smiles="CCCCCCO", standard_relation=None, standard_value="50",
             standard_units="nM", data_validity_comment="Outside typical range"),
        dict(activity_id=4, canonical_smiles="CCCCCCO", standard_relation=None, standard_value=None,
             standard_units="nM", data_validity_comment=None),
        dict(activity_id=5, canonical_smiles="CCCCCCO", standard_relation=None, standard_value="7",
             standard_units="%", data_validity_comment=None),
    ]
    raw = parse_chembl_activities(acts, "log10(mol/L)")
    assert raw["record_id"].tolist() == ["1", "2", "3", "4", "5"]  # sorted by activity_id
    assert raw["status"].tolist() == [
        "censored_value", "flagged_by_source", "", "missing_value", "unit_conversion_failed"
    ]
    assert raw.loc[2, "value"] == pytest.approx(-6.0)
    assert raw.loc[raw["status"] != "", "value"].isna().all()
    logd = parse_chembl_activities(
        [dict(activity_id=9, canonical_smiles="CCO", standard_value="1.8", standard_units=None)], "log10"
    )
    assert logd.loc[0, "value"] == 1.8 and logd.loc[0, "status"] == ""


class _FakeChembl:
    release = "ChEMBL_37"
    offline = True

    def __init__(self, records):
        self.records = records

    def fetch_all(self, endpoint, params):
        assert endpoint == "activity"
        return self.records


def test_chembl_source_hash_is_checked(tmp_path):
    recs = [dict(activity_id=1, canonical_smiles="CCCCCCO", standard_value="1.0")]
    with pytest.raises(ChecksumError, match="expected"):
        ed.load_records("logd74_az", tmp_path, client=_FakeChembl(recs))
    other = _FakeChembl(recs)
    other.release = "ChEMBL_99"
    with pytest.raises(ChecksumError, match="no pinned sha256"):
        ed.load_records("logd74_az", tmp_path, client=other)


def test_chembl_source_pinned_records_load(tmp_path, monkeypatch):
    recs = [
        dict(activity_id=2, canonical_smiles="CCCCCCO", standard_value="1.9"),
        dict(activity_id=1, canonical_smiles="CCCCCCO", standard_value="1.7"),
    ]
    src = dataclasses.replace(SOURCES["az_logd"], chembl_sha256={"ChEMBL_37": ed.chembl_records_sha256(recs)})
    monkeypatch.setitem(SOURCES, "az_logd", src)
    df = get_target("logd74_az", tmp_path, client=_FakeChembl(recs))
    assert df["value"].tolist() == [pytest.approx(1.8)]
    assert df.attrs["raw_sha256"].startswith("ChEMBL_37:")


# ----------------------------------------------------------------------------- download
class _Resp:
    def __init__(self, content):
        self.content = content

    def raise_for_status(self):
        pass


class _Session:
    def __init__(self, content):
        self.content = content
        self.calls = 0

    def get(self, url, timeout=None):
        self.calls += 1
        return _Resp(self.content)


def _src(data: bytes):
    return dataclasses.replace(SOURCES["aqsoldb"], sha256=sha256_bytes(data), filename="x.csv")


def test_fetch_downloads_verifies_and_caches(tmp_path):
    data = b"ID,SMILES\n"
    src = _src(data)
    sess = _Session(data)
    p = fetch(src, tmp_path, session=sess)
    assert p == tmp_path / "aqsoldb" / "x.csv" and p.read_bytes() == data
    assert fetch(src, tmp_path, session=sess) == p and sess.calls == 1  # served from cache
    assert fetch(src, tmp_path, offline=True) == p


def test_fetch_checksum_mismatch_is_not_stored(tmp_path):
    src = _src(b"expected content")
    with pytest.raises(ChecksumError, match="not stored"):
        fetch(src, tmp_path, session=_Session(b"tampered content"))
    assert not (tmp_path / "aqsoldb" / "x.csv").exists()
    assert list((tmp_path).rglob("*")) == []


def test_fetch_corrupted_cache_raises(tmp_path):
    src = _src(b"expected content")
    p = tmp_path / "aqsoldb" / "x.csv"
    p.parent.mkdir(parents=True)
    p.write_bytes(b"corrupted")
    with pytest.raises(ChecksumError, match="delete the file"):
        fetch(src, tmp_path, session=_Session(b"expected content"))
    assert p.read_bytes() == b"corrupted"  # not silently replaced


def test_fetch_offline_cache_miss(tmp_path):
    with pytest.raises(CacheMiss):
        fetch(_src(b"x"), tmp_path, offline=True)
    with pytest.raises(ValueError, match="not a file source"):
        fetch("az_logd", tmp_path, offline=True)


def test_registry_is_complete():
    for name, spec in TARGETS.items():
        src = SOURCES[spec.source]
        assert spec.name == name
        assert src.licence and src.doi and src.version
        assert bool(src.url and src.sha256) != bool(src.chembl_assay)
        assert spec.disagreement_sd > 0
        if src.key == "opera":
            assert spec.member in ed._OPERA_VALUE
    assert len({s.sha256 for s in SOURCES.values() if s.sha256}) == 4
    for spec in TARGETS.values():
        if spec.assay_range is not None:
            assert spec.assay_range[0] < spec.assay_range[1]


# ----------------------------------------------------------------------------- end to end
@pytest.fixture()
def aqsol_cache(tmp_path, monkeypatch):
    data = (FIX / "aqsoldb_excerpt.csv").read_bytes()
    src = dataclasses.replace(SOURCES["aqsoldb"], sha256=sha256_bytes(data))
    monkeypatch.setitem(SOURCES, "aqsoldb", src)
    p = tmp_path / "aqsoldb" / src.filename
    p.parent.mkdir(parents=True)
    p.write_bytes(data)
    return tmp_path


def test_get_target_offline_deterministic(aqsol_cache):
    a = get_target("logs_aqsoldb", aqsol_cache, offline=True)  # processes and caches
    b = get_target("logs_aqsoldb", aqsol_cache, offline=True)  # from the processed cache
    pd.testing.assert_frame_equal(a, b, check_exact=True)
    assert tuple(a.columns) == TARGET_COLUMNS
    assert a["inchikey"].is_monotonic_increasing and a["inchikey"].is_unique
    assert a["inchikey"].str.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]").all()
    assert not a["replicates_disagree"].any()
    assert a.attrs["unit"] == "log10(mol/L)" and a.attrs["source"] == "aqsoldb"
    full = get_target("logs_aqsoldb", aqsol_cache, offline=True, drop_disagreeing=False)
    assert len(full) - len(a) == full["replicates_disagree"].sum() == a.attrs["n_dropped_disagreeing"] > 0
    salts = get_target("logs_aqsoldb", aqsol_cache, offline=True, drop_disagreeing=False, keep_salt_forms=True)
    assert len(salts) > len(full)
    recs = ed.load_records("logs_aqsoldb", aqsol_cache, offline=True)  # cached
    fresh = process_records(parse_aqsoldb(pd.read_csv(FIX / "aqsoldb_excerpt.csv", dtype={"ID": str})), LOGS)
    pd.testing.assert_frame_equal(recs, fresh, check_exact=True)
    assert len(recs) == len(pd.read_csv(FIX / "aqsoldb_excerpt.csv"))
    assert recs.loc[recs["status"] != "ok", "value"].isna().all()


def test_coverage_full_key_and_first_block():
    t = {"x": pd.DataFrame({"inchikey": ["AAAAAAAAAAAAAA-BBBBBBBBBB-N", "CCCCCCCCCCCCCC-DDDDDDDDDD-N"]})}
    keys = ["AAAAAAAAAAAAAA-BBBBBBBBBB-N", "CCCCCCCCCCCCCC-EEEEEEEEEE-N", "ZZZZZZZZZZZZZZ-BBBBBBBBBB-N", ""]
    assert coverage(t, keys) == {"x": {"n": 1, "n_first_block": 2}}
    # two drugs (stereoisomers) sharing a connectivity block count twice: the column counts
    # drugs, not distinct blocks (regression)
    keys2 = ["AAAAAAAAAAAAAA-FFFFFFFFFF-N", "AAAAAAAAAAAAAA-GGGGGGGGGG-N", "AAAAAAAAAAAAAA-FFFFFFFFFF-N"]
    assert coverage(t, keys2) == {"x": {"n": 0, "n_first_block": 2}}


def _reject_all(smiles, **kw):
    from topoqspr.graph import Rejection, StandardizationResult

    return StandardizationResult(input_smiles=smiles, rejection=Rejection("unparsable", "probe"))


def test_processed_cache_invalidated_when_standardize_changes(aqsol_cache, monkeypatch):
    """Regression: the processed-records cache key ignored graph.standardize, so a changed
    standardiser silently reused stale InChIKeys/statuses."""
    a = get_target("logs_aqsoldb", aqsol_cache, offline=True)
    assert len(a) > 0
    fp = a.attrs["processing_fingerprint"]
    monkeypatch.setattr(ed, "standardize", _reject_all)
    b = get_target("logs_aqsoldb", aqsol_cache, offline=True)
    assert len(b) == 0 and b.attrs["processing_fingerprint"] != fp  # not served from the cache
    recs = ed.load_records("logs_aqsoldb", aqsol_cache, offline=True)
    assert (recs["status"] == "unparsable").sum() == recs["smiles"].ne("").sum()
    monkeypatch.setattr(ed, "standardize", standardize)
    c = get_target("logs_aqsoldb", aqsol_cache, offline=True)  # memo not stale either
    pd.testing.assert_frame_equal(a, c, check_exact=True)
    assert c.attrs["processing_fingerprint"] == fp
    assert len(list((aqsol_cache / "processed").glob("logs_aqsoldb-*.csv.gz"))) == 2


def test_processing_fingerprint_tracks_graph_source(tmp_path, monkeypatch):
    from topoqspr import graph

    fp = ed.processing_fingerprint()
    assert ed.processing_fingerprint() == fp  # deterministic
    edited = tmp_path / "graph.py"
    edited.write_text(Path(graph.__file__).read_text() + "\n# edited\n")
    monkeypatch.setattr(graph, "__file__", str(edited))
    assert ed.processing_fingerprint() != fp


class _PoolClient:
    """Minimal stand-in for ChemblClient used by approved_drug_pool."""

    release = "ChEMBL_37"

    def fetch_all(self, endpoint, params):
        assert endpoint == "drug_indication" and params["max_phase_for_ind"] == "4"
        return [
            dict(drugind_id=1, molecule_chembl_id="CHEMBL10", parent_molecule_chembl_id="CHEMBL1",
                 max_phase_for_ind="4.0", mesh_id="D1"),
            dict(drugind_id=2, molecule_chembl_id="CHEMBL1", parent_molecule_chembl_id="CHEMBL1",
                 max_phase_for_ind="4.0", mesh_id="D2"),
            dict(drugind_id=3, molecule_chembl_id="CHEMBL2", parent_molecule_chembl_id="CHEMBL2",
                 max_phase_for_ind="4.0", mesh_id="D1"),
            dict(drugind_id=4, molecule_chembl_id="CHEMBL3", parent_molecule_chembl_id="CHEMBL3",
                 max_phase_for_ind="4.0", mesh_id="D3"),
        ]

    def molecules(self, ids):
        def mol(cid, smi, parent, mtype="Small molecule", key="K"):
            return dict(molecule_chembl_id=cid, pref_name=cid, molecule_type=mtype, structure_type="MOL",
                        inorganic_flag=0, max_phase="4.0", atc_classifications=[],
                        molecule_hierarchy={"parent_chembl_id": parent},
                        molecule_structures={"canonical_smiles": smi, "standard_inchi_key": key})
        db = {
            "CHEMBL1": mol("CHEMBL1", "CN(C)C(=N)N=C(N)N", "CHEMBL1", key="XZWYZXLIPXDOLR-UHFFFAOYSA-N"),
            "CHEMBL10": mol("CHEMBL10", "CN(C)C(=N)N=C(N)N.Cl", "CHEMBL1", key="OETHQSJEHLVLGH-UHFFFAOYSA-N"),
            "CHEMBL2": mol("CHEMBL2", "CC(=O)Nc1ccc(O)cc1", "CHEMBL2", key="RZVAJINKPMORJF-UHFFFAOYSA-N"),
            "CHEMBL3": mol("CHEMBL3", "", "CHEMBL3", mtype="Protein", key=""),
        }
        return {i: db.get(i) for i in ids}

    def verify_release(self):
        pass

    def status(self):
        return {"chembl_db_version": self.release, "chembl_release_date": "2026-05-01"}


def test_approved_drug_pool_and_standardised_keys():
    pool = ed.approved_drug_pool(_PoolClient())
    inc = pool.included
    assert inc["parent_chembl_id"].tolist() == ["CHEMBL1", "CHEMBL2"]
    reasons = dict(zip(pool.manifest["molecule_chembl_id"], pool.manifest["reason"]))
    assert reasons["CHEMBL10"] == "collapsed_to_parent" and reasons["CHEMBL3"] == "not_small_molecule"
    keys = ed.standardised_keys(inc)
    assert keys["inchikey"].tolist() == ["XZWYZXLIPXDOLR-UHFFFAOYSA-N", "RZVAJINKPMORJF-UHFFFAOYSA-N"]
    assert pool.extra == {"indication_rows": 4, "n_mesh_ids": 3}


# ----------------------------------------------------------------------------- live
@pytest.mark.network
def test_live_aqsoldb_download_matches_pin(tmp_path):
    p = fetch("aqsoldb", tmp_path)
    df = pd.read_csv(p)
    assert len(df) == 9982
