"""Tests for topoqspr.chembl and topoqspr.cohort.

Default tests run offline from the recorded caches in tests/fixtures/cohort/ (ChEMBL_37):
Gout (D006073, 39 indication rows), Diabetes Mellitus, Type 2 (D003924, 376 rows), the
naive v0.x retrieval for "gout", and ATC M04 (antigout preparations).  Live tests are
marked ``network``.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import requests

from topoqspr.chembl import CacheMiss, ChemblClient, MeshClient
from topoqspr.cohort import (
    MANIFEST_COLUMNS,
    REASONS,
    build_atc_cohort,
    build_cohort,
    jaccard,
    naive_cohort,
    structure_reason,
    suggest_headings,
)

FIX = Path(__file__).resolve().parent / "fixtures" / "cohort"
SRC = Path(__file__).resolve().parents[1] / "src"
GOUT, ARTHRITIS_GOUTY, T2D = "D006073", "D015210", "D003924"
RELEASE = "ChEMBL_37"


@pytest.fixture(scope="module")
def client() -> ChemblClient:
    return ChemblClient(FIX / "chembl", offline=True)


@pytest.fixture(scope="module")
def mesh() -> MeshClient:
    return MeshClient(FIX / "mesh", offline=True)


@pytest.fixture(scope="module")
def gout(client):
    return build_cohort([GOUT], client)


@pytest.fixture(scope="module")
def t2d(client):
    return build_cohort([T2D], client)


def _row(cohort, chembl_id):
    m = cohort.manifest
    hit = m[m["molecule_chembl_id"] == chembl_id]
    assert len(hit) == 1, chembl_id
    return hit.iloc[0]


# ----------------------------------------------------------------------------- client
def test_offline_release_from_fixture(client):
    assert client.release == RELEASE
    assert client.status()["chembl_release_date"] == "2026-05-01"


def test_offline_cache_miss_raises(client):
    with pytest.raises(CacheMiss):
        client.drug_indications(["D000001"])
    with pytest.raises(CacheMiss):
        client.molecules(["CHEMBL999999999"])


def test_cache_keys_never_mix_releases(tmp_path):
    for rel in ("ChEMBL_36", "ChEMBL_37"):
        (tmp_path / rel).mkdir()
        (tmp_path / rel / "status.json").write_text(json.dumps({"chembl_db_version": rel}))
    a = ChemblClient(tmp_path, offline=True, release="ChEMBL_36")
    b = ChemblClient(tmp_path, offline=True, release="ChEMBL_37")
    pa, ka = a._cache_path("drug_indication", {"mesh_id": GOUT})
    pb, kb = b._cache_path("drug_indication", {"mesh_id": GOUT})
    assert pa != pb and ka != kb
    assert pa.parent.parent.name == "ChEMBL_36" and pb.parent.parent.name == "ChEMBL_37"
    # an entry recorded for 37 is invisible to a 36 client
    shutil.copytree(FIX / "chembl" / RELEASE, tmp_path / "ChEMBL_37", dirs_exist_ok=True)
    assert len(b.drug_indications([GOUT])) == 39
    with pytest.raises(CacheMiss):
        a.drug_indications([GOUT])


def test_cache_file_key_verified(tmp_path, client):
    shutil.copytree(FIX / "chembl", tmp_path / "chembl")
    c = ChemblClient(tmp_path / "chembl", offline=True)
    path, _ = c._cache_path("drug_indication", {"mesh_id": GOUT, "order_by": "drugind_id"})
    other, _ = c._cache_path("drug_indication", {"mesh_id": T2D, "order_by": "drugind_id"})
    shutil.copy(other, path)  # wrong content under this key
    with pytest.raises(RuntimeError, match="cache key mismatch"):
        c.drug_indications([GOUT])


class _Resp:
    def __init__(self, status, payload=None, headers=None):
        self.status_code = status
        self._payload = payload
        self.headers = headers or {}
        self.url = "http://test"

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")


class _Session:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls = []

    def get(self, url, params=None, timeout=None, headers=None):
        self.calls.append((url, dict(params or {})))
        r = self.responses.pop(0)
        if isinstance(r, Exception):
            raise r
        return r


def _online(tmp_path, responses, **kw):
    status = _Resp(200, {"chembl_db_version": "ChEMBL_99", "chembl_release_date": "2099-01-01"})
    sess = _Session([status, *responses])
    sleeps = []
    c = ChemblClient(tmp_path, session=sess, sleep=sleeps.append, backoff=0.5, **kw)
    return c, sess, sleeps


def test_retry_with_exponential_backoff(tmp_path):
    page = {"drug_indications": [{"drugind_id": 1}], "page_meta": {"total_count": 1, "next": None}}
    c, sess, sleeps = _online(
        tmp_path,
        [_Resp(503), requests.ConnectionError("boom"), _Resp(429, headers={"Retry-After": "3"}), _Resp(200, page)],
    )
    assert c.fetch_all("drug_indication", {"mesh_id": "D1"}) == [{"drugind_id": 1}]
    assert sleeps == [0.5, 1.0, 3.0]  # backoff * 2**attempt, Retry-After honoured
    # second call is served from the cache: no new HTTP request
    n = len(sess.calls)
    assert c.fetch_all("drug_indication", {"mesh_id": "D1"}) == [{"drugind_id": 1}]
    assert len(sess.calls) == n
    assert (tmp_path / "ChEMBL_99" / "status.json").exists()


def test_non_retryable_error_raises_immediately(tmp_path):
    c, sess, sleeps = _online(tmp_path, [_Resp(404)])
    with pytest.raises(requests.HTTPError):
        c.fetch_all("drug_indication", {"mesh_id": "D1"})
    assert sleeps == [] and len(sess.calls) == 2


def test_retries_exhausted(tmp_path):
    c, _, sleeps = _online(tmp_path, [_Resp(500)] * 3, max_retries=2)
    with pytest.raises(RuntimeError, match="failed after 3 attempts"):
        c.fetch_all("drug_indication", {"mesh_id": "D1"})
    assert sleeps == [0.5, 1.0]
    assert not list((tmp_path / "ChEMBL_99").glob("drug_indication/*"))  # failures never cached


def test_pagination_collects_all_pages(tmp_path):
    def page(ids, total, nxt):
        return _Resp(200, {"molecules": [{"i": i} for i in ids], "page_meta": {"total_count": total, "next": nxt}})

    c, sess, _ = _online(tmp_path, [page([1, 2], 5, "n"), page([3, 4], 5, "n"), page([5], 5, None)], page_limit=2)
    assert [r["i"] for r in c.fetch_all("molecule", {"x": "y"})] == [1, 2, 3, 4, 5]
    assert [p["offset"] for _, p in sess.calls[1:]] == [0, 2, 4]
    assert all(p["limit"] == 2 for _, p in sess.calls[1:])


def test_pagination_detects_short_result(tmp_path):
    resp = _Resp(200, {"molecules": [{"i": 1}], "page_meta": {"total_count": 3, "next": None}})
    c, _, _ = _online(tmp_path, [resp])
    with pytest.raises(RuntimeError, match="expected 3"):
        c.fetch_all("molecule", {"x": "y"})


def test_invalid_json_is_retried(tmp_path):
    class _Bad(_Resp):
        def json(self):
            raise ValueError("not json")

    page = {"molecules": [], "page_meta": {"total_count": 0, "next": None}}
    c, _, sleeps = _online(tmp_path, [_Bad(200), _Resp(200, page)])
    assert c.fetch_all("molecule", {"x": "y"}) == []
    assert sleeps == [0.5]


def test_verify_release_detects_switch(tmp_path):
    c, _, _ = _online(tmp_path, [_Resp(200, {"chembl_db_version": "ChEMBL_100"})])
    assert c.release == "ChEMBL_99"
    with pytest.raises(RuntimeError, match="changed during the run"):
        c.verify_release()


def test_sparql_pages_until_short_page(tmp_path):
    def page(n):
        return _Resp(200, {"results": {"bindings": [{"d": {"value": f"x{i}"}} for i in range(n)]}})

    sess = _Session([page(2), page(2), page(1)])
    m = MeshClient(tmp_path, session=sess)
    assert len(m._sparql("SELECT ?d WHERE {} ORDER BY ?d", page=2)) == 5
    assert [p["offset"] for _, p in sess.calls] == [0, 2, 4]


def test_release_pin_mismatch_raises(tmp_path):
    sess = _Session([_Resp(200, {"chembl_db_version": "ChEMBL_38"})])
    c = ChemblClient(tmp_path, session=sess, release="ChEMBL_37")
    with pytest.raises(RuntimeError, match="requested ChEMBL release"):
        c.release


# ----------------------------------------------------------------------------- structure rules
@pytest.mark.parametrize(
    "smiles, reason",
    [
        (None, "no_smiles"),
        ("", "no_smiles"),
        ("C1CC", "unparseable_smiles"),
        ("c1cccc1", "unparseable_smiles"),
        ("[Na+].[Cl-]", "no_carbon"),
        ("O=[Mg]", "no_carbon"),
        ("CCCC(=O)[O-].[Na+]", "metal_containing"),
        ("O=C([O-])c1cccnc1.[Cr+3]", "metal_containing"),
        ("C[As](C)C", "unsupported_element"),
        ("*CC", "unsupported_element"),
        ("CCO.CCN", "mixture"),
        ("CCO", None),
        ("OB(O)c1ccccc1", None),  # boron allowed (bortezomib-like)
        ("O=C1c2ccccc2[se]n1-c1ccccc1", None),  # selenium allowed (ebselen)
        ("N[C@@H](CC(=O)N1CCn2c(nnc2C(F)(F)F)C1)Cc1cc(F)c(F)cc1F", None),  # sitagliptin
    ],
)
def test_structure_reason(smiles, reason):
    assert structure_reason(smiles) == reason


# ----------------------------------------------------------------------------- controlled cohort
def test_manifest_columns_and_reasons(gout, t2d):
    for co in (gout, t2d):
        assert tuple(co.manifest.columns) == MANIFEST_COLUMNS
        codes = set(co.manifest["reason"].str.split(":").str[0])
        assert codes <= set(REASONS)
        assert set(co.manifest["decision"]) <= {"include", "exclude"}
        assert ((co.manifest["decision"] == "include") == (co.manifest["reason"] == "included")).all()
        assert (co.manifest["chembl_release"] == RELEASE).all()


def test_gout_matches_recorded_manifest(gout):
    assert gout.to_csv_bytes() == (FIX / "expected_gout.csv").read_bytes()


def test_gout_cohort_content(gout):
    inc = gout.included
    assert list(inc["pref_name"]) == [
        "COLCHICINE", "NAPROXEN", "SULFINPYRAZONE", "PROBENECID", "ALLOPURINOL",
        "TISOPURINE", "CINCHOPHEN", "BENZBROMARONE", "FEBUXOSTAT", "LESINURAD",
    ]
    assert (inc["phase_for_indication"] == 4.0).all()
    assert gout.reason_counts() == {
        "collapsed_to_parent": 2, "included": 10, "not_small_molecule": 2, "phase_below_min": 26,
    }


def test_every_candidate_appears_once(client, t2d):
    rows = client.drug_indications([T2D])
    candidates = {r["molecule_chembl_id"] for r in rows}
    m = t2d.manifest
    assert candidates <= set(m["molecule_chembl_id"])
    assert m["molecule_chembl_id"].is_unique
    extra = set(m["molecule_chembl_id"]) - candidates  # parents without their own row
    assert extra <= set(m["parent_chembl_id"])
    # every parent has exactly one decision row
    parents = m[m["molecule_chembl_id"] == m["parent_chembl_id"]]
    assert parents["parent_chembl_id"].is_unique
    assert set(parents["parent_chembl_id"]) == set(m["parent_chembl_id"])
    # counts in the sidecar add up
    side = t2d.sidecar("2026-10-02T00:00:00Z")
    assert sum(side["counts"]["by_reason"].values()) == len(m)
    assert side["indication_rows"] == len(rows) == 376


def test_sorted_by_parent_chembl_id(t2d):
    nums = [int(p[6:]) for p in t2d.manifest["parent_chembl_id"]]
    assert nums == sorted(nums)


def test_included_unique_inchikey_small_molecules(t2d):
    inc = t2d.included
    assert inc["standard_inchi_key"].is_unique and (inc["standard_inchi_key"] != "").all()
    assert (inc["molecule_type"] == "Small molecule").all()
    assert (inc["molecule_chembl_id"] == inc["parent_chembl_id"]).all()
    assert (inc["phase_for_indication"] >= 4).all()
    assert all(structure_reason(s) is None for s in inc["canonical_smiles"])
    assert len(inc) == 33


def test_sitagliptin_salts_collapse_to_parent(t2d):
    m = t2d.manifest
    fam = m[m["parent_chembl_id"] == "CHEMBL1422"]
    assert _row(t2d, "CHEMBL1422")["decision"] == "include"
    salts = fam[fam["molecule_chembl_id"] != "CHEMBL1422"]
    assert len(salts) >= 3
    assert (salts["reason"] == "collapsed_to_parent").all()
    assert salts["pref_name"].str.startswith("SITAGLIPTIN").all()
    # the cohort holds sitagliptin once, as the neutral parent structure
    assert (t2d.included["pref_name"].str.startswith("SITAGLIPTIN")).sum() == 1
    assert "." not in _row(t2d, "CHEMBL1422")["canonical_smiles"]


def test_phase_is_max_over_salt_family(client, t2d):
    # bromocriptine (CHEMBL493) has no T2D row of its own; phase 4 comes from the mesylate
    assert not [r for r in client.drug_indications([T2D]) if r["molecule_chembl_id"] == "CHEMBL493"]
    row = _row(t2d, "CHEMBL493")
    assert row["phase_for_indication"] == 4.0 and row["decision"] == "include"
    # saxagliptin: the family holds a phase-3 member, the parent takes the max (4)
    fam = t2d.manifest[t2d.manifest["parent_chembl_id"] == "CHEMBL385517"]
    assert 3.0 in set(fam["phase_for_indication"])
    assert _row(t2d, "CHEMBL385517")["phase_for_indication"] == 4.0


def test_exclusion_reasons_t2d(t2d):
    rc = t2d.reason_counts()
    assert rc["not_small_molecule"] == 13
    assert _row(t2d, "CHEMBL1200711")["reason"] == "mixture"  # omega-3-acid ethyl esters
    assert _row(t2d, "CHEMBL2110563")["reason"] == "no_smiles"  # cyanocobalamin
    assert _row(t2d, "CHEMBL1200711")["phase_for_indication"] == 4.0


def test_phase_threshold(client, gout):
    allp = build_cohort([GOUT], client, min_phase=None)
    p3 = build_cohort([GOUT], client, min_phase=3)
    assert len(gout.included) < len(p3.included) < len(allp.included)
    assert "phase_below_min" not in allp.reason_counts()
    assert set(gout.included["parent_chembl_id"]) <= set(p3.included["parent_chembl_id"])


def test_dailymed_only_filter(client, t2d):
    nodm = build_cohort([T2D], client, exclude_sole_evidence={"DailyMed"})
    dropped = set(t2d.included["pref_name"]) - set(nodm.included["pref_name"])
    assert dropped == {
        "BIOTIN", "CHLORPROPAMIDE", "FENOFIBRIC ACID", "PYRIDOXINE", "TOLAZAMIDE", "TOLBUTAMIDE", "TYROSINE",
    }
    assert set(nodm.included["pref_name"]) <= set(t2d.included["pref_name"])
    assert _row(nodm, "CHEMBL782")["reason"] == "no_qualifying_evidence"  # tolbutamide
    assert math.isnan(_row(nodm, "CHEMBL782")["phase_for_indication"])
    # an allow-list that keeps every observed source equals no filter
    allowed = {"ATC", "ClinicalTrials", "DailyMed", "EMA", "FDA", "USAN"}
    same = build_cohort([T2D], client, evidence_sources=allowed)
    assert same.to_csv_bytes() == t2d.to_csv_bytes()
    assert nodm.sidecar("t")["parameters"]["exclude_sole_evidence"] == ["DailyMed"]


def test_mesh_ids_must_be_confirmed_descriptors(client):
    with pytest.raises(ValueError, match="not MeSH descriptor IDs"):
        build_cohort(["gout"], client)
    with pytest.raises(ValueError, match="no drug_indication rows"):
        build_cohort(["D000001"], client)
    with pytest.raises(ValueError, match="empty"):
        build_cohort([], client)


def test_include_descendants(client, mesh, gout):
    with pytest.raises(ValueError, match="MeshClient"):
        build_cohort([GOUT], client, include_descendants=True)
    wide = build_cohort([GOUT], client, include_descendants=True, mesh_client=mesh)
    assert wide.parameters["mesh_ids"] == [GOUT, ARTHRITIS_GOUTY]
    assert wide.extra["mesh_descendants_added"] == {ARTHRITIS_GOUTY: "Arthritis, Gouty"}
    assert (wide.manifest["mesh_ids"] == f"{GOUT};{ARTHRITIS_GOUTY}").all()
    assert set(gout.included["parent_chembl_id"]) <= set(wide.included["parent_chembl_id"])
    # order of the requested IDs does not matter
    assert build_cohort([GOUT], client).to_csv_bytes() == build_cohort([GOUT, GOUT], client).to_csv_bytes()


# ----------------------------------------------------------------------------- suggestions
def test_suggest_headings(client, mesh):
    s = suggest_headings("gout", client)
    assert s.iloc[0]["mesh_id"] == GOUT
    assert list(s.columns) == ["mesh_id", "mesh_heading", "n_rows", "score", "matched_via"]
    s2 = suggest_headings("GOUT", client)
    assert s2.equals(s)  # case-insensitive
    bc = suggest_headings("breast cancer", client, mesh_client=mesh)
    assert bc.iloc[0]["mesh_id"] == "D001943"
    assert bc.iloc[0]["matched_via"] == "MeSH term: Breast Cancer"
    plain = suggest_headings("breast cancer", client)
    assert plain.empty or plain.iloc[0]["score"] < bc.iloc[0]["score"]


def test_suggest_headings_type2(client):
    s = suggest_headings("type 2 diabetes", client)
    assert s.iloc[0]["mesh_id"] == T2D


# ----------------------------------------------------------------------------- naive / ATC
def test_naive_cohort_is_v0x_behaviour(client, gout):
    n = naive_cohort("gout", client)
    assert n.extra["mesh_headings_matched"] == {GOUT: "Gout", ARTHRITIS_GOUTY: "Arthritis, Gouty"}
    inc = n.included
    # all phases kept, salt forms kept as separate molecules, no small-molecule filter
    assert (inc["phase_for_indication"] < 4).any()
    assert (inc["molecule_chembl_id"] != inc["parent_chembl_id"]).sum() == 3
    assert "CHEMBL1034" in set(inc["molecule_chembl_id"])  # diclofenac sodium
    assert n.reason_counts() == {"included": 37, "no_smiles": 8}
    nums = [int(x[6:]) for x in n.manifest["molecule_chembl_id"]]
    assert nums == sorted(nums)
    # controlled is a strict, much smaller subset at the parent level
    assert set(gout.included["parent_chembl_id"]) <= set(inc["parent_chembl_id"])
    j = jaccard(inc["standard_inchi_key"], gout.included["standard_inchi_key"])
    assert 0 < j < 0.5
    assert n.sidecar("t")["counts"]["parents"] is None


def test_atc_cohort(client):
    a = build_atc_cohort("M04", client)
    assert a.disease_label == "ATC M04 (ANTIGOUT PREPARATIONS)"
    names = set(a.included["pref_name"])
    assert {"ALLOPURINOL", "FEBUXOSTAT", "COLCHICINE", "PROBENECID"} <= names
    assert (a.manifest["phase_for_indication"].isna()).all()
    assert a.included["indication_ref_types"].str.startswith("ATC:M04").all()
    assert a.reason_counts() == {"included": 10, "not_small_molecule": 2}
    assert a.to_csv_bytes() == build_atc_cohort(["m04"], client).to_csv_bytes()
    with pytest.raises(ValueError):
        build_atc_cohort("diabetes", client)


def test_jaccard():
    assert jaccard({"a", "b"}, {"b", "c"}) == pytest.approx(1 / 3)
    assert math.isnan(jaccard([], [""]))


# ----------------------------------------------------------------------------- determinism / IO
def test_write_csv_and_sidecar(tmp_path, gout):
    csv_path, json_path = gout.write(tmp_path / "gout.csv", timestamp_utc="2026-10-02T12:00:00Z")
    assert csv_path.read_bytes() == gout.to_csv_bytes()
    side = json.loads(json_path.read_text())
    assert side["chembl_release"] == RELEASE
    assert side["created_utc"] == "2026-10-02T12:00:00Z"
    assert side["parameters"]["min_phase"] == 4.0
    assert side["counts"]["included"] == 10
    assert side["counts"]["by_reason"] == gout.reason_counts()
    assert side["topoqspr_version"]
    import hashlib

    assert side["manifest_sha256"] == hashlib.sha256(csv_path.read_bytes()).hexdigest()
    back = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    assert list(back.columns) == list(MANIFEST_COLUMNS)
    assert len(back) == len(gout.manifest)


def test_same_inputs_same_bytes_in_process(client):
    a = build_cohort([T2D], client).to_csv_bytes()
    b = build_cohort([T2D], ChemblClient(FIX / "chembl", offline=True)).to_csv_bytes()
    assert a == b


_SCRIPT = """
import sys
from topoqspr.chembl import ChemblClient
from topoqspr.cohort import build_cohort, naive_cohort, build_atc_cohort
c = ChemblClient(sys.argv[1], offline=True)
out = build_cohort(['D003924'], c).to_csv_bytes()
out += build_cohort(['D006073'], c, exclude_sole_evidence={'DailyMed'}).to_csv_bytes()
out += naive_cohort('gout', c).to_csv_bytes()
out += build_atc_cohort('M04', c).to_csv_bytes()
sys.stdout.buffer.write(out)
"""


def test_deterministic_across_hash_seeds(tmp_path):
    outs = []
    for seed in ("0", "4242"):
        env = {**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": str(SRC)}
        res = subprocess.run(
            [sys.executable, "-c", _SCRIPT, str(FIX / "chembl")], env=env, capture_output=True, check=True
        )
        outs.append(res.stdout)
    assert outs[0] == outs[1]
    assert len(outs[0]) > 1000


# ----------------------------------------------------------------------------- live
@pytest.mark.network
def test_live_status(tmp_path):
    c = ChemblClient(tmp_path)
    assert c.release.startswith("ChEMBL_")


@pytest.mark.network
def test_live_gout_matches_fixture(tmp_path):
    c = ChemblClient(tmp_path)
    if c.release != RELEASE:
        pytest.skip(f"live release {c.release} differs from fixture release {RELEASE}")
    shutil.copytree(FIX / "chembl" / RELEASE / "drug_indication__mesh_headings",
                    tmp_path / RELEASE / "drug_indication__mesh_headings")
    assert build_cohort([GOUT], c).to_csv_bytes() == (FIX / "expected_gout.csv").read_bytes()


@pytest.mark.network
def test_live_mesh_descendants(tmp_path):
    m = MeshClient(tmp_path)
    assert {"mesh_id": ARTHRITIS_GOUTY, "mesh_heading": "Arthritis, Gouty"} in m.descendants(GOUT)
    assert "D001943" in m.descriptors_for_text("breast cancer")
