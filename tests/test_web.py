"""Tests for the v1.0 web interface (src/topoqspr/web; docs/WEB.md). Offline.

Real data: only the quarantined pilot cohort (Migraine Disorders, D008881) from the
recorded ChEMBL_37 cache in tests/fixtures/pipeline/chembl. The job lifecycle runs the
real worker (child process) with a cut-down evaluation (ridge, 3 outer folds, one
target) so the module stays fast; no performance number is asserted or printed.
"""

from __future__ import annotations

import io
import json
import re
import time
import uuid
import zipfile
from pathlib import Path

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")  # fastapi.testclient

from fastapi.testclient import TestClient  # noqa: E402

from topoqspr.web import jobs as jb  # noqa: E402
from topoqspr.web.app import Settings, create_app, parse_smiles_text  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
CHEMBL = REPO / "tests" / "fixtures" / "pipeline" / "chembl"
STATIC = REPO / "src" / "topoqspr" / "web" / "static"
MIGRAINE = "D008881"
FAST = dict(families=("ridge",), outer_splits=3, inner_splits=2, seeds=(0,), n_jobs=1)


def _settings(data_dir: Path, **over) -> Settings:
    base = dict(data_dir=data_dir, chembl_cache_dir=str(CHEMBL), offline=True,
                rate_api_per_min=10_000, rate_cohort_per_min=1000, rate_indices_per_min=1000,
                rate_runs_per_hour=1000, worker=False, **FAST)
    base.update(over)
    return Settings(**base)


def _cohort(client, mesh=MIGRAINE, **over):
    payload = {"mesh_ids": [mesh], "min_phase": 4}
    payload.update(over)
    return client.post("/api/cohort", json=payload)


def _wait(client, job_id, timeout=900):
    t0 = time.monotonic()
    while time.monotonic() - t0 < timeout:
        job = client.get(f"/api/jobs/{job_id}").json()
        if job["status"] in ("done", "failed"):
            return job
        time.sleep(1.0)
    raise AssertionError(f"job {job_id} did not finish in {timeout} s: {job}")


# ------------------------------------------------------------------ lifecycle
@pytest.fixture(scope="module")
def data_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("web")


@pytest.fixture(scope="module")
def done_job(data_dir):
    """Run one job end to end with the real worker; returns (job_id, cohort preview)."""
    app = create_app(_settings(data_dir, worker=True))
    with TestClient(app) as c:
        prev = _cohort(c).json()
        r = c.post("/api/run", json={"cohort_token": prev["cohort_token"], "confirmed": True,
                                     "mode": "exploratory", "targets": {"computed": ["TPSA"]}})
        assert r.status_code == 202, r.text
        job = r.json()
        assert job["status"] == "queued" and jb.valid_id(job["job_id"])
        final = _wait(c, job["job_id"])
        assert final["status"] == "done", final
        return job["job_id"], prev


def test_health_and_meta(tmp_path):
    with TestClient(create_app(_settings(tmp_path))) as c:
        h = c.get("/health")
        assert h.status_code == 200 and h.json()["status"] == "ok"
        m = c.get("/api/meta").json()
        assert m["mode"]["name"] == "exploratory"
        assert m["mode"]["feature_sets"] == ["null", "size", "composition", "TI", "rdkit2d"]
        assert m["mode"]["permutation"] is False and m["mode"]["bootstrap"] is False
        assert "command-line protocol" in m["mode"]["note"]
        assert {t["name"] for t in m["measured_targets"]} >= {"mp_bradley", "logp_opera"}
        assert all(t["licence"] for t in m["measured_targets"])
        assert "CC BY-SA 3.0" in m["chembl"]["licence"]


def test_health_reports_dead_worker(tmp_path):
    app = create_app(_settings(tmp_path, worker=True))
    with TestClient(app) as c:
        app.state.web["worker"].stop()
        r = c.get("/health")
        assert r.status_code == 503 and r.json()["worker_alive"] is False


def test_headings_from_release_index(tmp_path):
    with TestClient(create_app(_settings(tmp_path))) as c:
        r = c.get("/api/headings", params={"q": "migraine"})
        assert r.status_code == 200, r.text
        data = r.json()
        assert data["release"] == "ChEMBL_37"
        ids = [h["mesh_id"] for h in data["headings"]]
        assert MIGRAINE in ids
        assert all(h["allowed"] for h in data["headings"])
        assert c.get("/api/headings", params={"q": "m"}).status_code == 422
        assert c.get("/api/headings", params={"q": "x" * 101}).status_code == 422
        assert c.get("/api/headings", params={"q": "mi\x00graine"}).status_code == 422


def test_cohort_preview_matches_pilot_manifest(tmp_path):
    with TestClient(create_app(_settings(tmp_path))) as c:
        r = _cohort(c)
        assert r.status_code == 200, r.text
        p = r.json()
        assert jb.valid_id(p["cohort_token"])
        assert p["chembl_release"] == "ChEMBL_37" and p["mesh_ids"] == [MIGRAINE]
        # the pilot cohort of docs/PIPELINE.md §8: 152 manifest rows, 33 included parents
        assert p["counts"]["manifest_rows"] == 152 and p["n_included"] == 33
        assert len(p["included"]) == 33 and len(p["excluded"]) == 152 - 33
        assert p["counts"]["by_reason"]["phase_below_min"] == 84
        assert all(e["reason"] for e in p["excluded"])
        assert p["runnable"] is True and "CC BY-SA 3.0" in p["data_licence"]
        man = c.get(f"/api/cohorts/{p['cohort_token']}/manifest.csv")
        assert man.status_code == 200 and man.headers["content-type"].startswith("text/csv")
        assert man.text.splitlines()[0].startswith("disease_label,mesh_ids,chembl_release")
        assert len(man.text.splitlines()) == 153
        # same inputs + release -> same manifest hash
        assert _cohort(c).json()["manifest_sha256"] == p["manifest_sha256"]


def test_cohort_validation_errors(tmp_path):
    with TestClient(create_app(_settings(tmp_path))) as c:
        assert _cohort(c, mesh="migraine").status_code == 422
        assert c.post("/api/cohort", json={"mesh_ids": []}).status_code == 422
        assert c.post("/api/cohort", json={"mesh_ids": [MIGRAINE], "bogus": 1}).status_code == 422
        assert _cohort(c, min_phase=2.5).status_code == 422
        assert _cohort(c, evidence_sources=["fda"]).status_code == 422
        r = _cohort(c, mesh="D999999")  # well-formed, but not in the release's index
        assert r.status_code == 400 and "D999999" in r.json()["detail"]
        assert c.post("/api/cohort", content=b"not json",
                      headers={"content-type": "application/json"}).status_code == 422


def test_job_lifecycle_and_results(done_job, data_dir):
    job_id, prev = done_job
    app = create_app(_settings(data_dir))
    with TestClient(app) as c:
        job = c.get(f"/api/jobs/{job_id}").json()
        assert job["status"] == "done" and job["progress"] == 1.0
        assert job["cohort"]["mesh_ids"] == [MIGRAINE] and job["cohort"]["n_included"] == 33
        res = c.get(f"/api/jobs/{job_id}/results")
        assert res.status_code == 200
        data = res.json()
        assert [t["name"] for t in data["targets"]] == ["TPSA"]
        t = data["targets"][0]
        assert [r["set"] for r in t["rows"]] == ["null", "size", "composition", "TI", "rdkit2d"]
        assert t["reference"] == "size"
        size = next(r for r in t["rows"] if r["set"] == "size")
        assert size["is_reference"] and size["d_r2_mean"] is None
        # n = 33 < min_n = 50: every row is labelled underpowered, as the protocol requires
        assert all(r["underpowered"] is True for r in t["rows"])
        assert all(r["n"] == 33 for r in t["rows"])
        ti = next(r for r in t["rows"] if r["set"] == "TI")
        assert set(ti) >= {"d_r2_ci_low", "d_r2_ci_high", "d_mae_mean", "ad_share_outside"}
        assert data["protocol"]["mode"] == "exploratory"
        assert data["chembl"]["release"] == "ChEMBL_37"
        assert data["chembl"]["manifest_verified"] is True
        assert data["manifest_sha256"] == prev["manifest_sha256"]
        # downloads
        z = c.get(f"/api/jobs/{job_id}/results.zip")
        assert z.status_code == 200 and z.headers["content-type"] == "application/zip"
        names = zipfile.ZipFile(io.BytesIO(z.content)).namelist()
        for n in ("run_record.json", "LICENCE-DATA.txt", "README.txt", "manifest.csv",
                  "summary.csv", "ad_summary.csv"):
            assert f"topoqspr-results/{n}" in names
        lic = zipfile.ZipFile(io.BytesIO(z.content)).read("topoqspr-results/LICENCE-DATA.txt")
        assert b"CC BY-SA 3.0" in lic and b"ChEMBL_37" in lic
        rec = c.get(f"/api/jobs/{job_id}/run_record.json").json()
        assert rec["config"]["evaluation"]["learners"] == ["selected"]
        assert rec["config"]["cohort"]["source"] == "manifest"
        man = c.get(f"/api/jobs/{job_id}/manifest.csv")
        assert man.status_code == 200 and len(man.text.splitlines()) == 153


def test_results_zip_is_byte_stable(done_job, data_dir):
    job_id, _ = done_job
    p = data_dir / "jobs" / job_id / "results.zip"
    with zipfile.ZipFile(p) as zf:
        assert {i.date_time for i in zf.infolist()} == {(2000, 1, 1, 0, 0, 0)}


def test_persistence_across_restart(done_job, data_dir):
    """A finished job survives a restart; a job left running is failed with the reason;
    a queued job stays queued and is picked up."""
    job_id, prev = done_job
    store = jb.JobStore(data_dir / "web.sqlite")
    running, queued = jb.new_id(), jb.new_id()
    req = {"created_utc": jb.now_utc(), "mode": "exploratory", "targets": {"computed": ["TPSA"]},
           "cohort": {}}
    store.create_job(running, prev["cohort_token"], req, "1.2.3.4", jb.now_utc())
    store.create_job(queued, prev["cohort_token"], req, "1.2.3.5", jb.now_utc())
    assert store.claim_next()["id"] == running  # simulate a crash while running
    assert store.claim_next() is None           # one running job at a time
    with TestClient(create_app(_settings(data_dir))) as c:
        r = c.get(f"/api/jobs/{running}").json()
        assert r["status"] == "failed" and "restarted" in r["error"]
        q = c.get(f"/api/jobs/{queued}").json()
        assert q["status"] == "queued" and q["queue_position"] == 1
        assert c.get(f"/api/jobs/{job_id}").json()["status"] == "done"
        assert c.get(f"/api/jobs/{job_id}/results").status_code == 200
    store.finish(queued, "failed", error="cleanup", message="failed")


def test_job_timeout_fails_with_reason(tmp_path):
    app = create_app(_settings(tmp_path, worker=True, job_timeout_s=0.5))
    with TestClient(app) as c:
        tok = _cohort(c).json()["cohort_token"]
        r = c.post("/api/run", json={"cohort_token": tok, "confirmed": True,
                                     "targets": {"computed": ["TPSA"]}})
        assert r.status_code == 202
        job = _wait(c, r.json()["job_id"], timeout=120)
        assert job["status"] == "failed" and "timed out" in job["error"]
        assert c.get(f"/api/jobs/{job['job_id']}/results").status_code == 409


# ------------------------------------------------------------------ caps and validation
def test_run_validation(tmp_path):
    with TestClient(create_app(_settings(tmp_path))) as c:
        tok = _cohort(c).json()["cohort_token"]
        ok = {"cohort_token": tok, "confirmed": True, "targets": {"computed": ["TPSA"]}}
        assert c.post("/api/run", json={**ok, "confirmed": False}).status_code == 400
        assert c.post("/api/run", json={**ok, "targets": {}}).status_code == 422
        assert c.post("/api/run", json={**ok, "targets": {"computed": ["MolWt2"]}}).status_code == 422
        assert c.post("/api/run", json={**ok, "targets": {"measured": ["nope"]}}).status_code == 422
        assert c.post("/api/run", json={**ok, "mode": "confirmatory"}).status_code == 422
        assert c.post("/api/run", json={**ok, "cohort_token": str(uuid.uuid4())}).status_code == 404
        assert c.post("/api/run", json={**ok, "cohort_token": "x" * 36}).status_code == 422
        r = c.post("/api/run", json={**ok, "targets": {"computed": ["TPSA"],
                                                       "measured": ["logp_opera"]}})
        assert r.status_code == 202
        cfg = json.loads((tmp_path / "jobs" / r.json()["job_id"] / "config.json").read_text())
        assert [t["provider"] for t in cfg["targets"]] == ["rdkit", "expdata"]
        assert cfg["evaluation"]["learners"] == ["selected"]
        assert cfg["feature_sets"] == list(jb.LADDER_SETS)


def test_web_mode_defaults():
    s = Settings()
    assert s.families == ("ridge", "pls", "knn") and s.outer_splits == 5 and s.seeds == (0,)
    assert s.max_molecules == 300 and s.max_queue >= 1 and 1 <= s.n_jobs <= 4


def test_molecule_cap(tmp_path):
    with TestClient(create_app(_settings(tmp_path, max_molecules=20))) as c:
        p = _cohort(c).json()
        assert p["runnable"] is False and "limit is 20" in p["not_runnable_reason"]
        r = c.post("/api/run", json={"cohort_token": p["cohort_token"], "confirmed": True,
                                     "targets": {"computed": ["TPSA"]}})
        assert r.status_code == 400 and "limit" in r.json()["detail"]


def test_queue_cap_and_one_active_job_per_ip(tmp_path):
    with TestClient(create_app(_settings(tmp_path, max_queue=1))) as c:
        tok = _cohort(c).json()["cohort_token"]
        body = {"cohort_token": tok, "confirmed": True, "targets": {"computed": ["TPSA"]}}
        assert c.post("/api/run", json=body).status_code == 202
        assert c.post("/api/run", json=body).status_code == 503
    with TestClient(create_app(_settings(tmp_path / "b", max_queue=5))) as c:
        tok = _cohort(c).json()["cohort_token"]
        body = {"cohort_token": tok, "confirmed": True, "targets": {"computed": ["TPSA"]}}
        assert c.post("/api/run", json=body).status_code == 202
        r = c.post("/api/run", json=body)
        assert r.status_code == 429 and "already" in r.json()["detail"]


def test_rate_limits(tmp_path):
    with TestClient(create_app(_settings(tmp_path, rate_cohort_per_min=2))) as c:
        assert _cohort(c).status_code == 200
        assert _cohort(c).status_code == 200
        r = _cohort(c)
        assert r.status_code == 429 and int(r.headers["retry-after"]) >= 1
    with TestClient(create_app(_settings(tmp_path / "b", rate_api_per_min=3))) as c:
        codes = [c.get("/api/meta").status_code for _ in range(4)]
        assert codes == [200, 200, 200, 429]
        assert c.get("/health").status_code == 200  # not rate limited
    with TestClient(create_app(_settings(tmp_path / "c", rate_runs_per_hour=1,
                                         max_active_per_ip=5))) as c:
        tok = _cohort(c).json()["cohort_token"]
        body = {"cohort_token": tok, "confirmed": True, "targets": {"computed": ["TPSA"]}}
        assert c.post("/api/run", json=body).status_code == 202
        assert c.post("/api/run", json=body).status_code == 429


def test_proxy_client_ip(tmp_path):
    with TestClient(create_app(_settings(tmp_path, rate_api_per_min=1, trust_proxy=True))) as c:
        assert c.get("/api/meta", headers={"x-forwarded-for": "9.9.9.9, 1.1.1.1"}).status_code == 200
        assert c.get("/api/meta", headers={"x-forwarded-for": "9.9.9.9, 2.2.2.2"}).status_code == 200
        assert c.get("/api/meta", headers={"x-forwarded-for": "8.8.8.8, 2.2.2.2"}).status_code == 429


def test_request_size_limit(tmp_path):
    with TestClient(create_app(_settings(tmp_path, max_body_bytes=1000))) as c:
        r = c.post("/api/indices", json={"text": "C" * 5000})
        assert r.status_code == 413
        assert c.post("/api/indices", json={"text": "CCO"}).status_code == 200


def test_mesh_allowlist(tmp_path):
    with TestClient(create_app(_settings(tmp_path, mesh_allowlist=("D000001",)))) as c:
        r = _cohort(c)
        assert r.status_code == 403 and "not enabled" in r.json()["detail"]
        hs = c.get("/api/headings", params={"q": "migraine"}).json()["headings"]
        assert hs and not any(h["allowed"] for h in hs)
        assert c.get("/api/meta").json()["mesh_allowlist"] == ["D000001"]
    with TestClient(create_app(_settings(tmp_path / "b", mesh_allowlist=(MIGRAINE,)))) as c:
        assert _cohort(c).status_code == 200


# ------------------------------------------------------------------ indices
def test_indices_zip(tmp_path):
    with TestClient(create_app(_settings(tmp_path))) as c:
        r = c.post("/api/indices", json={"text": "smiles,id\nCCO,ethanol\nc1ccccc1,benzene\n"
                                                 "[Na+].[Cl-],salt\nnot_a_smiles,bad\n"})
        assert r.status_code == 200, r.text
        assert r.headers["x-molecules"] == "4"
        zf = zipfile.ZipFile(io.BytesIO(r.content))
        rows = zf.read("topoqspr-indices/indices.csv").decode().splitlines()
        assert rows[0].startswith("id,input_smiles,status,reason,smiles,inchikey,")
        assert len(rows) == 5
        assert rows[1].startswith("ethanol,CCO,ok")
        fails = zf.read("topoqspr-indices/indices.failures.csv").decode()
        assert "bad,*,rejected:" in fails
        assert int(r.headers["x-molecules-ok"]) == 2


def test_indices_cap_and_errors(tmp_path):
    with TestClient(create_app(_settings(tmp_path, max_indices_molecules=3))) as c:
        r = c.post("/api/indices", json={"text": "C\nCC\nCCC\nCCCC\n"})
        assert r.status_code == 422 and "limit is 3" in r.json()["detail"]
        assert c.post("/api/indices", json={"text": ""}).status_code == 422
        assert c.post("/api/indices", json={"text": "CC a\nCCC a"}).status_code == 422


def test_parse_smiles_text():
    assert parse_smiles_text("CCO ethanol\n\n# c\nc1ccccc1", 10, 100) == [
        ("ethanol", "CCO"), ("mol2", "c1ccccc1")]
    assert parse_smiles_text("id,smiles\na,CCO\nb,CC", 10, 100) == [("a", "CCO"), ("b", "CC")]
    assert parse_smiles_text("SMILES\nCCO\n", 10, 100) == [("mol1", "CCO")]
    with pytest.raises(ValueError, match="longer than"):
        parse_smiles_text("C" * 101, 10, 100)
    with pytest.raises(ValueError, match="duplicate"):
        parse_smiles_text("C x\nCC x", 10, 100)


# ------------------------------------------------------------------ path safety
@pytest.mark.parametrize("bad", [
    "..", "..%2F..%2Fetc%2Fpasswd", "not-a-uuid", str(uuid.uuid1()),
    str(uuid.uuid4()).upper(), "00000000-0000-4000-8000-000000000000x",
])
def test_download_path_safety(tmp_path, bad):
    with TestClient(create_app(_settings(tmp_path))) as c:
        for suffix in ("", "/results", "/results.zip", "/run_record.json", "/manifest.csv"):
            r = c.get(f"/api/jobs/{bad}{suffix}")
            assert r.status_code == 404, (bad, suffix, r.status_code)
        assert c.get(f"/api/cohorts/{bad}/manifest.csv").status_code == 404
        assert c.get("/static/../web.sqlite").status_code == 404
        assert c.get("/static/..%2F..%2Fapp.py").status_code == 404


def test_unknown_but_valid_job_and_not_done(tmp_path):
    with TestClient(create_app(_settings(tmp_path))) as c:
        assert c.get(f"/api/jobs/{uuid.uuid4()}").status_code == 404
        tok = _cohort(c).json()["cohort_token"]
        jid = c.post("/api/run", json={"cohort_token": tok, "confirmed": True,
                                       "targets": {"computed": ["TPSA"]}}).json()["job_id"]
        for suffix in ("/results", "/results.zip", "/run_record.json", "/manifest.csv"):
            assert c.get(f"/api/jobs/{jid}{suffix}").status_code == 409


def test_valid_id():
    assert jb.valid_id(str(uuid.uuid4()))
    for bad in ("", "../x", str(uuid.uuid1()), str(uuid.uuid4()).upper(), None, 5):
        assert not jb.valid_id(bad)


# ------------------------------------------------------------------ front-end
def test_frontend_never_injects_html():
    js = (STATIC / "app.js").read_text(encoding="utf-8")
    for bad in ("innerHTML", "outerHTML", "insertAdjacentHTML", "document.write", "eval(",
                "new Function", "setTimeout(\"", "setTimeout('"):
        assert bad not in js, bad
    # no network call except this server's API: relative /api/ paths only
    assert not re.search(r"https?://", js)
    assert all(p.startswith("/api/") for p in re.findall(r'"(/api/[^"]*)"', js))
    html = (STATIC / "index.html").read_text(encoding="utf-8")
    assert not re.search(r"<script(?![^>]*\bsrc=\"/static/)", html)  # no inline scripts
    assert not re.search(r"\son[a-z]+\s*=", html)                    # no inline handlers
    assert not re.search(r"(src|href)=\"https?://", html)            # no external assets


def test_frontend_accessibility_basics():
    html = (STATIC / "index.html").read_text(encoding="utf-8")
    assert '<html lang="en">' in html and 'name="viewport"' in html
    ids = re.findall(r'<(?:input|select|textarea)[^>]*\bid="([^"]+)"', html)
    labelled = set(re.findall(r'<label[^>]*\bfor="([^"]+)"', html))
    assert ids and set(ids) <= labelled, set(ids) - labelled
    assert 'aria-live="polite"' in html


def test_index_page_and_security_headers(tmp_path):
    with TestClient(create_app(_settings(tmp_path))) as c:
        r = c.get("/")
        assert r.status_code == 200 and "TopoQSPR" in r.text
        csp = r.headers["content-security-policy"]
        assert "script-src 'self'" in csp and "connect-src 'self'" in csp
        assert r.headers["x-content-type-options"] == "nosniff"
        assert c.get("/static/app.js").status_code == 200
        assert c.get("/api/meta").headers["cache-control"] == "no-store"


def test_cors_restricted_to_configured_origin(tmp_path):
    origin = "https://topoqspr.example.org"
    with TestClient(create_app(_settings(tmp_path, allowed_origins=(origin,)))) as c:
        ok = c.get("/api/meta", headers={"origin": origin})
        assert ok.headers.get("access-control-allow-origin") == origin
        bad = c.get("/api/meta", headers={"origin": "https://evil.example"})
        assert "access-control-allow-origin" not in bad.headers
    with TestClient(create_app(_settings(tmp_path / "b"))) as c:  # default: same origin only
        r = c.get("/api/meta", headers={"origin": origin})
        assert "access-control-allow-origin" not in r.headers


def test_text_cells_never_leak_nan():
    """Regression: an empty 'reason' cell read by pandas is NaN; it must become None, not NaN."""
    import json as _json
    from topoqspr.web.jobs import _text
    assert _text(float("nan")) is None and _text("") is None and _text(None) is None
    assert _text(" small n ") == "small n"
    _json.dumps({"r": _text(float("nan"))}, allow_nan=False)


def test_indices_formats_entropies_view_and_catalog(data_dir):
    """ChemGraphX parity and beyond: file/graph formats, entropies, on-page view, formula catalogue."""
    from rdkit import Chem
    app = create_app(_settings(data_dir / "fmt"))
    with TestClient(app) as c:
        r = c.post("/api/indices", json={"text": "A = [[0,1,0,1],[1,0,1,0],[0,1,0,1],[1,0,1,0]]",
                                         "format": "adjacency_matrix", "view": True})
        assert r.status_code == 200, r.text
        v = r.json()
        row = dict(zip(v["columns"], v["rows"][0]))
        assert row["kind"] == "graph" and row["M1"] == 16
        assert abs(row["ENT_M1"] - 1.3862943611198906) < 1e-12  # ln 4 on the 4-cycle
        mb = Chem.MolToMolBlock(Chem.MolFromSmiles("c1ccccc1O"))
        r = c.post("/api/indices", json={"text": mb, "format": "mol", "view": True})
        assert dict(zip(r.json()["columns"], r.json()["rows"][0]))["M1"] == 30
        bad = c.post("/api/indices", json={"text": "0 1\n1 1", "format": "adjacency_matrix", "view": True}).json()
        assert bad["n_ok"] == 0 and "self-loops" in bad["failures"][0]["reason"]
        z = c.post("/api/indices", json={"text": "CCO ethanol", "format": "smiles"})
        assert z.status_code == 200 and z.headers["content-type"] == "application/zip"
        cat = c.get("/api/catalog").json()["items"]
        assert len(cat) == 104 + 91 and all(i["formula_latex"] and i["reference"] for i in cat)


def test_depiction_and_edge_partition(data_dir):
    app = create_app(_settings(data_dir / "dep"))
    with TestClient(app) as c:
        r = c.get("/api/depict", params={"smiles": "CC(=O)Oc1ccccc1C(=O)O.[Na+]"})
        assert r.status_code == 200 and r.headers["content-type"].startswith("image/svg+xml")
        assert "<svg" in r.text and "<script" not in r.text.lower()
        assert "default-src 'none'" in r.headers["content-security-policy"]
        assert c.get("/api/depict", params={"smiles": "[Na+].[Cl-]"}).status_code == 422
        v = c.post("/api/indices", json={"text": "CCO ethanol", "view": True}).json()
        row = dict(zip(v["columns"], v["rows"][0]))
        assert row["edge_partition"] == "(1,2)x2"


def test_chemgraphx_like_shell_and_time_taken(data_dir):
    """Layout follows ChemGraphX: Home/Instructions/Contact tabs, input-format buttons, a cite card,
    and a time-taken value with every indices result."""
    html = (STATIC / "index.html").read_text(encoding="utf-8")
    assert all(f'data-tab="{t}"' in html for t in ("home", "instructions", "contact"))
    assert html.count('class="fmt"') == 6 and "How to Cite" in html
    assert "<title>TopoQSPR</title>" in html
    app = create_app(_settings(data_dir / "shell"))
    with TestClient(app) as c:
        assert c.get("/static/logo.svg").status_code == 200
        v = c.post("/api/indices", json={"text": "CCO", "view": True}).json()
        assert isinstance(v["time_taken_s"], float) and v["time_taken_s"] > 0
        z = c.post("/api/indices", json={"text": "CCO"})
        assert float(z.headers["x-time-taken"]) > 0
