"""Fast, source-grounded adversarial audit checks (2026-10-12).

These are *regression and provenance diagnostics*, NOT a rerun of the original
nested-CV study. They never overwrite the frozen results. See
docs/INDEPENDENT_ADVERSARIAL_AUDIT_2026-10-12.md.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from topoqspr.pipeline import permutation_indices
from topoqspr.stats import corrected_ttest

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "results" / "expB" / "main_confirmatory"
SVR = ROOT / "results" / "expB" / "main_confirmatory_protocol_svr"
TARGETS = {"mp_bradley", "logp_opera", "logd74_az", "logs_aqsoldb"}


def _read_json(folder: Path, name: str):
    with (folder / name).open(encoding="utf-8") as f:
        return json.load(f)


def _read_csv(folder: Path, name: str):
    with (folder / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _bh(pvalues):
    """Small independent Benjamini-Hochberg reference (no project stats call)."""
    indexed = sorted(enumerate(pvalues), key=lambda item: item[1])
    out = [None] * len(pvalues)
    running = 1.0
    for j in range(len(indexed) - 1, -1, -1):
        index, p = indexed[j]
        running = min(running, p * len(indexed) / (j + 1))
        out[index] = running
    return out


def test_a1_record_controls_are_not_the_selected_model():
    """Do not accidentally claim the cheap controls retest all six families."""
    main = _read_json(MAIN, "run_record.json")
    amended = _read_json(SVR, "run_record.json")
    cfg = main["config"]
    cfg2 = amended["config"]

    assert cfg["protocol"]["learner"] == "selected"
    assert len(cfg["evaluation"]["families"]) == 6
    assert cfg["protocol"]["permutation"]["learner"] == "ridge"
    assert cfg["protocol"]["bootstrap"]["families"] == ["ridge"]
    assert cfg2["protocol"]["permutation"]["learner"] == "svr"
    assert cfg2["protocol"]["bootstrap"]["families"] == ["svr"]


def test_a5_registry_and_dirty_provenance_are_visible():
    """Mixed feature registries are already documented; do not conflate runs."""
    original = _read_json(MAIN, "run_record.json")
    amended = _read_json(SVR, "run_record.json")
    assert original["counts"]["feature_set_sizes"]["TI"] == 102
    assert amended["counts"]["feature_set_sizes"]["TI"] == 104
    assert original["code"]["dirty"] is False
    assert amended["code"]["dirty"] is True
    assert original["manifest_sha256"] == amended["manifest_sha256"]


def test_bh_matches_independent_reference_and_only_logp_is_positive():
    """Recompute four-member adjusted q values, not merely trust the CSV."""
    rows = _read_csv(MAIN, "confirmatory.csv")
    assert len(rows) == 4
    assert {r["target"] for r in rows} == TARGETS
    assert all(int(r["family_m"]) == 4 for r in rows)

    expected = _bh([float(r["d_r2_p"]) for r in rows])
    for r, q in zip(rows, expected):
        assert float(r["d_r2_q_bh"]) == pytest.approx(q, abs=1e-12)
    informative = {r["target"] for r in rows
                   if float(r["d_r2_q_bh"]) <= 0.05
                   and float(r["d_r2_mean"]) > 0}
    assert informative == {"logp_opera"}


def test_a1_original_bootstrap_instability_is_archived():
    """Archived pathological outcomes must not be discarded or called absent."""
    original = _read_csv(MAIN, "bootstrap.csv")
    amended = _read_csv(SVR, "bootstrap.csv")
    orig_logp = [float(r["d_r2"]) for r in original
                 if r["target"] == "logp_opera" and int(r["resample"]) >= 0]
    svr_logp = [float(r["d_r2"]) for r in amended
                if r["target"] == "logp_opera" and int(r["resample"]) >= 0]
    assert len(orig_logp) == 200
    assert len(svr_logp) == 200
    assert min(orig_logp) < -1000  # observed extreme extrapolation in Ridge
    assert min(svr_logp) > -1000   # narrower SVR control, NOT full selection
    assert all(np.isfinite(svr_logp))


def test_a3_permutation_is_a_bijection_but_not_a_null_proof():
    """A permutation of row indices alone does not prove exchangeability."""
    groups = np.array(["u", "u", "u", "s1", "s2", "s3"], dtype=object)
    changed_membership = False
    for seed in range(20):
        perm = permutation_indices(len(groups), groups,
                                   np.random.default_rng(seed))
        np.testing.assert_array_equal(np.sort(perm), np.arange(len(groups)))
        if set(perm[:3]) != {0, 1, 2}:
            changed_membership = True
    # Current implementation mixes the unique 3-member cluster and singletons.
    # This observation calls for a scientific null/exchangeability analysis;
    # it is *not* evidence that the Type I error is controlled.
    assert changed_membership


def test_a4_degenerate_corrected_ttest_behavior_is_explicit():
    """A statistical reporting convention, not a validated population claim."""
    out = corrected_ttest([0.1, 0.1, 0.1], n_train=8, n_test=2)
    assert out["p"] == 0.0
    assert out["ci_low"] == pytest.approx(0.1)
    assert out["ci_high"] == pytest.approx(0.1)
