"""Registered step adapters. Only these can run; a task cannot name arbitrary code or shell.

Each adapter declares its scope and which params are file paths (so the executor can allowlist
and hash them before anything runs). Adapters return a dict with `status` (ResultStatus value),
`findings` (list of str) and optional `evidence` (list of Evidence kwargs) and `data`.
"""
from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from researchos import evidence as ev
from researchos import inspectors
from researchos.safety import resolve_allowed
from researchos.schema import Independence, ResultStatus, Scope


@dataclass
class Adapter:
    name: str
    scope: Scope
    fn: Callable
    path_params: list[str] = field(default_factory=list)
    role: str = "verifier"
    doc: str = ""
    requires: list[str] = field(default_factory=list)   # tools that must be on PATH


REGISTRY: dict[str, Adapter] = {}


def adapter(name, scope, path_params=(), role="verifier", requires=()):
    def deco(fn):
        REGISTRY[name] = Adapter(name, Scope(scope), fn, list(path_params), role, (fn.__doc__ or "").strip(),
                                 list(requires))
        return fn
    return deco


def _read_csv(cfg, rel) -> list[dict]:
    p = resolve_allowed(rel, cfg.allow_roots, cfg.repo, cfg.private_paths)
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def _f(x) -> float:
    return float(x) if x not in ("", None) else math.nan


# ---------------------------------------------------------------------------------------------
@adapter("ingest_artifacts", "read-only", role="steward")
def ingest_artifacts(ctx, params, upstream):
    """Hash listed artifacts, check them against a SHA256SUMS file and their Git commit."""
    cfg = ctx.cfg
    d = Path(params["dir"])
    sums = {}
    if params.get("sums"):
        sums = ev.parse_sha256sums(resolve_allowed(d / params["sums"], cfg.allow_roots, cfg.repo, cfg.private_paths))
    rows, findings = [], []
    for name in params["files"]:
        rec = ev.check_artifact(cfg, str(d / name), sums.get(name) if sums else None, params.get("expect_commit"))
        rows.append(rec)
        if not rec["ok"]:
            findings.append(f"{rec['path']}: {rec.get('error')}")
    for name in sums:
        if name not in params["files"]:
            findings.append(f"{name}: listed in {params['sums']} but not ingested (kept outside Git: not verifiable here)")
    ok = all(r["ok"] for r in rows)
    return {"status": (ResultStatus.PASS if ok else ResultStatus.FAIL).value, "findings": findings,
            "data": {"artifacts": rows}}


@adapter("verify_external_points", "read-only", path_params=["summary", "paired", "recomputed"])
def verify_external_points(ctx, params, upstream):
    """Re-derive each archived contrast from the per-set summary with a separate CSV parser, and
    compare with the earlier independent recomputation from predictions."""
    S = _read_csv(ctx.cfg, params["summary"])
    P = _read_csv(ctx.cfg, params["paired"])
    R = _read_csv(ctx.cfg, params["recomputed"]) if params.get("recomputed") else []
    r2 = {(r["test"], r["set"]): _f(r["r2"]) for r in S}
    rec = {(r["test"], r["set"], r["reference"]): _f(r["d_r2"]) for r in R}
    rows, worst, worst_rec = [], 0.0, 0.0
    for r in P:
        k = (r["test"], r["set"]); kr = (r["test"], r["reference"])
        derived = r2[k] - r2[kr]
        diff = abs(derived - _f(r["d_r2"]))
        worst = max(worst, diff)
        rd = rec.get((r["test"], r["set"], r["reference"]))
        if rd is not None:
            worst_rec = max(worst_rec, abs(rd - _f(r["d_r2"])))
        rows.append({"test": r["test"], "set": r["set"], "reference": r["reference"], "d_r2_archived": _f(r["d_r2"]),
                     "d_r2_from_summary": derived, "abs_diff": diff,
                     "d_r2_independent_recompute": rd})
    tol = float(params.get("tol", 1e-12))
    ok = worst <= tol and (not R or worst_rec <= tol)
    findings = [f"{len(rows)} archived contrasts re-derived from per-set R²; max |diff| = {worst:.2e}",
                f"vs independent recomputation from predictions (earlier, separate script): max |diff| = {worst_rec:.2e}"
                if R else "no independent recomputation file supplied"]
    evid = [{
        "kind": "consistency_check", "independence": Independence.CONSISTENCY_CHECK,
        "summary": f"paired.csv contrasts equal per-set R² differences in summary.csv (max |diff| {worst:.1e})",
        "supports": "neutral", "reviewer": "researchos.adapters.verify_external_points",
        "method": "csv-module parse of two archived files; arithmetic only", "metrics": {"max_abs_diff": worst}}]
    if R:
        evid.append({
            "kind": "recomputation", "independence": Independence.INDEPENDENT_RECOMPUTATION,
            "summary": "point estimates match an earlier recomputation from predictions by a separate script "
                       "(docs/research/scripts/e3_recompute_and_failures.py)",
            "supports": "neutral", "reviewer": "docs/research/scripts/e3_recompute_and_failures.py",
            "method": "independent code path over predictions.csv.gz (held outside Git)",
            "metrics": {"max_abs_diff": worst_rec}})
    return {"status": (ResultStatus.PASS if ok else ResultStatus.FAIL).value, "findings": findings,
            "data": {"rows": rows}, "evidence": evid}


@adapter("learner_family_audit", "read-only", path_params=["external_summary", "internal_summaries"])
def learner_family_audit(ctx, params, upstream):
    """Compare two feature sets cell by cell: R² and the learner family each one selected."""
    a, b = params["set_a"], params["set_b"]
    comps = []
    ext = {(r["test"], r["set"]): r for r in _read_csv(ctx.cfg, params["external_summary"])}
    for test in sorted({k[0] for k in ext}):
        ra, rb = ext.get((test, a)), ext.get((test, b))
        if ra and rb:
            comps.append({"context": f"external/{test}", "r2_a": _f(ra["r2"]), "r2_b": _f(rb["r2"]),
                          "family_a": {ra["family"]: 1}, "family_b": {rb["family"]: 1}, "n": int(ra["n"]),
                          "uncertainty_on_difference": None, "fit": "single fit"})
    for rel in params.get("internal_summaries", []):
        rows = {r["set"]: r for r in _read_csv(ctx.cfg, rel)}
        ra, rb = rows.get(a), rows.get(b)
        if ra and rb:
            comps.append({"context": f"internal/{ra['endpoint']}_{ra['estimand']}", "r2_a": _f(ra["r2_fold_mean"]),
                          "r2_b": _f(rb["r2_fold_mean"]), "family_a": json.loads(ra["family_counts"]),
                          "family_b": json.loads(rb["family_counts"]), "n": int(ra["n"]),
                          "uncertainty_on_difference": None, "fit": "repeated nested CV (fold-mean)"})
    findings = []
    for c in comps:
        c["gap_a_minus_b"] = c["r2_a"] - c["r2_b"]
        c["same_learner"] = c["family_a"] == c["family_b"] and len(c["family_a"]) == 1
        if not c["same_learner"]:
            how = ("both mixed learner families across folds (selected per fold)" if c["family_a"] == c["family_b"]
                   else f"{a} used {c['family_a']}, {b} used {c['family_b']}")
            findings.append(f"{c['context']}: {how} -> R² gap {c['gap_a_minus_b']:+.3f} "
                            "is not a fixed-learner representation effect")
    findings.append(f"no archived interval exists for the {a} − {b} contrast in any cell")
    status = ResultStatus.FLAG if any(not c["same_learner"] for c in comps) else ResultStatus.PASS
    return {"status": status.value, "findings": findings, "data": {"comparisons": comps, "set_a": a, "set_b": b}}


@adapter("classify_claim", "read-only", role="coordinator")
def classify_claim(ctx, params, upstream):
    """Classify a claim from upstream comparisons. Automation may downgrade a claim; it never
    upgrades one to `supported` (that needs the owner)."""
    comps = upstream[params["from_step"]]["data"]["comparisons"]
    margin = float(params.get("margin", 0.01))
    against, confounded, favour = [], [], []
    for c in comps:
        # the claim says set_b (edge types) carries all of set_a's (all indices) gain: a > b + margin contradicts it
        gap = c["gap_a_minus_b"]
        (against if gap > margin else favour).append(c["context"])
        if not c["same_learner"]:
            confounded.append(c["context"])
    if against and confounded:
        verdict, why = "inconclusive", ("contradicted descriptively in " + ", ".join(against) +
                                        "; every contrast is learner-confounded and none has an interval")
    elif against:
        verdict, why = "refuted", "contradicted under matched learners in " + ", ".join(against)
    else:
        verdict, why = "inconclusive", "consistent descriptively but no interval on the difference"
    return {"status": ResultStatus.FLAG.value if verdict != "supported" else ResultStatus.PASS.value,
            "findings": [f"claim {params['claim_id']}: {verdict} ({why})",
                         f"in favour: {favour}; against (gap > {margin}): {against}; learner-confounded: {confounded}"],
            "data": {"claim_id": params["claim_id"], "recommended_status": verdict, "reason": why,
                     "narrowed_wording": params.get("narrowed_wording")},
            "evidence": [{"kind": "classification", "independence": Independence.CONSISTENCY_CHECK,
                          "summary": why, "supports": "contradicts" if against else "neutral",
                          "reviewer": "researchos.adapters.classify_claim", "method": f"rule: gap > {margin} R²",
                          "metrics": {c["context"]: round(c["gap_a_minus_b"], 4) for c in comps}}]}


@adapter("adversarial_review", "read-only", role="reviewer")
def adversarial_review(ctx, params, upstream):
    """Deterministic objections from the recorded comparisons (no language model)."""
    comps = upstream[params["from_step"]]["data"]["comparisons"]
    obj = []
    if any(not c["same_learner"] for c in comps):
        obj.append(("high", "Learner confounding: each set picked its own family by inner CV, so an R² gap is not a "
                            "representation effect. A fixed-learner comparison is needed."))
    if any(c["fit"] == "single fit" for c in comps):
        obj.append(("high", "The external result is one fit per set. The bootstrap covers test sampling only, not "
                            "training variability or learner selection."))
    if all(c["uncertainty_on_difference"] is None for c in comps):
        obj.append(("high", "No interval on the a − b difference exists anywhere; the gap's sign is not established."))
    obj.append(("medium", "The question was asked after seeing the results (post hoc). Any follow-up test needs a "
                          "protocol frozen before fitting and must not be tuned to remove the gap."))
    obj.append(("medium", "Source-disjoint is not independent: OPERA lacks a reference field, so shared original "
                          "measurements cannot be excluded."))
    obj.append(("medium", "Alternative explanation to test, not assume: non-degree indices (distance, Hosoya Z, "
                          "entropy) may extrapolate better to unseen formulas with very different sizes."))
    obj.append(("low", "Composition alone is weak on unseen external formulas (R² ≈ 0.23); differences between "
                       "two composition-augmented sets there may be unstable."))
    return {"status": ResultStatus.FLAG.value, "findings": [f"[{s}] {t}" for s, t in obj],
            "data": {"objections": [{"severity": s, "text": t} for s, t in obj]},
            "evidence": [{"kind": "review", "independence": Independence.CONSISTENCY_CHECK,
                          "summary": f"{len(obj)} rule-based objections recorded", "supports": "neutral",
                          "reviewer": "researchos.adapters.adversarial_review", "method": "deterministic rules"}]}


@adapter("fold_level_contrast", "analysis-only", path_params=["folds"])
def fold_level_contrast(ctx, params, upstream):
    """Corrected repeated-CV t interval for set_a − set_b from an archived folds.csv (no refit).
    Learner choice is still per set, so this bounds the selected-model gap, not a representation effect."""
    from statistics import mean, variance
    rows = _read_csv(ctx.cfg, params["folds"])
    a, b = params["set_a"], params["set_b"]
    fa = {(r["repeat"], r["fold"]): r for r in rows if r["set"] == a}
    fb = {(r["repeat"], r["fold"]): r for r in rows if r["set"] == b}
    keys = sorted(fa)
    if not keys or set(keys) != set(fb):
        return {"status": ResultStatus.ERROR.value, "findings": ["folds of the two sets do not match"]}
    d = [_f(fa[k]["r2"]) - _f(fb[k]["r2"]) for k in keys]
    n1 = mean(_f(fa[k]["n_train"]) for k in keys); n2 = mean(_f(fa[k]["n_test"]) for k in keys)
    J = len(d); m = mean(d); s2 = variance(d)
    se = math.sqrt((1 / J + n2 / n1) * s2)
    t975 = {14: 2.1448, 4: 2.7764, 9: 2.2622, 24: 2.0639}.get(J - 1)
    if t975 is None:
        return {"status": ResultStatus.ERROR.value, "findings": [f"no t quantile tabulated for {J - 1} df"]}
    lo, hi = m - t975 * se, m + t975 * se
    return {"status": ResultStatus.PASS.value,
            "findings": [f"{a} − {b}: ΔR² {m:+.4f} [{lo:+.4f}, {hi:+.4f}] (corrected t, {J} folds; selected learners)"],
            "data": {"mean": m, "lo": lo, "hi": hi, "folds": J},
            "evidence": [{"kind": "metric", "independence": Independence.INDEPENDENT_RECOMPUTATION,
                          "summary": f"fold-level corrected-t interval {m:+.4f} [{lo:+.4f}, {hi:+.4f}]",
                          "supports": "neutral", "reviewer": "researchos.adapters.fold_level_contrast",
                          "method": "Nadeau–Bengio corrected t on archived fold R²", "metrics": {"mean": m, "lo": lo, "hi": hi}}]}


@adapter("fixed_learner_experiment", "compute", path_params=["protocol_inputs"], role="experiment")
def fixed_learner_experiment(ctx, params, upstream):
    """Never executed by v0.1. The executor blocks compute steps; this body exists only so a
    mistaken direct call also fails closed."""
    raise PermissionError("compute steps are not executable in ResearchOS v0.1")


@adapter("llm_review", "analysis-only", role="reviewer", requires=["claude"])
def llm_review(ctx, params, upstream):
    """Optional language-model review. v0.1 never sends repository content to an external model;
    it reports the adapter as unavailable and the deterministic reviewer stands in."""
    return {"status": ResultStatus.NOT_RUN.value,
            "findings": ["LLM review adapter is not enabled in v0.1 (no external calls); deterministic review used"]}


def describe_compute(cfg, params) -> dict:
    """Dry-run description of the fixed-learner experiment: protocol draft, resources, input checks."""
    checks = []
    for rel in params.get("protocol_inputs", []):
        checks.append(ev.check_artifact(cfg, rel))
    return {"protocol_draft": params.get("protocol_draft", ""), "resource_estimate": params.get("resource_estimate", ""),
            "input_checks": checks}


def git_head(cfg) -> str:
    return inspectors.git(cfg.repo, "rev-parse", "HEAD").strip()
