# TopoQSPR — Publication and Adversarial Audit Recommendations (2026-10-12)

**Repository:** `Singati2/TopoQSPR`  
**Review reference:** `main` at `ab766682825c3cb73a5d066313fed5a927222ed2` (inspected 2026-10-12, America/New_York).  
**Status:** Recommendations and issue triage only — **not** a claim that experiments have been rerun or that every defect is independently reproduced.  
**Scope:** v1.0 pipeline, archived Experiment B and amended SVR controls, code/statistics, release provenance, and publication claims.  
**Preservation rule:** Do not rewrite the frozen protocols, source records, archived output files, original results, or negative findings.

## Executive judgment

TopoQSPR has progressed from a prototype into a substantially more reproducible graph-theoretic QSPR platform, with controlled drug manifests, verified topological indices, nested evaluation, scaffold experiments, published-study replication, and archived model results. **It is not ready for a final journal submission until the inferential contract, comparison against composition, and experiment provenance are reconciled.**

A defensible main scientific statement today is narrower than "topological indices improve prediction": the frozen, selected-learner Experiment B finds a positive TI-versus-two-feature-size contrast for **logP**, but it does not establish superiority over elemental composition, Morgan/ECFP fingerprints, or strong RDKit baselines.

### Archived, frozen Experiment B — verified from the checked-in `results/expB/main_confirmatory/confirmatory.csv`

| Measured endpoint | n | Paired mean Delta R2 (TI − size) | Corrected 95% interval | BH q (4 tests) | BH-positive? |
|---|---:|---:|---|---:|---|
| `mp_bradley` | 711 | +0.071005 | [+0.004302, +0.137708] | 0.075828 | No |
| `logp_opera` | 505 | +0.152927 | [+0.076329, +0.229526] | 0.001553 | Yes |
| `logd74_az` | 369 | +0.038623 | [-0.030688, +0.107935] | 0.348587 | No |
| `logs_aqsoldb` | 351 | -0.012533 | [-0.073262, +0.048195] | 0.673940 | No |

These are **the archived results**, not newly executed confirmations. Significance of the original paired test is not a proof of increased information beyond chemical composition, nor evidence of equivalence on the three non-positive endpoints.

## Prioritized action items, with falsification tests and acceptance criteria

### P0 — Establish incremental graph information beyond composition

**Classification:** Important missing empirical contrast, not a confirmed code defect.  
**Evidence:** `docs/INDEPENDENT_ADVERSARIAL_AUDIT_2026-10-12.md` A2; `docs/FEYNMAN_REVIEW.md`; `src/topoqspr/pipeline.py`; the frozen `confirmatory.csv` compares TI against `size`, not composition.

1. Run paired `TI_composition − composition` and `TI − composition` comparisons on exactly matched outer folds, with full inner selection and all training-dependent transformations fitted on training-only data.
2. Include Morgan/ECFP and RDKit 2D baselines with comparable tuning budgets and the same held-out compounds.
3. Archive per-fold and per-molecule predictions, split/group identities, seeds, selected learner, corrected-CV intervals, multiplicity treatment, and pre-specified or expressly exploratory status.
4. Evaluate random/grouped and scaffold-grouped generalization separately. Where molecular graphs or near-identical structures share identifiers, do not split them between train/test.
5. Report both positive and negative contrasts. Do not infer `TI adds nothing` from a CI that includes zero. Do not promote a post hoc analysis to prospective confirmation.

**Acceptance:** Every manuscript statement differentiates "beyond two molecule-size counts" from "beyond element composition" and "beyond established fingerprints"; supporting contrasts and uncertainty are reproducible.

### P0 — Repair the degenerate corrected t-test inferential contract

**Classification:** Confirmed archived implementation/reporting convention; whether this affects a published numerical conclusion remains to be quantified.  
**Evidence:** `src/topoqspr/stats.py`, `tests/test_stats.py`, and PR #2 (`tests/test_adversarial_audit_20261012.py`) treat constant nonzero fold differences as `p=0` with a zero-width confidence interval.

**Small counterexample:** `corrected_ttest([0.1, 0.1, 0.1], n_train=8, n_test=2)` currently gives `p=0`, CI `[0.1, 0.1]`. Repeated CV folds are dependent; identical observed differences alone do not demonstrate a point-identified population effect or justify an unconditional significance statement.

1. Define the inferential contract for zero estimated variance; strongly consider returning a finite descriptive mean but `p` and inferential CI as undefined, with explicit reason, instead of manufacturing certainty. Document any chosen alternative mathematically.
2. Add failing-before/passing-after tests for all-zero differences, constant nonzero differences, tiny numerical variance, missing or failed folds, and normal nondegenerate cases.
3. Review PR #2's `test_a4_degenerate_corrected_ttest_behavior_is_explicit` against the corrected, justified contract, rather than preserving an invalid expectation merely because CI currently passes.
4. Recalculate affected result tables and report before/after deltas if this case occurred in any archived inference. **Never quietly overwrite the original archived results.**

**Acceptance:** The inferential interface is explicit, mathematically defensible, robust at degeneracy, and manuscript claims trace to the applicable implementation version.

### P0 — Test grouped-permutation exchangeability, not just index bijection

**Classification:** Statistically plausible invalid null procedure for some group structures; relevance to archived real-data tests needs independent verification.  
**Evidence:** `src/topoqspr/pipeline.py` grouped permutation fallback; `docs/INDEPENDENT_ADVERSARIAL_AUDIT_2026-10-12.md` A3; PR #2's adversarial permutation diagnostic.

A returned vector that is a permutation of `y` does not establish exchangeability. A unique-size correlated group can be broken up by row-level mixing, destroying dependence.

1. Run the documented synthetic grouped-null simulation and report seed, group sizes, empirical Type I error, confidence interval, and sensitivity to several intra-group correlations.
2. Replay the **actual saved group identities and permutations** for all original and amended confirmatory cells; quantify how often the problematic fallback was exercised.
3. State and check the exchangeability assumptions for any purported permutation p-value. Under unsupported patterns, fail closed or use a defensible group-respecting null with a clear estimand.
4. If the fallback never occurred in the frozen analyses, state exactly that narrow fact; do not claim universal permutation validity.

**Acceptance:** The method has an explicit null validity argument or a fail-closed restriction; archived affected/unaffected cells are listed with evidence.

### P1 — Resolve selected-learner versus fixed-family controls

**Classification:** Interpretation/protocol scope limitation.  
**Evidence:** `results/expB/main_confirmatory/run_record.json` describes a selected learner from six families, while the original bootstrap/permutation control is Ridge-only. `results/expB/main_confirmatory_protocol_svr/` contains amended SVR-only controls; its logP bootstrap interval is approximately [+0.046, +0.189].

1. Label the main paired selected-model test and Ridge-only/SVR-only sensitivity controls with separate names and separate versioned settings in text, tables, and captions.
2. If feasible, repeat resampling for the **entire learner-selection process**, or do an independent reduced-cost replication expressly labeled as such.
3. Preserve original Ridge instability (including extreme negative outcomes) rather than replacing it in the narrative with the SVR run.
4. Never treat a fixed SVR-family bootstrap as validating the full selected-family procedure.

**Acceptance:** Every inferential number clearly identifies the model-selection process being evaluated.

### P1 — Reconcile 102/104 indices and archived build provenance

**Classification:** Provenance/release issue; not necessarily an index-implementation error.  
**Evidence:** Original `run_record.json` says 102 TI features and clean working tree; the amended SVR control reports 104 and `dirty=true`. The README/roadmap and published-feature descriptions must not blur these versions.

1. Build a claim-provenance ledger connecting **each** manuscript value to its data, registry size, code commit, clean/dirty status, study config, split hashes, and output CSV.
2. Rerun the amended SVR control from a clean, pinned tree, or identify its precise source tree as unrecoverable if missing diff files prevent reconstruction. Document changes in numerics, if any.
3. Verify the status of local-only audit commits against remote `main` and open PR #2; do not describe a local branch as pushed if it is not.
4. Align abstract, figures, table captions, `paper/main.tex` macros, and package metadata with the precise experimental software version. Do not retroactively replace frozen 102-index results with 104-index output.

**Acceptance:** No result claims to have been produced by code or indices different from those recorded in its run record.

### P1 — Reproduce numerical, mathematical, and software claims independently

**Classification:** Verification work still required, not a known aggregate failure.  
**Evidence:** GitHub Actions runs for main and PR #2 were reported successful, but an automated test suite is not the same as an independent full-study regeneration.

1. Run the full test suite and focused adversarial tests from an explicit SHA and pinned environment; archive commands, dependency versions, exit statuses, skips, timestamps, and output hashes.
2. Independently regenerate the reported primary TI index values using algebraic small-graph checks, published examples, and separate reference implementations (including ChemGraphX where shared). Flag structural equivalence, degeneracy, and numerical overflow or extrapolation.
3. Reproduce at least one complete endpoint's raw-input-to-manifest-to-indices-to-nested-CV-to-table pipeline and compare individual predictions and inferential outputs, not merely rounded headline R2.
4. Include published-paper replication evidence and negative results: distinguish reproducing an in-sample published fit from externally valid prediction; quantify optimism and small-n uncertainty.
5. Verify licensing attribution for any redistribution of ChEMBL-derived data and complete the release checklist (versioned tag, exact dependencies, manuscript metadata, and archived DOI where applicable).

**Acceptance:** An external reader can independently recover the principal scientific table and identify exactly what could not be repeated.

## Manuscript strategy

**Primary scientific narrative:** Controlled, disease-anchored cohort construction plus mathematically verified graph descriptors, reproducible nested evaluation, comparison against meaningful chemical baselines, and a critical empirical analysis of when graph indices help or fail.

**Preserve:** Null and size controls; the single BH-positive logP size contrast; three non-positive comparisons; scaffold limitations; model-extrapolation errors; replication-study optimism findings; numerical index verification; missing or non-estimable claims.

**Avoid:** Claims that all 104 indices contributed to the original 102-index analysis; universal TI superiority; equivalence from non-significance; treating overlapping CV folds as independent samples; calling an exploratory amendment a pre-registration; calling every passing pytest run a validated scientific inference.

**Journal path:** A software/methodology article for *Journal of Cheminformatics* may be justified once the statistical audit and release are complete; *SoftwareX* remains a reasonable alternative if methodological novelty is too limited. The exact venue, fit and limits should be checked against current author guidelines before submission.

## Immediate order of work for Claude / Codex

1. Inspect HEAD, branches, PR #2, original and amended run records; reconcile local-only audits without destructive merges.
2. Fix or explicitly resolve the degenerate inference contract, with a counterexample and regression tests.
3. Independently assess permutation Type I behavior and the actual frozen group arrangements.
4. Complete composition-conditional and fingerprint-baseline paired comparisons as **exploratory new analyses**, unless genuinely prospectively frozen.
5. Reconcile descriptor versions, dirty-run reproducibility, source-to-claim tables, and literal manuscript wording.
6. Run focused tests, full suite, and at least one independent end-to-end reconstruction; publish a claim-by-claim PASS / LIMITED / UNRESOLVED / FAIL report.

**Do not alter**: original frozen experimental results, protocol commitments, or negative findings. Every confirmed bug needs a minimal falsifying example, the path and function affected, test proving the correction, and an explicit numerical effect or a verified no-change result.

## Related existing repository material

- [Independent adversarial audit](./INDEPENDENT_ADVERSARIAL_AUDIT_2026-10-12.md)
- [Claude follow-up prompt](./CLAUDE_FOLLOWUP_PROMPT_2026-10-12.md)
- [Feynman review](./FEYNMAN_REVIEW.md)
- [Experiment plan and amendments](./EXPERIMENT_PLAN.md)
- [Checkpoint](./CHECKPOINT.md)
- [Earlier publication recommendations](../RECOMMENDATIONS.md)
- [PR #2: adversarial audit tests](https://github.com/Singati2/TopoQSPR/pull/2)

*This file records concrete recommendations and archived evidence as of the review. It is not a new experiment, independent rerun, acceptance decision, or alteration to a pre-registered analysis.*
