# ResearchOS decision record: MP external unseen-formula gap: TIa_composition vs edgetypes_composition

- Task: `T-a5cea8c5230c` · state **blocked** · attempt 0 · spec hash `a5cea8c5230c1f50`
- Repository: `~/moleculens-researchos` · branch `feature/researchos-orchestrator` · HEAD `797f72ca1faf0c88c7d797364ca885ef9c9c40bf` · working tree clean
- Generated 2026-10-16T15:30:54Z by ResearchOS 0.1.0 (deterministic; no language model, no training)
- Approved scopes: read-only

## Question

Why does TIa_composition apparently outperform edgetypes_composition on the source-disjoint unseen-formula external test (R2 0.471 vs 0.423)? What can and cannot be inferred without new training?

## Steps

| step | adapter | scope | status | input hash |
|---|---|---|---|---|
| ingest | ingest_artifacts | read-only | pass | `e46593b2dfe2` |
| learners | learner_family_audit | read-only | flag | `16bf54240a45` |
| points | verify_external_points | read-only | pass | `fcf5f87e6105` |
| classify | classify_claim | read-only | flag | `fa48042815b8` |
| fold_contrast | fold_level_contrast | analysis-only | awaiting approval | `-` |
| review | adversarial_review | read-only | flag | `b80af6fd6aa6` |
| fixed_learner | fixed_learner_experiment | compute | not_run | `compute` |

## Evidence inputs (sha256 and commit provenance)

| path | sha256 | last commit | verified |
|---|---|---|---|
| `docs/research/out/E3/external/summary.csv` | `a13dfa7ef906201b…` | `4a8cb1d0f` | yes |
| `docs/research/out/E3/external/paired.csv` | `b4232c7d61091e5b…` | `4a8cb1d0f` | yes |

## Findings

**ingest** (pass)
- info.json: listed in SHA256SUMS but not ingested (kept outside Git: not verifiable here)
- predictions.csv.gz: listed in SHA256SUMS but not ingested (kept outside Git: not verifiable here)

**learners** (flag)
- external/known_formula: TIa_composition used {'extratrees': 1}, edgetypes_composition used {'svr': 1} -> R² gap -0.026 is not a fixed-learner representation effect
- external/unseen_formula: TIa_composition used {'extratrees': 1}, edgetypes_composition used {'svr': 1} -> R² gap +0.048 is not a fixed-learner representation effect
- internal/mp_bradley_E3a: TIa_composition used {'extratrees': 3, 'svr': 12}, edgetypes_composition used {'svr': 15} -> R² gap -0.006 is not a fixed-learner representation effect
- internal/mp_bradley_E3b: TIa_composition used {'svr': 15}, edgetypes_composition used {'extratrees': 5, 'svr': 10} -> R² gap -0.001 is not a fixed-learner representation effect
- internal/logp_opera_E3a: both mixed learner families across folds (selected per fold) -> R² gap -0.003 is not a fixed-learner representation effect
- internal/logp_opera_E3b: TIa_composition used {'ridge': 1, 'svr': 14}, edgetypes_composition used {'extratrees': 1, 'svr': 14} -> R² gap -0.020 is not a fixed-learner representation effect
- no archived interval exists for the TIa_composition − edgetypes_composition contrast in any cell

| context | R² TIa_composition | R² edgetypes_composition | gap | learner TIa_composition | learner edgetypes_composition | fit |
|---|---|---|---|---|---|---|
| external/known_formula | 0.4577 | 0.4838 | -0.026 | {"extratrees": 1} | {"svr": 1} | single fit |
| external/unseen_formula | 0.4708 | 0.4232 | +0.048 | {"extratrees": 1} | {"svr": 1} | single fit |
| internal/mp_bradley_E3a | 0.6567 | 0.6627 | -0.006 | {"extratrees": 3, "svr": 12} | {"svr": 15} | repeated nested CV (fold-mean) |
| internal/mp_bradley_E3b | 0.6285 | 0.6297 | -0.001 | {"svr": 15} | {"extratrees": 5, "svr": 10} | repeated nested CV (fold-mean) |
| internal/logp_opera_E3a | 0.8279 | 0.8311 | -0.003 | {"extratrees": 9, "svr": 6} | {"extratrees": 9, "svr": 6} | repeated nested CV (fold-mean) |
| internal/logp_opera_E3b | 0.7660 | 0.7859 | -0.020 | {"ridge": 1, "svr": 14} | {"extratrees": 1, "svr": 14} | repeated nested CV (fold-mean) |

**points** (pass)
- 14 archived contrasts re-derived from per-set R²; max |diff| = 0.00e+00
- vs independent recomputation from predictions (earlier, separate script): max |diff| = 0.00e+00

**classify** (flag)
- claim C-EDGE10: inconclusive (contradicted descriptively in external/unseen_formula; every contrast is learner-confounded and none has an interval)
- in favour: ['external/known_formula', 'internal/mp_bradley_E3a', 'internal/mp_bradley_E3b', 'internal/logp_opera_E3a', 'internal/logp_opera_E3b']; against (gap > 0.01): ['external/unseen_formula']; learner-confounded: ['external/known_formula', 'external/unseen_formula', 'internal/mp_bradley_E3a', 'internal/mp_bradley_E3b', 'internal/logp_opera_E3a', 'internal/logp_opera_E3b']

**review** (flag)
- [high] Learner confounding: each set picked its own family by inner CV, so an R² gap is not a representation effect. A fixed-learner comparison is needed.
- [high] The external result is one fit per set. The bootstrap covers test sampling only, not training variability or learner selection.
- [high] No interval on the a − b difference exists anywhere; the gap's sign is not established.
- [medium] The question was asked after seeing the results (post hoc). Any follow-up test needs a protocol frozen before fitting and must not be tuned to remove the gap.
- [medium] Source-disjoint is not independent: OPERA lacks a reference field, so shared original measurements cannot be excluded.
- [medium] Alternative explanation to test, not assume: non-degree indices (distance, Hosoya Z, entropy) may extrapolate better to unseen formulas with very different sizes.
- [low] Composition alone is weak on unseen external formulas (R² ≈ 0.23); differences between two composition-augmented sets there may be unstable.

**fixed_learner** (not_run)
- compute step: requires_owner_approval and is never launched by v0.1
- see `researchos task dry-run` for the protocol draft and input checks

## Claim `C-EDGE10`: status **inconclusive** (exploratory)

> The gain is carried by 10 numbers (9 degree-pair edge counts + max degree).

Scope: as worded in PUBLICATION_READINESS.md at 8e5bc2e

| evidence | independence | numerical? | supports | method | summary |
|---|---|---|---|---|---|
| `E-e6d16df5fff8` | qualitative_llm_review | no (opinion) | supports | reading summary tables; no computation | an assistant summarised edgetypes_composition >= TIa_composition as 'carries the whole gain' |
| `E-8597377f995f` | consistency_check | yes | neutral | csv-module parse of two archived files; arithmetic only | paired.csv contrasts equal per-set R² differences in summary.csv (max |diff| 0.0e+00) |
| `E-52f5e07f1ae2` | independent_recomputation | yes | neutral | independent code path over predictions.csv.gz (held outside Git) | point estimates match an earlier recomputation from predictions by a separate script (docs/research/scripts/e3_recompute_and_failures.py) |
| `E-6088b08b4a31` | consistency_check | yes | contradicts | rule: gap > 0.01 R² | contradicted descriptively in external/unseen_formula; every contrast is learner-confounded and none has an interval |
| `E-df58fac85c08` | consistency_check | yes | neutral | deterministic rules | 7 rule-based objections recorded |

## Proposed narrowed wording (owner decides)

> Within the cross-validated cells, edge-type counts plus maximum degree match the full index set (MP E3a 0.663 vs 0.657; E3b 0.630 vs 0.629). On unseen external formulas the full set scored higher (0.471 vs 0.423) with different selected learners and no interval on the difference; whether the non-degree indices add information there is unresolved.

## Awaiting the owner (nothing below has been executed)

**fold_contrast** (analysis-only): path outside allowlist: ~/topoqspr_runs/research_2026-10-13_E3/cells/mp_bradley_E3b/folds.csv -> ~/topoqspr_runs/research_2026-10-13_E3/cells/mp_bradley_E3b/folds.csv; awaiting owner approval for scope 'analysis-only'

**fixed_learner** (compute): awaiting owner approval for scope 'compute'; compute is never launched by ResearchOS v0.1; 1 protected process(es) running; compute would wait

Protocol draft (not approved, not frozen):

> DRAFT (not frozen, not approved). Fixed-learner comparison on the E3 external design.
> 1. Data, splits and features exactly as E3-ext at protocol commit fc69883 (Bradley PHYSPROP-free training; OPERA PHYSPROP MP test; known/unseen formula subsets).
> 2. Sets: composition, TIa_composition, edgetypes_composition.
> 3. For EACH learner family in {ridge, elasticnet, pls, knn, svr, extratrees} separately: tune that family's default grid by the same grouped inner CV, refit, predict the test set. No selection across families.
> 4. Estimand per family: dR2(TIa_composition - edgetypes_composition) on unseen formulas, formula-cluster bootstrap (2000 resamples, seed fixed before running). Report all six families.
> 5. Decision (owner to confirm the margin before freezing): 'edge types carry the gain' only if the 90% CI lies inside +/-0.05 for the families that are competitive (R2 within 0.05 of the best); otherwise report the gap as unresolved or as information in the non-degree indices.
> 6. Optional second arm (more costly): MP E3b repeated nested CV with the family fixed, seeds 0-2.
> 7. Labelled post hoc: the question arose after seeing E3 results.

Resource estimate: Grounded in E3 archive timings: the external stage (10 sets, all six families, 14 workers) took 42 min wall on 2026-10-14, so arm 1 (3 sets) is roughly 15-25 min at 14 workers. Arm 2: E3b took 62,780 s for 10 sets (shared CPU), so 2-3 sets is roughly 3-5 h. During E3 the coordinating process used about 1.1 GB RSS and system memory use stayed below about 17 GB of 62 GB (observed 2026-10-14). Requires the machine to be free of the owner's jobs.

Input checks: `docs/research/scripts/e3_isomer.py` ok; `docs/research/06b_AMENDMENT_REVISION_BEFORE_RUN.md` ok

## Untrusted content

No instruction-like text found in the task file's prose.

## Limitations of this record

- ResearchOS re-derives and cross-checks archived numbers; it does not re-run models. `consistency_check` evidence shows archived files agree with each other, not that the analysis is right.
- `independent_recomputation` here refers to an earlier separate script over predictions held outside Git; ResearchOS verified the agreement, not the predictions themselves.
- The approval gate stops accidental or model-originated execution inside ResearchOS. It is not a security boundary against a local user with shell access.
- No language model was consulted unless an evidence row says `qualitative_llm_review`.

## Audit trail (latest first)

- `2026-10-16T15:30:54Z` executor task.transition T-a5cea8c5230c {"from": "verifying", "reason": "runnable steps done; awaiting owner: ['fold_contrast', 'fixed_learner']", "to": "blocked"}
- `2026-10-16T15:30:54Z` coordinator claim.status C-EDGE10 {"from": "supported", "reason": "contradicted descriptively in external/unseen_formula; every contrast is learner-confounded and none has an interval", "to": "inconclusive"}
- `2026-10-16T15:30:54Z` executor task.transition T-a5cea8c5230c {"from": "running", "reason": "all runnable steps finished", "to": "verifying"}
- `2026-10-16T15:30:54Z` executor step.result T-a5cea8c5230c/fixed_learner {"input_hash": "compute", "status": "not_run"}
- `2026-10-16T15:30:54Z` executor step.result T-a5cea8c5230c/review {"input_hash": "b80af6fd6aa6295f84fb40c56b68ed6986c030e18d883ade4ffaed7ad95336d8", "status": "flag"}
- `2026-10-16T15:30:54Z` executor step.result T-a5cea8c5230c/classify {"input_hash": "fa48042815b8e06c247ef8acb2e9108b37a715efae81cb0e96ec55a1dc2555d2", "status": "flag"}
- `2026-10-16T15:30:54Z` executor step.result T-a5cea8c5230c/points {"input_hash": "fcf5f87e61057494c9a271397a21363fd642b266b20f8f9024f452a70d4613ad", "status": "pass"}
- `2026-10-16T15:30:54Z` executor step.result T-a5cea8c5230c/learners {"input_hash": "16bf54240a452004cc4e816ecb144199b599fd70a307362ae2f9b77bdbd4bbfe", "status": "flag"}
- `2026-10-16T15:30:54Z` executor step.result T-a5cea8c5230c/ingest {"input_hash": "e46593b2dfe22fffab3ae1bbb7aba829806647453cf351352ce64bcc80e2926a", "status": "pass"}
- `2026-10-16T15:30:54Z` executor task.transition T-a5cea8c5230c {"from": "queued", "reason": "start", "to": "running"}
- `2026-10-16T15:30:54Z` executor task.transition T-a5cea8c5230c {"from": "approved", "reason": "queued", "to": "queued"}
- `2026-10-16T15:30:54Z` executor task.transition T-a5cea8c5230c {"from": "proposed", "reason": "approved scopes ['read-only']", "to": "approved"}
- `2026-10-16T15:30:54Z` cli task.propose T-a5cea8c5230c {"title": "MP external unseen-formula gap: TIa_composition vs edgetypes_composition"}
- `2026-10-16T15:30:54Z` owner-import claim.create C-EDGE10 {"status": "supported"}
