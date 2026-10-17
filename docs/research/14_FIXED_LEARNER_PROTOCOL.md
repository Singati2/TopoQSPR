# 14. Fixed-learner protocol: full index set vs edge-type counts (DRAFT, not frozen, not approved)

**Status:** draft for the owner. Nothing in this protocol has been run. It becomes frozen only when
the owner approves it and the commit that freezes it is recorded in §9.
**Answers:** doc 12 P0 step 3, as modified in doc 13 (all six families, not only SVR and extra trees).
**Date:** 2026-10-16.

## 1. Question

On melting point, does the full TI set add information beyond the 9 degree-pair edge counts plus
maximum degree, once both representations use the **same learner**?

`TIa_composition` = 104 registry indices + element counts. `edgetypes_composition` = 9 edge-type
counts m_ij + Δ + element counts. Both as built by `docs/research/scripts/e3_isomer.py` at the E3 commit.

## 2. What is already known (P0 step 1, no training; `out/P0/`, doc 13 §2)

| context | R² full | R² edge | Δ (full − edge) [95% CI] | selected families (full / edge) |
|---|---:|---:|---|---|
| MP E3a (15 folds, corrected t) | 0.657 | 0.663 | −0.006 [−0.015, +0.003] | SVR 12 + ET 3 / SVR 15 |
| MP E3b | 0.629 | 0.630 | −0.001 [−0.014, +0.011] | SVR 15 / SVR 10 + ET 5 |
| External, known formula (1,014; 211 formulas) | 0.458 | 0.484 | −0.026 [−0.070, +0.007] | ET / SVR |
| External, unseen formula (7,327; 4,382 formulas) | 0.471 | 0.423 | **+0.048 [+0.025, +0.072]** | ET / SVR |

The only interval that excludes zero is external unseen-formula, which is also the only context
where the two selected learners differ in a single fit. The fixed-learner experiment is designed
to separate the two explanations there.

## 3. Estimand

For each learner family F ∈ {ridge, elasticnet, pls, knn, svr, extratrees} (the six E3 families):

Δ_F = R²(TIa_composition, F) − R²(edgetypes_composition, F)

with hyperparameters tuned **inside the training data only** (inner 3-fold group CV, as in E3).
This is **feature × fixed learner**, a different estimand from E3's best-selected model per
representation. It is exploratory. It does not replace or revise any E3 number.

**Primary contrast:** external unseen-formula Δ_F for F = the family E3 selected for *either*
representation there (extra trees and SVR), reported as **two co-primary contrasts**. The other four
families are reported in full, as secondary results.

## 4. Data, splits and code (identical to E3)

| item | value |
|---|---|
| training data | Bradley MP, PHYSPROP-free parents (external) and the E3 inclusion set (internal), rebuilt by `e3_isomer.build` |
| external test | OPERA PHYSPROP MP parents not in the training constitutions, known/unseen split as `data/E3/external/test_split.csv` (sizes 1,014 / 7,327, asserted) |
| internal cells | MP E3a (groups = constitution), MP E3b (groups = formula); seeds 0, 1, 2 × 5 outer folds, as E3 |
| tuning | `evaluate.nested_cv` / `select_and_refit` with `families=[F]` (one family), same grids and inner splits as E3 |
| feature hashes | must equal `data/E3/run_info.json` and `external/info.json` (run aborts otherwise) |
| logP | not included: edge types already match or beat the full set there (P0), so it cannot explain the gap |

New code: one script, `docs/research/scripts/p0_fixed_learner.py`, which calls the E3 functions
with `families=[F]` and writes predictions plus folds in the E3 layout. It contains no new modelling code.

## 5. Inference

- **External:** paired formula-cluster bootstrap, 2,000 resamples, seed 2026, percentile 95% and 90%
  intervals, as in E3.
- **Internal:** corrected resampled t over 15 seed × fold scores, as in E3.
- **Multiplicity:** BH over all 18 Δ_F (6 families × 3 contexts: external unseen, E3a, E3b). Co-primary
  contrasts are reported with and without adjustment. External known-formula (n = 1,014) is reported,
  but it is outside the family because the 211 formulas give a wide interval.

## 6. Decision rule (pre-specified margin δ = 0.02 R²)

δ = 0.02 is the smallest gap that would justify a manuscript sentence that the non-degree indices add
information. It is less than half the gap observed for the selected models (0.048).

For each co-primary contrast (extra trees, SVR):

| outcome | rule | wording allowed |
|---|---|---|
| (a) equivalent | 90% CI inside [−0.02, +0.02] | "with a fixed learner, the 10 edge-type numbers match the full set within 0.02 R²" |
| (b) full set better | 95% CI lower bound > 0 **and** point estimate ≥ 0.02 | "the non-degree indices add information on unseen formulas (exploratory)" |
| (b′) edge types better | 95% CI upper bound < 0 **and** point estimate ≤ −0.02 | report as found |
| (c) unresolved | anything else | "unresolved"; no attribution claim |

- **If the two co-primary contrasts disagree:** report both. The conclusion is (c) unless the
  disagreement is itself the finding, i.e. the gap appears for one learner only. That would be
  evidence of a learner × feature interaction, and it is reported as such.
- **No re-runs to change the outcome.** A crash caused by a code bug is fixed, logged, and then
  everything is re-run. No design change follows any result.

## 7. Budget and resources (estimated from E3 wall times, 16 CPUs)

| stage | basis | estimate |
|---|---|---|
| external, 6 families × 2 sets | E3 external stage, 10 sets × 6 families: 42 min | 15–25 min |
| MP E3a, 6 × 2 | E3 MP E3a, 10 sets × 6 families: 59,499 s | ≈ 3–3.5 h |
| MP E3b, 6 × 2 | 62,780 s for 10 sets | ≈ 3.5 h |

Peak memory in E3 was about 1.1 GB RSS. GPU is not used.

- **Run order:** external first. It answers the primary question, and the internal cells only add context.
- **Stop rule:** the internal cells run even if external is decisive, because they were declared here.
- **Machine sharing:** the run waits while any protected owner job is running (ResearchOS
  `processes()` list), and uses `n_jobs ≤ 8` when sharing the machine.

## 8. Possible results and what each means

- **(a) on both:** N2 can be stated with a fixed learner: the edge counts suffice. The 0.048 gap was
  learner selection.
- **(b) on both:** the readiness wording changes to "the full set adds about X on unseen formulas".
  N2 is dropped and replaced by that result.
- **Interaction:** report it. It weakens any representation-only claim.
- **(c):** the claim stays "unresolved", exactly as in `PUBLICATION_READINESS.md` at `9dce4d3`.

## 9. Owner approval checkpoint

- [ ] δ = 0.02 accepted, or changed to: ____
- [ ] co-primary families (extra trees, SVR) accepted
- [ ] budget accepted; the run may start when no protected job is running
- [ ] frozen at commit: ____ (filled in by the commit that freezes it; the run log must cite it)
