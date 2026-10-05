# GO / NO-GO checkpoint for the main experiment

Required by `RECOMMENDATIONS.md` and `docs/EXPERIMENT_PLAN.md` §11. The main experiment starts only
when this report says **GO**. Status as of 2026-10-03, branch `v1.0-dev`.

## Decision: **NO-GO (not yet)**

The software, the index library, the cohort and data layers, the confirmatory protocol and the
experiment runners are implemented and adversarially verified. The main experiment is blocked only
by five decisions that the authors must take (below). Nothing below required
looking at main-panel model results; none were computed.

## 1. Code and definition issues

| Item | Status | Evidence |
|---|---|---|
| v0.x defects (Hosoya Z, wrong numerators, DS1 = DS2, Balaban J on salts, zero-filled failures, non-nested CV, shared model instance, non-deterministic cohort) | fixed in the v1.0 package | `AUDIT_TOPOQSPR.md`, `INDEX_VERIFICATION.md`, regression tests |
| Cross-version float reproducibility (Balaban J used built-in `sum()`) | fixed (`math.fsum`) | Docker (3.12) vs venv (3.10) differed by ≤ 2e-16 before; byte-identical `indices` output after (re-checked 2026-10-03, PIPELINE.md §5) |
| Intended definitions of the v0.x SS_/reverse variants | open, optional | would need the group's original scripts; v1.0 uses published definitions only |
| BMA / TMA | dropped | not defined in the paper the v0.x code attributed them to |

## 2. Indices

- 102 indices, each traced to a primary source (`INDEX_VERIFICATION.md`).
- Verified six ways: (1) closed forms; (2) brute-force matchings; (3) networkx; (4) Mordred (overlapping
  indices); (5) ChemGraphX (all 41 shared indices agree on 43 graphs, max relative difference
  1.3e-15, `CHEMGRAPHX_CROSSCHECK.md`); and (6) values printed in the primary papers
  (IJQC 2023 R_k paper: 21 values; SAR QSAR 2023 hybrid paper: 50/51 values).
- Errata found in a primary source: three printed values in the IJQC 2023 paper contradict its
  own definitions (recorded, not "fixed" in our code).

## 3. Cohorts

- Controlled cohorts from confirmed MeSH IDs, deterministic manifests, salt → parent, reasons for
  every exclusion (`COHORT_DEFINITION.md`); independently rebuilt by a script that does not use
  `topoqspr` (identical parent sets).
- Naive vs controlled overlap (Jaccard) 0.05–0.21 on four diseases.
- **Pilot (Migraine Disorders, quarantined):** 152 manifest rows → 33 parents; all standardise;
  pipeline ran end to end (6.4 h on a loaded 16-core machine).

## 4. Leakage audit

| Check | Result |
|---|---|
| Tuning and model selection inside outer folds | yes (`evaluate.nested_cv`; instrumented test) |
| Preprocessing inside the pipeline | yes (imputation, constant removal, scaling fitted on training folds) |
| Stereoisomers (identical graphs) split across folds | **fixed**: folds grouped by InChIKey connectivity block; reproduced and tested (citalopram/escitalopram, cetirizine/levocetirizine) |
| Duplicate outer seeds narrowing CIs | fixed: rejected at config validation |
| Pure-noise targets | nested R² −0.06 to −0.16 (3 data seeds); in-sample RF R² 0.84–0.85 (`tools/noise_sanity.py`) |

## 5. Protocol to freeze

| Element | Proposed (EXPERIMENT_PLAN) | Status |
|---|---|---|
| Primary endpoint | paired ΔR²(TI − size) and ΔMAE, corrected repeated k-fold t-test CI, BH (BY sensitivity) | fixed in plan; statistics code verified against independent implementations |
| Confirmatory targets | measured: mp_bradley (711 drugs), logp_opera (505), logd74_az (369), logs_aqsoldb (351) on the pooled approved set (1,821 parents) | **decision 2** |
| Disease panel | §3.3 rule leaves only Hypertension (85 parents); no disease × measured target reaches 50 | **decision 1** (recommended: per-disease descriptive only) |
| Baselines | null, size, size+rings, composition, TI, TI by family, RDKit 2D minus target-derivation family; same learners and tuning for all | implemented; null, size and composition enforced in config validation (opt-out flag recorded in the run record) |
| Practical significance | ΔR² ≥ 0.05 | **decision 3** |
| Evidence filter for cohorts | none / exclude DailyMed-only / ATC | **decision 3** |
| Y-randomised control, bootstrap sensitivity, AD coverage reporting, confirmatory table (BH/BY, fixed m, underpowered excluded) | in plan | **implemented** in the pipeline `protocol` stage (Phase 2, adversarially verified); pilot smoke run completed. Bootstrap defaults to on for pooled cohorts only, as planned |
| Attribution (SHAP) | was §9 of the plan | **not in v1.0** (no attribution code); deferred, plan §9 updated before freezing |
| Experiment runners A (cohort sensitivity), C (optimism chain), D (replication harness), size dependence | in plan | **implemented** (`experiments.py`, `tools/run_experiment.py`), with a quarantine guard that refuses modelling on main-panel diseases and the pooled set before freezing |
| Replication data | 10 studies extracted from open-access full text (210 drugs, 209 resolved, chondroitin sulfate unresolved; 1,045 property cells re-verified); reproducibility: exact 4, approximate 5, not reproducible 1 | done (`REPLICATION_DATA.md`) |
| Compute budget | pilot: 6.4 h for n = 33 (12 feature sets × 11 learners × 4 targets); protocol extrapolation: K = B = 200 needs ≈ 5.6 h permutations + 2.5 h bootstrap per target family even at n = 10–15; pooled n is 350–700 | **decision 4**: reduced learner list and K/B fixed at freeze, or a dedicated machine |

## 6. Replication candidates (Experiment D)

Blind, audited screen of the 60 newest of 239 candidates: 19 eligible (`REPLICATION_SHORTLIST.md`).
First five in date order proposed. **Decision 5** (co-authors confirm). All 19 model
predicted/computed properties; median n = 20.

## 7. Threats to validity

Internal:
- Small per-disease n; per-disease results can only be descriptive.
- Measured targets pool laboratories and protocols (AqSolDB, PHYSPROP); replicate disagreement
  rules exclude the worst cases but cannot remove inter-lab noise.
- Stereo-strict matching of measured values to drugs loses some drugs (first-block matching would
  add 8–26 % depending on the endpoint, e.g. AqSolDB 351 → 441); chosen for correctness.
- Corrected t-test CIs are approximate (Bates et al. 2023); bootstrap sensitivity planned.

External:
- Approved drugs are a narrow, size-restricted chemical space; conclusions may not transfer to
  screening libraries.
- Topological indices of the hydrogen-suppressed graph ignore elements and bond orders by
  construction; results speak to this descriptor class only.
- ChEMBL indication data change between releases; all results are tied to ChEMBL 37.

## 8. Reproducibility status

| Item | Status |
|---|---|
| Package, CLI, run record (commit, versions, config hash, manifest hash, seeds) | done (`PIPELINE.md`) |
| Lock file, Dockerfile, CI | done; image builds; CI runs offline tests |
| Byte-identical rerun | tested on fixtures; cross-version float issue fixed |
| Tagged release + Zenodo DOI | at submission |
| OSF preregistration | after GO, before main results |

## Decisions needed from the authors

Status after Advik's feedback (`docs/ADVIK_RECOMMENDATIONS.md`, 2026-10-03):

1. Per-disease panel: **decided** — pooled set for modelling, per-disease analyses descriptive (option (a)).
2. Confirmatory measured targets: **decided** — mp_bradley, logp_opera, logd74_az, logs_aqsoldb. Any later
   change after freezing is a dated amendment.
3. Practical-significance threshold and cohort evidence filter: **decided by the lead author on delegation
   (2026-10-04)**: ΔR² ≥ 0.05 (Ash et al. 2025: pre-specified practical significance); evidence filter =
   exclude indications supported only by DailyMed label text (removes label-mining noise such as
   lisinopril for diabetes; T2D 33 → 26, hypertension 85 → 69), with "none" and "FDA/EMA/ATC only" as
   Experiment A sensitivity analyses ("FDA/EMA/ATC only" alone is too strict: asthma 8). The pooled
   modelling set is almost unaffected because a drug needs only one approved indication.
4. Compute plan: **proposed (lead author, 2026-10-05), from a synthetic-data benchmark** (n = 500,
   p = 120, random features, n_jobs = 4; per 5-fold outer repeat): ridge 13 s, PLS 5 s, kNN 5 s, SVR 4.5 s,
   elastic net 10.5 s; extra trees 194 s, gradient boosting 371 s, random forest 438 s, XGBoost 628 s.
   - Confirmatory learners: ridge, elastic net, PLS, kNN, SVR **plus extra trees** (one tree ensemble is kept
     because the audit found tree models change whether indices beat size; extra trees is the cheapest).
   - Confirmatory feature sets: null, size, composition, TI, TI+composition, RDKit 2D; TI-family sets with the
     fast learners only (exploratory).
   - Y-randomisation K = 200 and bootstrap B = 200 (the plan's defaults) with ridge as the control learner
     and the full 5 × 5 repeated nested CV per draw; draws run in parallel (runtime.n_jobs = 14). The
     benchmark above ran under load 30–50; the machine is now idle (load ≈ 3 on 16 cores).
   - Two runs: (1) confirmatory, pooled cohort, sets null/size/composition/TI/TI_composition/rdkit2d,
     learners selected + the six families, protocol stage on; (2) exploratory, TI family sets with the
     five fast learners, no protocol stage.
5. Replication studies: **decided by the lead author on delegation (2026-10-04)**: Advik suggested studies "that show the best possible results"; choosing by
   reported performance would be selection on the outcome (and favours the most inflated claims), so it is
   not used. Rule, blind to reported R²: the studies that can be reproduced *exactly* from their published
   data (`REPLICATION_DATA.md` rows 2, 3, 4, 8): Kabeer et al. 2026 (antidepressants, n = 25), Vimala &
   Ravi Sankar 2026 (anticancer, n = 23), Almulla et al. 2026 (antifungal, n = 20), Godlin & Radha 2026
   (major depressive disorder, n = 21). Note: Kabeer and Almulla share two co-authors (Mufti, Ganati).
   Self-contradictory printed tables in rows 2 and 4 are reported as findings, and the re-analysis uses
   values recomputed from the stated definitions. This replaces the earlier date-order rule; changed before
   freezing and before any re-analysis.

Outside the experiment:
- Licence: the original author grants full flexibility on licensing (Advik, 2026-10-03); a one-line
  confirmation from his own GitHub account would make the MIT relicensing record complete.
- Before writing the paper (Advik): the website must work reliably on the v1.0 code (v1.0 has no web
  interface yet; the live site runs the v0.x prototype), and the index definitions are proof-checked
  against further independent papers (reverse degree, neighbourhood degree sum).
- Paper format: follow the structure and style of ChemGraphX (J. Cheminform. 2026) as closely as the
  target journal allows; SoftwareX word limit still to confirm.
