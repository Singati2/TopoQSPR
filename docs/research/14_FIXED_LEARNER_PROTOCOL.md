# 14. Fixed-learner protocol, revision 2: four arms (DRAFT, not frozen, not approved)

**Status:**
- Draft for the owner. **Nothing in this protocol has been run on real labels.** Only permuted-label
  smoke and timing runs have been made (§10).
- The protocol is frozen only when the owner approves it. The freezing commit is then recorded in §11.

**Revision history:**
- r1 (861b3ce): two arms (full vs edge types) and a δ = 0.02 decision rule.
- r2 (this file): four arms and a corrected decision rule, after the owner's Colab-first brief and
  the four-arm review (2026-10-17). The point-by-point countercheck is in doc 18.

**Date:** 2026-10-17.

**Label:**
- This is an **exploratory, post hoc** experiment. It was designed after the E3 results and the
  P0 contrast had been seen.
- It is not prospective and not confirmatory, and it cannot revise any E3 number. Its estimands
  (feature set × *fixed* learner) differ from E3's (best-selected model per representation).

## 1. Questions

1. **Composition-controlled (manuscript question):** with the learner fixed, does topology add
   predictive value for melting point beyond element composition?
2. **Edge sufficiency:** with the learner fixed, does the registered index set add predictive value
   once the degree-pair edge counts are already supplied?
3. **Descriptive:** how do the compact edge-count set and the full set compare, and do edge counts
   add anything once the full registry is supplied?

All answers are **predictive-value contrasts under the stated learners and pipelines**. They do not
prove that any feature carries intrinsic, causal or information-theoretic content.

## 2. Arms

Every arm is built by the unchanged E3 code (`e3_isomer.featurise`). The new union arm is assembled in
`fixed_learner.arm_frames`.

| arm | E3 name | columns (p, MP internal) | construction |
|---|---|---|---|
| C | `composition` | element counts | as E3 |
| C+E | `edgetypes_composition` | 9 counts m_ij for degree pairs (i ≤ j ≤ 4, without (1,1)) + max degree + C | raw counts, as E3 |
| C+T | `TIa_composition` | asinh of the 104 registry indices + C | as E3 |
| C+E+T | new | the E block (raw), then the T block (asinh), then C | concatenation; same per-column values as in the other arms |

**Nesting.**
- C ⊂ C+E ⊂ C+E+T and C ⊂ C+T ⊂ C+E+T hold as column sets.
- C+E and C+T are **not** nested. The degree-based indices are linear in the m_ij only *before* the
  asinh transform, and the transform breaks that linearity.
- So `C+T − C+E` compares two different sets. It cannot show that the indices add information beyond
  the edge counts. That question is `C+E+T − C+E` (co-primary 2).
- Dimensions, numeric rank, constant and duplicate columns per arm are in `feature_manifest.json`,
  written by `fixed_learner.py build`. Doc 18 §2 summarises them.

**Preprocessing.** Each step is fitted on training rows only, inside the scikit-learn pipeline
(`topoqspr.models`):
1. median imputation;
2. removal of training-constant and exact-duplicate columns;
3. standardisation, for the scale-sensitive families;
4. target standardisation, for elastic net and SVR.

The asinh transform is a fixed elementwise function, so it has no fitted parameters.

## 3. Learners and tuning

- **Families:** ridge, elastic net, PLS, kNN, SVR and extra trees. Each is run **alone**
  (`families=[F]`), with the E3 grids (`topoqspr.models.FAMILIES`).
- **The grids are the same for every arm.** The only dependence on dimension is the PLS constraint
  n_components ≤ p, which leaves the grid unchanged for every arm here (p ≥ 11).
- **Tuning:** inner grouped 3-fold CV inside each training set, minimising MSE, with the same splits
  for every arm (seeded by 1000·seed + fold, as E3). The tuning budget is matched across arms.
- **Co-primary learners: extra trees and SVR.**
  - They are the two families E3 selected in almost every MP fold for these representations.
  - That choice is itself outcome-informed. To counter it, the other four families are run and
    reported in full, as secondary.
  - Ridge or kNN cannot be substituted after the results are seen.
- **External context:**
  - Extra trees is fitted with training seeds 0–4. The seed changes the forest and the inner splits.
  - The other five families are deterministic given the inner splits, so they are fitted with seed 0
    only. Seed 0 is primary.

## 4. Contexts, splits and data

| context | training | evaluation | groups |
|---|---|---|---|
| mp_E3a | Bradley isomer set (6,725 molecules, 758 formulas), E3 inclusion list asserted | E3 outer folds, seeds 0–2 × 5 | constitution |
| mp_E3b | same | same, formula-grouped | formula |
| external | PHYSPROP-free Bradley isomer set (2,252) | OPERA/PHYSPROP parents: known formula 1,014 (211 formulas), **unseen formula 7,327 (4,382)** | train constitutions excluded from test |

- **logP is not included.** The question is about melting point, and P0 showed edge types match or
  beat the full set on logP.
- **Fail-closed checks**, done in `build` and `run` before any fit:
  - inclusion lists, labels and formula assignments equal the E3 archive;
  - feature hashes equal the archived ones (external: `external/info.json`; MP internal:
    `feature_sha256_mp_internal.json`, recorded post hoc in the E3 environment, see doc 18);
  - the external split equals `test_split.csv`, in order and in sizes;
  - every unit checks that its train and test groups are disjoint;
  - training and test columns are identical.
- **External test labels are never loaded by the fitting code.** Only `analyse` reads them, from the
  restored public dataset.

## 5. Estimands and contrasts

For each family F, context c and arm pair (a, b): Δ = R²(a, F, c) − R²(b, F, c).

R² is defined identically in every arm: 1 − SSE/SST, with SST around the mean of the evaluation set.
Internal R² is the mean over 15 folds of fold R². External R² is pooled R² per test subset.

R² depends on the variance of y in each evaluation set, which differs between contexts (doc 18 §3).
So **MAE and RMSE are reported for both arms next to every ΔR²**.

| role | contrast |
|---|---|
| **co-primary 1** | C+T − C |
| **co-primary 2** | C+E+T − C+E |
| secondary | C+T − C+E (not nested), C+E+T − C+T (redundancy), C+E − C |

## 6. Inference

- **Internal: corrected resampled t** (Nadeau & Bengio 2003), over the 15 seed × fold scores, as in E3.
  - The 15 scores overlap and are not independent replicates. The correction partly allows for this,
    but its coverage under formula grouping is unverified.
  - Notebook 04 estimates that coverage by simulation. Its result is reported next to the intervals,
    whatever it is.
  - A seed-level summary (3 seed means) is reported as a descriptive check.
- **External:** a paired formula-cluster bootstrap. The 2,000 resamples (seed 2026) resample the test
  **formulas**, not molecules.
  - The training fit is **held fixed** during the bootstrap, so the interval covers only test-set
    sampling.
  - Training variability is reported separately: the SD of Δ over extra-trees seeds 0–4. The
    category is also reported for each seed.
- **One-sided bounds** are the ends of the two-sided 90% intervals.

## 7. Margin δ and decision rule

**About δ = 0.02 R²:**
- It is a **candidate smallest effect of interest**, proposed *after* the E3 and P0 effects were seen.
  It is not a demonstrated chemical threshold.
- Its calibration (doc 18 §3): on the evaluation populations, 0.02 R² corresponds to about
  1.0–1.5 °C RMSE at the observed operating points (MP internal 1.5, external known 1.05, unseen 1.4).
- Bradley replicate spread is median 0.7 °C, mean 2.8 °C and 90th percentile 5.3 °C (4,665 parents
  with ≥ 2 values). So 0.02 is about one typical disagreement between replicate measurements.
- The **primary** analysis uses δ = 0.02. Results at δ = 0.01 and 0.03 are reported as a sensitivity
  check only. They cannot change the primary statements.

**Statements, per (learner, context, contrast):**

| statement | requirement |
|---|---|
| S1 "adds predictive value" | H0: Δ ≤ 0 rejected, one-sided α = 0.05 (co-primary: Holm-adjusted) |
| S2 "adds at least δ" | H0: Δ ≤ δ rejected, i.e. one-sided 95% lower bound > δ (co-primary: Holm-adjusted). **A point estimate ≥ δ with a 95% CI excluding 0 is not S2.** |
| S3 "equivalent within ±δ" | TOST: the whole 90% CI inside [−δ, +δ] |
| S4 "worse" / "worse by at least δ" | mirror of S1/S2 |
| inconclusive | none of the above |

`fixed_learner.classify` implements these categories from the unadjusted intervals. The Holm
columns decide the co-primary statements.

**Multiplicity:**
- **Co-primary family**, Holm over 12 one-sided tests: 2 estimands × 2 learners (extra trees, SVR) × 3
  contexts (mp_E3a, mp_E3b, external unseen formula, seed 0). This is done separately for S1 and S2.
- Everything else (other families, secondary contrasts, external known formula, other seeds) gets BH
  q-values and is descriptive.
- S3 is an intersection-union test and needs no adjustment.

**Wording for the manuscript:**
- A statement is reported as **consistent** only if it holds for both co-primary learners in all
  three contexts.
- Otherwise it is **mixed by learner/context**: each cell is listed, and no general claim is made.
- A learner × feature interaction (a gap for one learner only) is reported as a finding.

**Feasibility, stated before running:**
- **External unseen formula:** the P0 bootstrap SE of a full-vs-edge difference was 0.012, which gives
  a 90% half-width of about 0.020. S3 there would need an estimate within about ±0.0003 of zero, so in
  practice external equivalence cannot be shown.
  - The external equivalence question is kept, but labelled **exploratory and underpowered**. The
    estimate and full CI are reported, and an inconclusive result is called inconclusive.
- **External S2** needs an estimate above about 0.04.
- **Internal:** with 90% half-widths of about 0.008 (P0), S2 and S3 can both be reached.

**Possible outcomes.** Every outcome is allowed and is reported as found: positive, negative, mixed or
inconclusive.
- No re-runs to change an outcome. A crash from a code bug is fixed, logged, and then the affected
  context is re-run in full.
- No change to the design after any real-label result.

## 8. Robustness (notebook 04, reported with the results)

- **(a)** Formula-cluster vs molecule-level bootstrap SE on the external test.
- **(b)** Coverage of the corrected-t interval under formula grouping. Ridge, 400 replicates; SVR, 100.
- **(c)** Training-seed variability of the archived E3 external pipeline (seeds 1–4). This needs the
  same approval as §11.
- **(d)** The pre-declared within-formula R² (negative for every TI representation) next to the
  post hoc centred metric, both labelled.

## 9. What this experiment cannot show

- Causal or information-theoretic novelty of any index.
- Superiority over RDKit 2D or ECFP. E3 already shows that TI adds nothing practically meaningful
  beyond RDKit 2D, and slightly hurts ECFP.
- Independence of the external labels. The external test is source-disjoint, not
  measurement-independent.
- Generalisation to other endpoints.

## 10. Compute budget and execution

**Work:** 760 units (`fixed_learner.units("full")`):
- 720 internal: 2 contexts × 6 families × 4 arms × 15 folds;
- 40 external: 4 arms × (5 families + 5 extra-trees seeds).

**Projection from the permuted-label bench run** (`fixed_learner.py run --mode bench`, then
`budget`). The exact table is in doc 18 §5.

The CPU-hours below are core-hours. Each unit runs on one core; units run in parallel over workers.

| where | workers | projected wall time |
|---|---|---|
| this workstation (16 cores, n_jobs = 8 when sharing) | 8 | see doc 18 §5 |
| Colab free tier (2 vCPU, sessions ≤ 12 h, may be cut earlier) | 2 | see doc 18 §5; needs several sessions |

**How to run it:**
- **Resume:** every unit is written atomically with a sha256 record, and a re-run skips only verified
  units. Set `USE_DRIVE=True` on Colab so units survive a disconnect.
- **Shards:** `--shard i/N` splits the units deterministically, so several sessions or machines can
  share the run. Results do not depend on which shard computed a unit.
- **Order:** external first (minutes), then mp_E3b, then mp_E3a. The internal cells run even if the
  external result is decisive, because they are declared here.
- **Machine sharing:** the run does not start while a protected owner job runs, and uses n_jobs ≤ 8
  when sharing.

## 11. Owner approval checkpoint

- [ ] four arms, co-primary estimands (C+T − C, C+E+T − C+E) and co-primary learners (ET, SVR) accepted
- [ ] statements S1–S4, Holm over the 12 co-primary tests, BH elsewhere accepted
- [ ] δ = 0.02 as a candidate threshold (post hoc, calibrated in doc 18 §3); sensitivity 0.01/0.03 descriptive only
- [ ] compute budget (doc 18 §5) and where it runs (workstation / Colab shards / HPC)
- [ ] notebook 04 seeds 1–4 (training variability) included
- [ ] frozen at commit: ____ (the commit that ticks these boxes; `fixed_learner.py run --mode full`
      refuses to start without `OUT/APPROVED_FULL_RUN` holding this SHA)
