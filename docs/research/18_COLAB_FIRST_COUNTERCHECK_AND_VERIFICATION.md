# 18. Colab-first reproducibility: countercheck of the four-arm review, and verification report

**Date:** 2026-10-17.

**Branch:** `repro/colab-first`, from `research-2026-10-13` at `861b3ce`. Separate worktree; nothing
merged; `main` untouched.

**Inputs checked:**
- the owner's Colab-first brief;
- `TOPOQSPR_REVISED_FOUR_ARM_PROTOCOL_REVIEW_2026-10-17.md`;
- docs 14 (r1) and 16;
- the E3 code, the published E3 data and the committed outputs.

**Labels:**
- `accepted`: adopted as written.
- `modified`: adopted with a change, and the reason is given.
- `rejected`: not adopted, with evidence.

None was rejected.

## 1. What was executed, and where

| item | where | result |
|---|---|---|
| unit tests `tests/test_repro_colab.py` (9) | clean Python 3.12.15 venv (pinned set) | 9 passed |
| OPERA restore from the pinned zip | same | 3/3 restored columns match the recorded hashes |
| notebook 01 logic, `reproduce_e3.py` | same | **571 PASS, 0 FAIL**: every committed E3/P0 number recomputable from artifacts reproduced at \|diff\| ≤ 1e-9 |
| E3 feature rebuild from an empty cache | same | MP internal: all hashes = reference (§2); external train: C, C+E, C+T, T, size, distance, ECFP hashes = archive; **RDKit-2D hashes differ** (§4.1) |
| `fixed_learner.py build` | same | 25 PASS, 0 FAIL (inclusion, labels, hashes, split order and sizes, group disjointness) |
| `fixed_learner.py run --mode smoke` (permuted labels, 16 units, 2 workers) | same | 16/16 verified in 10.5 min; resume re-run skips 16/16 |
| `fixed_learner.py run --mode bench` (permuted labels, 48 units, 6 workers) | same | 48/48; budget §5 |
| `robustness.py --mode smoke` | same | 3 PASS, 0 FAIL; **seed-0 refit of the archived external pipeline reproduces the archived E3 predictions exactly** (max \|diff\| 0.0 °C, composition/SVR) |
| corrected-t coverage simulation (ridge 400, SVR 100 replicates) | same | §3.3 |
| patent set, metadata only | same | CC0, md5 and size as pinned; **no download** |
| website claims check (`fix/website-claims`) | same | first run found "11 ML models" (10 implemented) → fixed in `3598c36` → 12 PASS; branch smoke test PASSED |
| notebooks 00–05 end-to-end with nbclient | same venv as a Jupyter kernel, local clone at the code SHA | all six OK at `104899a`; §6 |
| **Google Colab** | — | **not executed by the assistant: Colab unverified** |

**Not run:**
- any real-label fixed-learner unit;
- the patent download or census;
- the notebook-04 training seeds 1–4.

Each needs owner approval (doc 14 §11, doc 15 §7).

## 2. Line-by-line countercheck of the four-arm review

### §1 Representations

| recommendation | verdict | evidence |
|---|---|---|
| C = composition alone | accepted | `n_C … n_I, n_other, n_H`, p = 11. `n_P` is constant in both isomer sets (no P-containing formula has ≥ 5 constitutions), so the effective p is 10. The constant column is dropped inside every fold. |
| C+E = C + 9 edge counts + max degree | accepted | `e12 e13 e14 e22 e23 e24 e33 e34 e44 maxdeg` + C, p = 21, rank 19 (constant + one linear dependency). These are raw counts, as in E3. |
| C+T = `TIa_composition` | accepted | asinh(104 indices) + C, p = 115, standardised numeric rank 114. **3,264 missing cells** (internal) and 48 (external train) in the T block; median-imputed inside each training fold. |
| C+E+T | accepted | Built as E (raw), T (asinh), C. p = 125, rank 123. No exact duplicate columns in any arm. |
| Do not assume subsets after transforms | accepted | Column-wise, C ⊂ C+E ⊂ C+E+T and C ⊂ C+T ⊂ C+E+T. **C+E and C+T are not nested.** The raw registry has rank 38 and its degree-based indices are linear in m_ij (03a, Rada 2026). After asinh, the T block has near-full numeric rank (≈ 104), so the edge counts are no longer a linear function of the T columns. `C+T − C+E` therefore cannot answer the question "beyond edge counts". |
| machine-readable manifest | accepted | `fixed_learner.py build` writes `feature_manifest.json`: names, p, rank, constant columns, duplicates, missing cells and sha256 per arm. |

### §1 Contrasts

| recommendation | verdict | evidence |
|---|---|---|
| primary C+T − C | accepted | co-primary 1 |
| edge sufficiency C+E+T − C+E | accepted | co-primary 2; the only nested test of the question "registry beyond edge counts" |
| C+T − C+E as non-nested | accepted | secondary |
| C+E+T − C+T redundancy | accepted | secondary |
| (added) C+E − C | modified (addition) | Shows how much of the C+T gain the 10 edge numbers alone recover. No extra fits needed. |
| "predictive-value contrasts", not information-theoretic | accepted | doc 14 r2 §1, §9 |
| keep RDKit 2D / ECFP benchmarks | **modified** | They stay in the paper at full strength, from E3 (selected learners). They are **not** re-fitted in the fixed-learner run: ECFP (2,048 columns) dominated E3 compute, and the question here is internal to topology. Doc 14 §9 states that the experiment cannot show superiority over them. |

### §2 Model selection and leakage

| recommendation | verdict | evidence |
|---|---|---|
| six families, E3 splits fixed | accepted | `make_outer_splits` with the E3 seeds; inclusion, labels and order asserted against the archive |
| inner grouped CV, matched budget | accepted | Same grids for every arm. The only p-dependent rule (PLS n_components ≤ p) is inert at p ≥ 11. Same inner splits. |
| all fitted preprocessing inside training folds | accepted (verified in code) | imputer, constant/duplicate drop, scaler and target scaler are pipeline steps (`topoqspr.models._prep`). asinh is a fixed elementwise transform. |
| save train/test group ids, params, predictions, timings, failures per unit | **modified** | Each unit stores predictions with molecule ids; sha256 of the sorted train and test group sets; group counts; tuned params; the inner tuning table; seconds; status and reason. Full id lists are not duplicated per unit (760 × up to 6.7k ids): the splits are reproducible from the seeds, and `run_unit` re-asserts disjointness before every fit. |
| not prospective or confirmatory | accepted | doc 14 r2 header |

### §3 Statistics

| recommendation | verdict | evidence |
|---|---|---|
| ΔR² with interval, absolute R², MAE, RMSE | accepted | `analyse` writes all of them per contrast |
| δ = 0.02 is a post hoc candidate; calibrate | accepted | §3.1 |
| sensitivity δ alternatives, frozen primary | accepted | 0.01 / 0.03 descriptive only |
| TOST = whole 90% CI inside ±δ; too wide = inconclusive | accepted | r1 already used 90%. **My earlier suggestion to drop external equivalence was wrong.** The question is kept and labelled exploratory/underpowered (doc 14 §7). |
| superiority at δ needs one-sided LCB > δ | accepted, **r1 corrected** | r1 rule (b) ("95% LCB > 0 and point ≥ δ") was not a test of Δ > δ. r2 S2 requires LCB > δ. `classify` implements it, and a unit test pins the r1 case to "not shown". |
| bootstrap formulas; state that training is fixed | accepted | as E3. The training fit is held fixed, and this is stated in every output row. |
| training-seed and resampled-training variability | **modified** | External: extra-trees seeds 0–4 (the other families are deterministic given the inner splits). Internal: 3 repeats × 5 folds already resample training sets. No extra external training-set bootstrap (cost; noted as a limitation). |
| audit corrected-t coverage under grouping | accepted, done | §3.3 |
| repeat-level summaries | accepted | 3 seed-level means per contrast |
| co-primary estimands and learners; joint rule; positive/negative/mixed/inconclusive | accepted | Holm over 12 one-sided tests (2 × 2 × 3), separately for S1 and S2; BH elsewhere. "Consistent" only if a statement holds for both learners in all three contexts. ET/SVR are flagged as outcome-informed. |

### §4 Provenance

| recommendation | verdict | evidence |
|---|---|---|
| source independence ≠ measurement independence | accepted | Bradley aggregates PHYSPROP, Alfa Aesar, Karthikeyan and others. The external training set is the PHYSPROP-free Bradley subset; the test set is OPERA (= PHYSPROP). **3,253 of the 8,341 external test constitutions also occur in the Bradley isomer set used for the internal cells** (6,725 molecules, which include PHYSPROP-derived records), though not in the PHYSPROP-free training subset. That supports "source-disjoint training", not "independent measurements". Wording stays `source_disjoint_data`. |
| duplicates, formula overlap, known/unseen | accepted | parent-level unique ids; train/test constitutions disjoint (asserted); the known/unseen split reproduces the archive exactly |
| patent overlap, salts, decomposition, units | accepted | The census removes every constitution in Bradley (all rows, including `donotuse`) or OPERA MP/logP; salts → parent standardisation; decomposition/range exclusions are frozen in doc 15 §5 before values are read; °C range checks in notebook 00. |
| census before replication | accepted | `patent_census.py` refuses to download without `--approve-download`. The confirmation is deliberately unimplemented until GO plus a frozen protocol. |
| keep the negative within-formula R² next to the centred one | accepted | notebook 04 (d) writes both, labelled PREDECLARED / POSTHOC |

### §5 Colab

- **accepted**: see §6 for status.
- **Deviation:** the brief's notebook 02 name `02_fixed_learner_three_arm` became `02_fixed_learner_four_arm`.
- **Superseded (kept unchanged; the hash-locked data README names it):** the older `docs/research/colab/P0_full_vs_edgetypes.ipynb` cloned a branch, not a
  SHA. Its job is now done by notebook 01, which runs P0 and compares it with the committed values.

## 3. Evidence behind the protocol decisions

### 3.1 Margin calibration (δ = 0.02 R²)

| evaluation population | n | SD(y) °C | R² (C+T, E3) | RMSE °C | Δ RMSE for ΔR² = 0.01 / 0.02 / 0.03 |
|---|---:|---:|---:|---:|---|
| MP internal (E3a) | 6,725 | 88.4 | 0.657 | 51.7 | 0.75 / **1.51** / 2.26 |
| external known formula | 1,014 | 77.0 | 0.458 | 56.7 | 0.52 / **1.05** / 1.57 |
| external unseen formula | 7,327 | 101.8 | 0.471 | 74.0 | 0.70 / **1.40** / 2.10 |

- ΔRMSE ≈ ΔR²·Var(y) / (2·RMSE) at the operating point.
- **Replicate spread** of Bradley parents with ≥ 2 values (n = 4,665): median SD 0.71 °C, mean
  2.79 °C, 75th percentile 1.77 °C, 90th percentile 5.3 °C.
- The spread is heavy-tailed, so "uncertainty about the spread" is best summarised by these quantiles
  rather than a single SD.
- **Reading:** δ = 0.02 is about 1.0–1.5 °C of RMSE. That is above the typical replicate disagreement
  (median) and below the mean, so it is a defensible "practically negligible" threshold. It is still
  a choice, made after seeing E3, and doc 14 says so.

### 3.2 Feasibility of each statement (proxy: the P0 full-vs-edge intervals, before any fixed-learner fit; fixed-learner SEs will differ)

| context | SE(ΔR²) | 90% half-width | S3 (equivalence ±0.02) feasible? | S2 (> 0.02) needs Δ̂ ≳ |
|---|---:|---:|---|---:|
| MP E3a | ≈ 0.0045 | ≈ 0.008 | yes, if \|Δ̂\| ≲ 0.012 | 0.028 |
| MP E3b | ≈ 0.006 | ≈ 0.010 | yes, if \|Δ̂\| ≲ 0.010 | 0.030 |
| external unseen | 0.012 | 0.020 | **practically no** (\|Δ̂\| < 0.0003) | 0.040 |

### 3.3 Coverage of the corrected resampled t under formula grouping

Design: `robustness.coverage`. The 758 MP formulas are the population. In each replicate:
- half of them form the study sample (≈ 3,350 molecules);
- E3b-type grouped 3 × 5 CV is run, comparing C+T with C under a fixed learner;
- the target is the Monte Carlo expected fold-level ΔR² for the same train/test sizes, over random
  formula-disjoint draws from the whole population.

This is a finite-population approximation, because study samples overlap.

| learner | replicates | target ΔR² (MC SE) | mean estimate (SD) | **coverage, corrected t** (MC SE) | coverage, naive t | mean CI width corrected / naive |
|---|---:|---|---|---:|---:|---|
| ridge (α = 1) | 400 | 0.1139 (0.0009) | 0.1129 (0.0099) | **0.998** (0.002) | 0.858 | 0.067 / 0.031 |
| SVR (C = 3) | 100 | 0.1057 (0.0020) | 0.0999 (0.0124) | **0.990** (0.010) | 0.850 | 0.085 / 0.039 |

Output: `out/repro_2026-10-17/coverage_corrected_t.csv`.

**Reading:**
- In this design the corrected resampled t is **conservative**: its intervals are about 1.6–1.7 times
  as wide as nominal, with coverage ≈ 0.99 at a nominal 0.95.
- The naive t, which treats the 15 folds as independent, **under-covers** (≈ 0.85).
- **Caveat: overlapping study samples.** Two half-population samples share formulas, which lowers the
  replicate-to-replicate SD by about √0.5 relative to independent samples.
- **Allowing for it:** the corrected half-width is still ≈ 2.4 independent-sample SDs, so it stays
  conservative. The naive half-width becomes ≈ 1.1 SD, so its under-coverage is worse than shown.
- **Implications:**
  1. The E3 intervals are not anti-conservative. The positive MP findings are not artefacts of the
     inference.
  2. Internal equivalence (S3) and S2 are *harder* to reach than nominal widths suggest. §3.2 already
     uses corrected-t widths, so its feasibility statements stand.
  3. Doc 14 r2 keeps the corrected t, and reports this simulation next to the results.
- **Scope:** fixed ridge/SVR and C+T vs C only. Coverage for other learners and contrasts is assumed,
  not shown.


## 4. Failed checks, defects found, and risks

1. **RDKit-2D descriptors are not bit-identical across Python 3.10 and 3.12** with the same RDKit
   2026.03.6.
   - Only QED differs, by ≤ 1.1e-16 (one ulp), in 693 of 3,000 test molecules.
   - All 216 other descriptors are identical.
   - So the `rdkit2d_asinh*` feature hashes cannot be matched on Colab, though the values agree to
     machine precision.
   - The four fixed-learner arms are unaffected: their hashes match. Re-fitting the RDKit-2D arms on
     Colab (not planned) would need a tolerance check, not a hash check.
2. **BH family size was misstated as 42.** The committed q-values come from 41 secondary contrasts
   (27 internal + 14 external). The primary contrast is outside the family, so 42 is the total
   contrast count.
   - Recomputed q-values match exactly.
   - Wording fixed in 08, 16 and `PUBLICATION_READINESS`.
3. **The MP internal feature hashes were never archived.** `run_info.json` keeps only the last
   endpoint run.
   - They were recorded post hoc (`data/E3/feature_sha256_mp_internal.json`) in the E3 environment,
     after checking that inclusion and labels equal the archive.
   - The clean 3.12 rebuild reproduces them.
4. **The website claimed "11 ML models"; 10 are implemented.** Fixed on `fix/website-claims` (`3598c36`).
   Not merged, not deployed.
5. **Bugs found by the smoke runs, all fixed before the code commit:**
   - npz object arrays were not loadable with `allow_pickle=False`;
   - `NpzFile` cannot be pickled to joblib workers;
   - the coverage simulation lacked imputation for missing T values.
6. **`e3_isomer.physprop_free_bradley_keys` hard-coded the workstation cache path.** It now reads
   through `expdata.fetch` (same file, sha256-checked). Results are unchanged: the external hashes
   match the archive.
7. **Colab risks not testable here:**
   - Colab's Python version may move;
   - pip wheels for the pinned set may disappear;
   - Colab may preload a numpy that conflicts.

   Mitigation: all science runs in `python -I` subprocesses, so the kernel's preloaded packages are
   never used, and the setup cell stops on any version mismatch.
8. **I edited a hash-locked published file.** I added a note to `data/E3/README.md`, which is listed
   in `SHA256SUMS`.
   - Notebook 00 failed on it, as designed.
   - The edit was reverted, and the note moved to `data/E3/ADDENDUM_2026-10-17.md`.
9. **Compute:** the full run is 14.3 CPU-hours on this workstation (§5). Colab vCPUs are slower;
   assume 1.5–2×. A free Colab run therefore needs 2–3 sessions with Drive checkpointing.

## 5. Compute budget (bench, permuted labels)

CPU-hours by family and context. One unit = one core. Measured with `nice 10` while the owner's jobs
were running.

| family | external | mp_E3a | mp_E3b |
|---|---:|---:|---:|
| extra trees | 1.00 (5 seeds) | 3.40 | 3.40 |
| SVR | 0.04 | 1.68 | 1.68 |
| elastic net | 0.05 | 1.15 | 1.15 |
| kNN | 0.01 | 0.22 | 0.22 |
| PLS | 0.00 | 0.08 | 0.08 |
| ridge | 0.00 | 0.07 | 0.07 |
| **total** | **1.1** | **6.6** | **6.6** |

**Total: 14.3 CPU-hours.**
- By arm: C 1.2, C+E 1.8, C+T 5.6, C+E+T 5.7.
- mp_E3a is projected from the E3b fold.
- Timings on real labels may differ by tens of percent, mostly for SVR.

**Wall time:**
- workstation with 8 workers: ≈ 2 h;
- Colab free tier with 2 workers: ≈ 10–15 h, over 2–3 sessions;
- external alone: ≈ 1.1 CPU-hours, about 35 min on 2 workers.

## 6. Verification status (local vs Colab, kept separate)

**Code SHA checked out by the notebooks:** `104899afec01c366b9c31d11ba311aadbdd432c0`.

**Local test setup:**
- clean Python 3.12.15 venv with the pinned set;
- run as a Jupyter kernel through `nbclient`;
- repository cloned at that SHA from a local path, with the same fetch-by-SHA code and branch
  fallback as on GitHub;
- `nice 10` on a shared 16-core machine.

**Check reports:** `out/repro_2026-10-17/local_nbclient_104899a/`.

| notebook | local nbclient, empty work folder | local nbclient, persisted folder (Drive-style resume) | checks | Google Colab |
|---|---|---|---|---|
| 00 audit | **OK**, 171 s (fresh downloads) | OK, 22 s | 32 PASS / 0 FAIL | **unverified** |
| 01 E3/P0 reproduction | OK, 179 s (first run, at b2053e0) | OK, 92 s | 571 / 0 | **unverified** |
| 02 fixed learner (smoke) | OK, 3,321 s incl. feature build (first run, at b2053e0) | OK, 53 s (16/16 units verified and skipped) | build 25 / 0; smoke 16/16 units | **unverified** |
| 03 patent (metadata) | OK, 29 s | OK, 17 s | 4 / 0 | **unverified** |
| 04 robustness (smoke) | OK, 117 s | OK, 131 s | 3 / 0 | **unverified** |
| 05 website | OK, 200 s | OK, 790 s (machine loaded) | 12 / 0; smoke PASSED | **unverified** |

**History of these runs:**
- The first full run, at `b2053e0`, failed notebook 00 on two real defects: the edited hash-locked README (§4 item 8), and a private-path pattern that matched its own definition. Both
  were fixed in `104899a`, and all six were run again.
- **Full mode was never executed.** Nor were the census download and the training seeds 1–4.

**Colab acceptance checklist for the owner (about 30 min of attention):**
1. Open notebook 00 from the README link. Choose Runtime → Run all on a fresh CPU runtime.
   - Expect `== 00_environment_and_data_audit: 32 PASS, 0 FAIL ==`.
   - Note the `python` and `VERSION MISMATCHES: none` lines.
2. Do the same for notebook 01. Expect `571 PASS, 0 FAIL`.
3. Notebook 02 with `USE_DRIVE=True`: run it, disconnect the runtime during step 2, reconnect, and
   run all again.
   - Step 1 should be skipped.
   - Step 2 should report "N verified done".
   - The final `TOTAL 16/16`.
4. Notebooks 03, 04 and 05: run all. Each should end with 0 FAIL.
5. Send back the six `checks_*.json` files, or a screenshot of the last cell of each. Only then
   does this report change "Colab unverified" to "Colab verified (date, Python version)".

## 7. Owner decisions requested

- [ ] Approve or edit doc 14 r2 (§11 checklist), then name where the full run goes.
- [ ] Approve the patent census download (doc 15 §7).
- [ ] Run notebooks 00 and 01 once on a fresh Colab runtime and report the PASS/FAIL tables. That is
      the remaining step before "Colab verified" can be claimed.
- [ ] `fix/website-claims` (`3598c36`): review and decide whether to merge or deploy.
