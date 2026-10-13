# 04 — Code and statistics audit (Gate B)

Builds on `docs/audit-2026-10-12/AUDIT_FINDINGS.md` (F0–F18) and PR #2 items A1–A9. Each
row: location, test or counterexample, observed output, consequence, commit (branch
`audit-2026-10-12` → `research-2026-10-13`; nothing pushed).

## Recomputation from archived predictions (mission §6)

`scripts/recompute_from_predictions.py` (no topoqspr imports) on `main_confirmatory` (102 TI,
clean, commit 659fac3 → 61f63d8) and `main_confirmatory_protocol_svr` (104 TI, dirty, fe92837 → 1a504a9):

| Check | Observed (both runs) |
|---|---|
| Molecule seen twice in one (target, learner, set, repeat) | 0 |
| Molecule in different folds across sets/learners within a repeat | 0 |
| Distinct outer seeds | 0, 1, 2, 3, 4 |
| Fold R², MAE, RMSE vs `folds.csv` (R² with test-fold-mean SST) | max abs diff 0.0 |
| Pooled R² (mean over repeats) vs `summary.csv` | ≤ 1.1e-13 |
| TI − size corrected t mean, CI, p, BH q, BY q vs `confirmatory.csv` | ≤ 1.1e-16 (frozen), 0.0 (SVR) |
| BH-informative endpoints | logP only (both runs) |
| ΔR² ≥ 0.05 (practical) | MP and logP |

Frozen values: MP +0.0710 [+0.0043, +0.1377] q_BH 0.0758; logP +0.1529 [+0.0763, +0.2295]
q 0.00155; logD +0.0386 [−0.0307, +0.1079] q 0.349; logS −0.0125 [−0.0733, +0.0482] q 0.674.
These match the mission's quoted values exactly. Hand-checkable fixture: the corrected-t hand
example in `tests/test_stats.py::test_corrected_ttest_hand_example` (passes).

## Findings register (severity)

| ID | Location | Verdict | Severity | Test / evidence | Fix commit |
|---|---|---|---|---|---|
| A4/F2 corrected t zero variance | `stats.py:66` | CONFIRMED bug; no archived cell affected | MEDIUM | `test_corrected_ttest_identical_nonzero_differences_are_degenerate`, edge-case tests | `33e8198` |
| A3/F17 permutation fallback | `pipeline.py` `permutation_indices` | CONFIRMED in general (0.48 at 0.05, worst case); **0/6400** group draws split in the replayed archive; calibrated for the archived structure | MEDIUM | replay and calibration; `test_protocol_withholds_p_when_fallback_fires` | `90ca304` |
| A5/F16 dirty SVR run | run record | CONFIRMED; exact tree unrecoverable; `src/`, `studies/` unchanged across the window | MEDIUM | `git diff 1a504a9 ab2c30f -- src/ studies/` empty | diff archiving `90ca304` |
| A1/F12 controls ≠ selected learner | run records | LIMITATION (documented) | MEDIUM (wording) | PR #2 `test_a1_*` | manuscript |
| F1 Ridge/elastic-net explosion | archived predictions | CONFIRMED mechanism: Hosoya Z (TI sets), RDKit Ipc (rdkit2d); amended SVR bounded | HIGH for the paper's wording | `audit-2026-10-12/out/catastrophic_main_confirmatory.csv` | manuscript |
| ReZ1 ≡ n registered | `indices.py` registry vs `docs/INDEX_VERIFICATION.md` | CONFIRMED doc/registry inconsistency | LOW | 03 M3 F5 | not fixed (would change the 104 registry; author decision) |
| Hosoya Z float64 in tables | `index_table` | CONFIRMED; no effect (none of 9 affected molecules has a target) | LOW | 03 M5 | doc wording |
| HCM2 on edgeless graphs | `indices.py` | DEFINITION DIFFERENCE (NaN vs empty-sum 0), unreachable after standardisation | LOW | 03a | none needed |
| R_k for k > Δ → NaN | `indices.py` | documented domain restriction (literal formula defined) | LOW | 03b | doc |
| Float atom-order invariance | indices | values for isomorphic graphs agree to ≤ 5e-16, not bitwise | LOW | audit F9 | tests should use rel. tolerance |
| CSV float parsing | readers of archived CSVs | CONFIRMED (tooling) | LOW | audit F15 | scripts fixed |
| Group leakage | `evaluate.nested_cv` | VERIFIED none (E2: 0/23,945; archived: fold membership consistent) | — | — | — |
| v0.x defects (Hosoya > 18 bonds, Balaban on salts, numerators, shared estimator, full-data tuning) | v1.0 package | **Not re-attacked in this pass.** Covered by the regression tests listed in `docs/CHECKPOINT.md` §1 (1111 tests pass) | — | — | — |
| Website / API destructive tests | `topoqspr.web` | **NOT RUN** (deployment is v0.x anyway, F7) | — | — | — |
| Mutation and property-based tests (shuffle y, permute IDs, thread counts) | — | **NOT RUN** beyond the existing suite's determinism tests | — | — | — |

## Statistical inference notes

* Independent unit: molecule (stereo variants grouped). Repeated CV folds are dependent. The
  Nadeau–Bengio correction is a heuristic and its coverage is not exact (Bates et al. 2023). All
  intervals are approximate.
* Multiplicity in the frozen family is BH/BY over 4 endpoints, verified. E2 used BH over 4
  (cluster split). E3 uses one primary test plus a BH family over secondary contrasts (06b R5).
* "No detectable gain" ≠ equivalence. Equivalence needs a 90% CI inside ±0.05, which no
  composition contrast has reached so far.
