# 07 — Result reconciliation: E3 (exploratory amendment AMD-2026-10-13-E3), 2026-10-14, final 2026-10-16

**Status: E3 complete** (chain ended 2026-10-16T07:14). Every cell in the run order ran: MP E3a,
MP E3b, the external test, the negative control, S-noflt, logP E3a and logP E3b. None is NOT RUN.
Under 06b R5 the secondary BH family covers all cells actually run: 28 internal and 14 external
contrasts, 42 in total. **The q values below are final.**

**Run.** Commit `fc69883` (clean tree at start). Data: Bradley MP, PHYSPROP-free external
training set, OPERA PHYSPROP MP external test. Inclusion: n = 6,725 parents in 758 formula
groups, each group with at least 5 constitutions. Seeds 0–2 under the pre-declared smoke rule,
with 5 outer folds per seed, giving 15 folds. Learners: the six package families, selected by
inner CV. Wall time: MP E3a 59,499 s, MP E3b 62,780 s, logP E3a 27,229 s. logP: OPERA logP, n = 3,994 parents in 495 formula groups (same inclusion rule). Archive: `~/topoqspr_runs/research_2026-10-13_E3/`.
Small outputs and hashes: `docs/research/out/E3/`.

All results are labelled **exploratory**. E3 is a dated amendment that was written and committed
before any real-target fit. It was not preregistered with a third party.

## 1. Verification

| Check | Result |
|---|---|
| Independent recomputation of the fold metrics, corrected-t CIs and p values from `predictions.csv.gz`, with no topoqspr import (`scripts/e3_recompute_and_failures.py`) | MP and logP, E3a and E3b: max abs diff vs archive **0.0** (fold R² and every interval end) |
| External test: point estimates | exact (diff 0.0) |
| External test: bootstrap CIs, independent seed (99 vs 2026) | CI ends within **0.007** of the archive |
| Group integrity | assertion in-run: identical folds across sets and no group straddling folds (0 straddles), for both E3a (constitution groups) and E3b (formula groups) |
| Worker failures, failed folds | 0 failed folds in any set. No Traceback, MemoryError or BrokenProcessPool in the log |
| Resource contention | Throughput was reduced while the owner's separate jobs shared the CPU. This affects wall time only, not results (deterministic seeds) |

## 2. Primary test (single test, α = 0.05; 06 §5)

MP E3a (known formula, unseen constitution), ΔR² = R²(TIa_composition) − R²(composition):

| ΔR² | 95% CI | 90% CI | p | ΔMAE | folds positive |
|---|---|---|---|---|---|
| **+0.085** | [+0.065, +0.105] | [+0.069, +0.102] | < 1e-4 | **4.8 °C** [3.6, 6.0] lower | 15/15 |

Decision rule (06 §5): p < 0.05 and ΔR² > 0, so "topology adds measurable within-formula
information for MP" is **supported**. The point estimate is ≥ 0.05, so it is "practically
meaningful" by the frozen threshold. The 95% lower bound (0.065) is also ≥ 0.05.

## 3. Secondary contrasts: MP (BH over the full 42-contrast family, final)

Internal (repeated nested CV, corrected t):

| Contrast | E3a ΔR² [95% CI] | q | E3b ΔR² [95% CI] | q |
|---|---|---|---|---|
| TIa_composition − composition | (primary) | — | +0.106 [0.087, 0.125] | <1e-4 |
| edgetypes_composition − composition | +0.091 [0.071, 0.111] | <1e-4 | +0.107 [0.085, 0.129] | <1e-4 |
| distance_composition − composition | +0.068 [0.054, 0.081] | <1e-4 | +0.070 [0.051, 0.089] | <1e-4 |
| ECFP_TIa − ECFP | −0.010 [−0.021, +0.001] | 0.114 | **−0.015 [−0.029, −0.002]** | 0.048 |
| rdkit2d_asinh_TIa − rdkit2d_asinh | +0.002 [−0.002, +0.005] | 0.439 | +0.000 [−0.005, +0.006] | 0.948 |
| ECFP − TIa_composition | +0.137 [0.118, 0.156] | <1e-4 | +0.140 [0.101, 0.180] | <1e-4 |
| rdkit2d_asinh − TIa_composition | +0.175 [0.150, 0.199] | <1e-4 | +0.190 [0.153, 0.226] | <1e-4 |

Equivalence (90% CI inside ±0.05). Adding TI to RDKit 2D is **practically equivalent** to RDKit 2D
alone in E3a and E3b. Adding TI to ECFP is practically equivalent or worse: the 90% CIs are
inside ±0.05 and the E3b 95% CI is below 0.

External, source-disjoint (Bradley PHYSPROP-free → OPERA PHYSPROP). This is a single fit per set.
Intervals come from a formula-cluster bootstrap and cover test sampling only.

| Contrast | known formula (n = 1,014; 211 formulas) | q | unseen formula (n = 7,327; 4,382 formulas) | q |
|---|---|---|---|---|
| TIa_composition − composition | +0.072 [0.019, 0.138] | 0.010 | **+0.240 [0.203, 0.281]** | <1e-3 |
| edgetypes_composition − composition | +0.098 [0.037, 0.172] | <1e-3 | +0.192 [0.166, 0.221] | <1e-3 |
| distance_composition − composition | +0.052 [0.018, 0.100] | 0.010 | −0.028 [−0.058, +0.002] | 0.101 |
| ECFP_TIa − ECFP | −0.010 [−0.041, +0.015] | 0.550 | **−0.023 [−0.041, −0.005]** | 0.015 |
| rdkit2d_asinh_TIa − rdkit2d_asinh | +0.003 [−0.002, +0.009] | 0.298 | +0.002 [−0.002, +0.007] | 0.439 |
| ECFP − TIa_composition | +0.198 [0.114, 0.316] | <1e-3 | +0.095 [0.074, 0.117] | <1e-3 |
| rdkit2d_asinh − TIa_composition | +0.289 [0.209, 0.407] | <1e-3 | +0.202 [0.180, 0.224] | <1e-3 |

Absolute R² on unseen external formulas: composition 0.23 (barely above size, 0.22),
TIa_composition 0.47, ECFP 0.57, RDKit 2D 0.67. The large +0.24 mainly reflects composition
failing to extrapolate to 4,382 new formulas. It is not a strong TI model.

## 3b. logP (secondary endpoint; OPERA logP, n = 3,994, 495 formulas)

Fold-mean R²: composition 0.811 / 0.767 (E3a / E3b), TIa_composition 0.828 / 0.766, edgetypes_composition
0.831 / 0.786, ECFP 0.910 / 0.875, RDKit 2D 0.955 / 0.913. 0 failed folds.

| Contrast | E3a ΔR² [95% CI] | q | E3b ΔR² [95% CI] | q |
|---|---|---|---|---|
| TIa_composition − composition | +0.016 [0.001, 0.032], 15/15 folds | 0.056 | −0.001 [−0.023, +0.020] | 0.926 |
| edgetypes_composition − composition | +0.020 [0.009, 0.031] | 0.003 | +0.019 [−0.015, +0.052] | 0.321 |
| distance_composition − composition | +0.018 [0.008, 0.027] | 0.004 | +0.004 [−0.014, +0.022] | 0.677 |
| ECFP_TIa − ECFP | +0.002 [−0.001, +0.005] | 0.298 | +0.004 [−0.008, +0.017] | 0.550 |
| rdkit2d_asinh_TIa − rdkit2d_asinh | −0.001 [−0.003, +0.001] | 0.434 | −0.002 [−0.020, +0.015] | 0.836 |
| ECFP − TIa_composition | +0.082 [0.073, 0.091] | <1e-4 | +0.109 [0.060, 0.158] | <1e-3 |
| rdkit2d_asinh − TIa_composition | +0.127 [0.114, 0.139] | <1e-4 | +0.147 [0.087, 0.206] | <1e-3 |

**logP reading.** Composition alone explains most logP variance (R² 0.77–0.81). Topology's gain is
not significant after BH (E3a q = 0.056; E3b ≈ 0). It is **practically equivalent to zero**: the
90% CIs of TIa_composition − composition lie inside ±0.05 in both cells (E3a [0.004, 0.029]).
Within known formulas a small edge-type signal is detectable (+0.020, q = 0.003), but it does not
transfer to unseen formulas. Adding topology to RDKit 2D or ECFP is equivalent to zero. The
pre-declared within-formula R² is negative for every TI set and for composition (E3a TIa_composition −0.67).

## 4. Controls (reported, not tested; 06b R5)

| Control | ΔR² TIa_composition − composition | Reading |
|---|---|---|
| Negative control: y permuted within formula (rng 777), E3a, seed 0 | −0.028 [−0.052, −0.004], 0/5 folds positive | No spurious gain. With within-formula signal removed, the 104 extra features only add variance. Expected ≈ 0, observed slightly < 0 |
| S-noflt: replicate filter off, E3a, seed 0 (n = 6,916) | +0.088 [0.052, 0.123], 5/5 folds | The primary result does not depend on the replicate filter |

## 5. Pre-declared within-formula R² and the post hoc centred metric

The pre-declared within-formula R² (06 §5) is **negative for every TI set and for composition**
(E3a: composition −0.69, TIa −1.11, TIa_composition −0.34, edgetypes_composition −0.32). It is
positive only for the labelled representations (ECFP 0.19, RDKit 2D 0.34). Its denominator is centred on the within-formula mean, and its numerator includes each
prediction's formula-level bias. It therefore measures bias plus ranking, not ranking alone. It is
reported as declared.

**Post hoc (not pre-declared; `scripts/e3_within_formula_posthoc.py`):** a centred within-formula R²
(predictions and targets both centred per fold × formula):

| set | E3a | E3b |
|---|---|---|
| composition | 0.000 | 0.000 |
| distance_composition | 0.127 | 0.111 |
| TIa_composition | 0.189 | 0.184 |
| edgetypes_composition | 0.211 | 0.181 |
| ECFP | 0.534 | 0.516 |
| rdkit2d_asinh | 0.588 | 0.574 |

Unlabeled topology explains about 18–21% of within-formula MP variance. Labelled representations
explain 52–59%.

## 6. Failure analysis (06b R6, E3a, descriptive)

- 356 molecules have |error| > 100 °C under TIa_composition.
- The worst formulas are amino-acid-like: C3H7NO2 (MAE 133 °C), C6H13NO2, C6H11NO2, C4H9NO2.
  These are zwitterions whose MP records are often decomposition temperatures. That is label
  chemistry an unlabeled graph cannot represent.
- 2,096 molecules belong to pairs with isomorphic unlabeled graphs. Every unlabeled invariant must
  confuse these pairs. For them the TI gain over composition is smaller: 3.5 °C, against 5.4 °C for
  the rest. ECFP's advantage over TIa_composition is larger: 11.7 °C, against 9.5 °C.

## 7. Reconciliation with earlier experiments

| Earlier result | E3 | Consistent? |
|---|---|---|
| E2 (drug cohort, cluster split): MP TIa_composition − composition +0.069 [−0.017, +0.154], q 0.45 | E3a +0.085 [0.065, 0.105] | **Yes.** Same point estimate, and E2's interval contains E3's. E2 was underpowered (711 drugs with MP). E3 has 6,725 molecules. E2's "no detectable gain" was absence of evidence, not evidence of absence |
| E2 logP/logD/logS: no detectable gain | logP E3a +0.016 (q 0.056), E3b −0.001; practically equivalent to 0 | **Yes**, and now with equivalence shown for logP |
| E2: ECFP and RDKit 2D beat TIa_composition on MP | Same direction and similar size (+0.14 / +0.18) | Yes |
| 03 M1: 81 indices = linear in 9 edge counts | edgetypes_composition ≥ TIa_composition in E3a, E3b and external (known); slightly lower on external unseen (0.42 vs 0.47) | Yes. The non-degree indices add little or nothing; on new external formulas they add some |
| 03 M3: distance indices discriminate most (95.0% vs 87.4% of pairs) | distance_composition predicts worst of the three TI sets, and adds nothing on unseen external formulas | Discrimination ≠ prediction (C3), now on a measured endpoint |

## 8. What is not reported

- E5 (small-data advantage): NOT RUN (00 R10).
- GNN baselines: NOT RUN (00 R21).
