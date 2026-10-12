# Result reconciliation (audit 2026-10-12)

Every number below was recomputed from archived CSVs or from the audit's own run. Commit IDs
are given as *old → new* (history rewrite of 2026-10-12, `sha_map_2026-10-12.tsv`).

## 1. Runs and registries

| Run (archived under `results/expB/`) | Code commit | Split | TI registry | Control learner | Status |
|---|---|---|---|---|---|
| `main_confirmatory` | 659fac3 → 61f63d8 | random, grouped by InChIKey connectivity | **102** (frozen tag) | ridge | frozen confirmatory |
| `main_confirmatory_protocol_svr` | fe92837 → 1a504a9 | random, grouped | **104** | SVR | dated amendment 2026-10-05 |
| `main_scaffold` | 2c5e1d1 → 325ecb2 | Bemis–Murcko | **104** | — | frozen secondary |
| `exploratory_families` | 963a362 → 70f3648 | random, grouped | 104 (per family) | — | exploratory |
| audit E1 (`~/topoqspr_runs/audit_2026-10-12_sensitivity`) | audit branch | random **and** scaffold | **102** (the frozen run's feature matrix) | — | exploratory, this audit |

The audit run reproduced the archived fold-mean R² to 1e-9 for every raw set of
`main_confirmatory` (size, composition, TI, TI+composition, rdkit2d) and for the
registry-independent sets of `main_scaffold` (size, composition, rdkit2d), all four endpoints.
Output hashes (SHA-256): `summary.csv` 4a21b4ac…, `paired_tests.csv` f6b04197…,
`folds.csv` 4b35aa73…; copies in `out/sensitivity_*.csv`.

## 2. Sample sizes

| | Value | Source |
|---|---|---|
| Standardised approved parents (pooled) | 1821 | `molecules.csv` (all `status = ok`) |
| …with ≥ 1 measured endpoint | **959** | `targets.csv` |
| MP / logP / logD / logS | 711 / 505 / 369 / 351 | `targets.csv` |
| Connectivity groups | 33 groups, 67 parents | `molecules.csv` |
| Disease cohorts (Exp. A, controlled) | asthma 36, breast 35, hypertension 85, T2D 33 | `~/topoqspr_runs/expA/*/molecules_controlled.csv` |
| Largest disease × endpoint cell | hypertension × MP = **47** (< 50) | cohort ∩ `targets.csv` |

The pooled set is the only modelling set. No disease-specific predictive result exists (B5).

## 3. Headline and sensitivity numbers (fold-mean ΔR², corrected repeated-CV t, 95%)

### 3a. TI − size (the frozen headline)

| Endpoint | Frozen random, 102 | Amended random, 104 | Scaffold, 104 (frozen) | **Scaffold, 102 (audit)** |
|---|---|---|---|---|
| MP | +0.071 [+0.004, +0.138] | +0.071 [+0.005, +0.138] | +0.102 [+0.051, +0.152] | +0.102 [+0.051, +0.153] |
| logP | **+0.153 [+0.076, +0.230]** | +0.150 [+0.076, +0.224] | +0.098 [−0.050, +0.246] | **+0.100 [−0.050, +0.249]** |
| logD | +0.039 [−0.031, +0.108] | +0.039 [−0.030, +0.109] | +0.000 [−0.085, +0.086] | +0.002 [−0.086, +0.090] |
| logS | −0.013 [−0.073, +0.048] | −0.007 [−0.059, +0.045] | +0.004 [−0.057, +0.065] | +0.004 [−0.058, +0.066] |

With the registry held fixed (column 1 vs column 4), the logP gain falls by about a third under
scaffold folds and its interval includes zero. The 102 → 104 registry change moves nothing by
more than 0.006. Paper sentence `main.tex:235` ("the same size but wider") is wrong (F5).

### 3b. TI + composition − composition (the decisive question; never tested in the frozen run)

| Endpoint | Split | Raw TI: ΔR² | **asinh(TI): ΔR²** | **asinh(TI): ΔMAE** (+ = better) |
|---|---|---|---|---|
| MP | random | −0.04 [−0.60, +0.53] | +0.054 [−0.013, +0.121] | +2.09 °C [−0.00, +4.18] |
| MP | scaffold | −11.0 (explodes) | +0.092 [−0.008, +0.191] | +3.16 °C [+0.04, +6.28] |
| logP | random | −0.94 (explodes) | −0.028 [−0.232, +0.175] | +0.03 [−0.05, +0.11] |
| logP | scaffold | −14.9 (explodes) | +0.020 [−0.061, +0.102] | +0.02 [−0.06, +0.10] |
| logD | random | −0.026 [−0.104, +0.051] | −0.017 [−0.071, +0.037] | −0.03 [−0.07, +0.01] |
| logD | scaffold | −0.013 [−0.075, +0.049] | −0.014 [−0.086, +0.058] | −0.01 [−0.05, +0.04] |
| logS | random | −568 (explodes) | +0.009 [−0.055, +0.074] | +0.01 [−0.07, +0.08] |
| logS | scaffold | −527 (explodes) | +0.020 [−0.106, +0.145] | +0.04 [−0.05, +0.13] |

**Result (exploratory).** With extrapolation removed by a fit-free monotone transform, no
ΔR² interval excludes zero for any endpoint or split. Point estimates are ≤ 0.02 in absolute
value except melting point (+0.05 random, +0.09 scaffold). The one interval that excludes zero
is MP's scaffold ΔMAE (+3.2 °C, lower bound +0.04). That is 1 of 16 uncorrected intervals, so
it is a lead, not a finding. Raw TI+composition explodes in 5 of 8 cells. The scaffold split
adds two cells (MP and logP) that the frozen random split did not show.

### 3c. Does asinh help in general? (no tuning was done after seeing these)

* TI alone vs size, asinh: MP improves (+0.118 random, +0.144 scaffold vs +0.071/+0.102 raw).
  logP gets **worse** on random folds (+0.037 [−0.485, +0.560], worst fold R² −1.87; raw TI
  +0.153). The transform is not a uniform fix; with TI alone, kNN/SVR selection becomes
  unstable on logP.
* rdkit2d, asinh vs raw: fixes logS (random R² 0.632 vs −739; scaffold 0.569 vs 0.333) and is
  neutral elsewhere (|ΔR²| ≤ 0.05).
* Composition-free conclusion unchanged: rdkit2d (raw or asinh) > composition > TI for every
  endpoint and split. TI never reaches composition.

## 4. Representation floor (Proposition 1, `scripts/representation_collisions.py`)

| Endpoint | n | TI classes | TI R² ceiling | TI+comp R² ceiling | size R² ceiling |
|---|---|---|---|---|---|
| MP | 711 | 653 | 0.940 | 0.996 | 0.359 |
| logP | 505 | 471 | 0.993 | 1.000 | 0.491 |
| logD | 369 | 356 | 0.997 | 1.000 | 0.455 |
| logS | 351 | 322 | 0.976 | 0.999 | 0.637 |

Collisions are real but rare. The observed CV R² of TI (0.15–0.40) sits far below these
ceilings, so the bottleneck is learning and generalisation, not unlabeled-graph identity.

## 5. Controls (frozen vs amended), recomputed from `confirmatory.csv`

| Endpoint | Frozen ridge: perm. uninformative? | Frozen ridge bootstrap CI (ΔR²) | Amended SVR bootstrap CI |
|---|---|---|---|
| MP | **yes** (p = 0.254) | [−607 253, +0.078] | [+0.016, +0.135] |
| logP | no | [−1 662 029, +0.190] | [+0.046, +0.189] |
| logD | no | [−5 005, +0.029] | [−0.031, +0.090] |
| logS | **yes** (p = 0.970) | [−46 373, +0.093] | [−0.032, +0.056] |

The frozen control is uninterpretable because ridge extrapolates on the out-of-domain
molecules. The amended SVR control is well behaved but uses a different estimand than the
`selected` headline (F12).

## 6. Experiment D (TI_study − size, nested CV)

Godlin (n = 21): 13 of 27 fits have a 95% interval inside ±0.2. Kabeer (n = 25): 1 of 24.
Almulla (n = 20): 0 of 27. 0 of 78 significant at 0.05. Index values: 471 of 540 agree and
52 of 66 compounds have all values matching (`expD_macros.tex`). The authors attribute the
disagreements to different chemical graphs (printed vs resolved structures); this audit read the
macros but did **not** re-verify that attribution value by value (UNVERIFIED here).
