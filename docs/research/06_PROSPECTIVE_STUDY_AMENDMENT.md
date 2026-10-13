# 06 — Amendment AMD-2026-10-13-E3 (exploratory): formula-controlled isomer benchmark

**Status:** locked by the git commit that adds this file, **before** any model is fitted to
real targets in this study. **Exploratory**, not confirmatory, not preregistered. Base:
branch `research-2026-10-13` (from local `audit-2026-10-12`, not pushed).

**Disclosure of prior information.** Melting point is chosen as the primary endpoint for two
reasons. It has the most isomer data (Gate C census: 758 formula groups with ≥ 5 constitutions,
6587 molecules). The audit's E2 run on approved drugs also showed MP as the only endpoint with a
consistently positive, non-significant ΔR²(TI+composition − composition). The second reason
makes the choice **outcome-informed**. The Bradley set contains the approved drugs used in E2,
so this study is **not** an independent replication of E2. Only its molecules outside the
drug set are new.

## 1. Data (target-blind selection rule)

* Endpoints: `mp_bradley` (primary) and `logp_opera` (secondary), via
  `topoqspr.expdata.get_target(name, offline=True)` with default rules (standardised parent,
  full InChIKey aggregation, median value, parents with replicate SD above the target's
  disagreement threshold dropped). Source files are the SHA-256-pinned cached copies.
* Structures: the record's `canonical_smiles`, re-standardised by
  `experiments.prepare_structures` (min 3 heavy atoms). Rejected parents are counted and excluded.
* **Inclusion (target-blind):** parents whose molecular formula (RDKit `CalcMolFormula` of
  the standardised parent) has ≥ 5 distinct constitutions (InChIKey first block) among the
  endpoint's measured parents. Stereo variants of one constitution stay as separate rows,
  grouped together in CV.

## 2. Feature sets (one matrix per endpoint; nothing fitted outside folds)

| name | definition |
|---|---|
| `size` | n_heavy, n_bonds (`features.size_features`) |
| `composition` | 11 counts incl. n_H (`features.composition_features`) |
| `TIa` | asinh of all indices of the **current registry (104)**, `experiments.index_table`; registry names and hash recorded |
| `TIa_composition` | TIa + composition |
| `rdkit2d_asinh` | asinh of `features.rdkit2d_features(target=None)` (no exclusion: MP and logP are not computed by any descList entry; Crippen MolLogP is a fitted logP estimator and is kept, F11) |
| `rdkit2d_asinh_TIa` | rdkit2d_asinh + TIa |
| `ECFP` | Morgan count, radius 2, 2048 bins (as in E2) |
| `ECFP_TIa` | ECFP + TIa |

NaN feature values are median-imputed inside the training fold by the model pipeline, as in the package.

## 3. Estimands and splits (outer 5-fold; inner 3-fold; groups respected in both)

* **E3a, known formula, unseen constitution:** groups = constitution. Formulas appear in train
  and test, constitutions never do.
* **E3b, unseen formula:** groups = molecular formula. No test formula is seen in training. No
  target-derived formula mean is used anywhere.
* E3c (similarity transfer) is **NOT RUN** in this amendment.

## 4. Learners and budget

`topoqspr.evaluate.nested_cv` with the six package families {ridge, elastic net, PLS, kNN,
SVR, extra trees}, default grids, `random_state=0`, `n_jobs=14`. Outer seeds: **0–4 if the
smoke test's projected total wall time ≤ 10 h, else 0–2** (rule fixed now). Run order: MP E3a,
MP E3b, logP E3a, logP E3b. If the budget runs out, later cells are reported NOT RUN, never
selected.

## 5. Hypotheses, metrics, decision rules

**Primary (single test, α = 0.05):** MP, E3a. ΔR² = R²(`TIa_composition`) − R²(`composition`),
fold-mean over paired outer folds. Corrected repeated-CV t-test (Nadeau–Bengio), two-sided.
"Topology adds measurable within-formula information for MP" is written only if p < 0.05 and
ΔR² > 0. "Practically meaningful" needs ΔR² ≥ 0.05, the frozen plan's threshold.
"Practically equivalent" needs the 90% CI inside [−0.05, +0.05]. Otherwise the result is
"inconclusive".

Why this contrast isolates topology: composition is constant within a formula. Under E3a the
composition model can predict only a formula-level value, so any gain of TIa_composition is
within-formula (constitutional) information.

**Secondary (BH over all secondary contrasts at q = 0.05; descriptive otherwise):** for every
endpoint × estimand cell actually run:
S1 `TIa_composition − composition` (except the primary cell); S2 `ECFP_TIa − ECFP`;
S3 `rdkit2d_asinh_TIa − rdkit2d_asinh`; S4 `ECFP − TIa_composition`; S5
`rdkit2d_asinh − TIa_composition`. Plus fold-mean R² and MAE for all sets, ΔMAE for every
contrast, and a descriptive **within-formula R²** = 1 − Σ(y − ŷ)² / Σ(y − ȳ_formula,test)²
over test formulas with ≥ 2 test members. The latter uses test labels only in the metric's
denominator, never in predictions.

**Negative control (run once, after the main cells):** MP E3a, y permuted **within formula
groups** (`numpy.random.default_rng(777)`), seed 0, sets `composition` and `TIa_composition`.
Expected ΔR² ≈ 0 (within-formula permutation destroys exactly the information under test).

## 6. Smoke test, checks, stopping

* Smoke: MP E3a with y fully permuted (`default_rng(12345)`), seed 0, all sets. Used only for
  code paths and timing. Numbers discarded unread.
* Integrity asserts in the run: no group straddles outer folds (checked from predictions);
  identical folds across sets; feature matrices hashed.
* One real run. Code-bug crash → fix, log, rerun everything. No design change after any
  real-target output.

## 7. Outputs

`~/topoqspr_runs/research_2026-10-13_E3/`: `cells/<endpoint>_<estimand>/{summary,paired,folds}.csv`,
`predictions.csv.gz`, `inclusion.csv` (ids, formula, constitution; no values),
`feature_hashes.json`, `run_info.json`. Entrypoint `docs/research/scripts/e3_isomer.py`.
Results: `07_RESULT_RECONCILIATION.md`.
