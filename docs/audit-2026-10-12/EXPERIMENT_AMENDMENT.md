# Amendment AMD-2026-10-12-E2 (exploratory): composition-conditional gain of topological indices, with ECFP and RDKit baselines, on random, scaffold and similarity-cluster splits

**Status:** locked by the git commit that adds this file, **before** the run script is executed
on real targets. This is a prospective plan for a new *exploratory* analysis. It is **not** a
preregistration. It was written after the audit's E1 run (RESULT_RECONCILIATION.md §3b) had
already shown the random- and scaffold-split results of the primary contrast. Only the
similarity-cluster split and every ECFP comparison are new, unseen information. Requested by
`docs/CLAUDE_FOLLOWUP_PROMPT_2026-10-12.md` §C (remote `main` @ `ab76668`).

Base commit at lock time: `8b9d95b` (local branch `audit-2026-10-12`).

## 1. Data (fixed, no recomputation of indices)

Archived frozen run `~/topoqspr_runs/main_confirmatory` (ChEMBL_37 pooled approved parents;
TI registry = **102** indices of `plan-freeze-v1`). SHA-256:

| file | sha256 |
|---|---|
| features.csv | 3298dc27d9d2c26c45374aafef319f37cca8b77e2335969f71ce2319dcd06ec2 |
| molecules.csv | 1dd8e419946f898627a6315d60e113c44cf805f33c2c8c72f00bcb51243ec6f6 |
| targets.csv | 0a0dc4c593689e15f9d1197fdedb58296e1cca483fb75c56666f2960c4910fc3 |

Read with `float_precision="round_trip"`. Endpoints and rows exactly as archived: mp_bradley
711, logp_opera 505, logd74_az 369, logs_aqsoldb 351 (rows with a value and `status = ok`). No
molecule is removed. Environment: venv Python 3.10.12, RDKit 2026.03.6, scikit-learn 1.7.2;
`OMP/OPENBLAS/MKL_NUM_THREADS=1`, `PYTHONHASHSEED=0`.

## 2. Feature sets

| name | definition |
|---|---|
| `size` | `size.n_heavy`, `size.n_bonds` (archived) |
| `composition` | the 11 archived `composition.*` counts |
| `TIa_composition` | asinh(each of the 102 `ti.*` columns), elementwise, plus `composition` raw |
| `TI_composition` | raw 102 `ti.*` plus `composition` (reported for continuity; not tested) |
| `ECFP` | Morgan count fingerprint, radius 2, 2048 folded bins, `rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048).GetCountFingerprintAsNumPy`, on `Chem.MolFromSmiles(molecules.smiles)` (standardised parent), chirality not used, raw counts |
| `rdkit2d_asinh` | asinh of the 217 archived `rdkit2d.*` columns |

asinh has no fitted parameters, so nothing is learned from test rows. It is a model-input
transform; no descriptor definition changes. ECFP is computed per molecule from its own
structure only.

## 3. Splits (outer 5-fold × seeds 0–4; inner 3-fold; groups respected in both)

* `random`: groups = InChIKey connectivity group (else the molecule), as archived.
* `scaffold`: Bemis–Murcko over all accepted molecules (acyclic → singleton), with
  connectivity-group members taking their first member's scaffold, as archived.
* `cluster` (**new**): Butina clustering (`rdkit.ML.Cluster.Butina.ClusterData`,
  `isDistData=True`, default `reordering=False`) of all 1821 accepted molecules (target-free) on
  Tanimoto distance (1 − similarity) of Morgan radius-2, 2048-bit **bit** fingerprints from the same
  generator. Distance threshold **0.4** (similarity 0.6), as written in
  SCIENTIFIC_RESEARCH_PLAN.md E2 before any result existed. Connectivity-group members take
  their first member's cluster. Cluster labels are saved, and the largest-cluster share is reported.

Folds come from `topoqspr.evaluate.nested_cv(..., groups=...)` (package `balanced_group_kfold`).

## 4. Learners

In-fold selection over {ridge, elastic net, PLS, kNN, SVR, extra trees} with the package grids
(`nested_cv`, default families, `random_state=0`, applicability domain kNN as archived). Same
for every set and split.

## 5. Hypotheses, statistics, decision rules

**Primary family H1 (4 tests, cluster split only):** for each endpoint,
ΔR² = R²(`TIa_composition`) − R²(`composition`), fold-mean over 25 paired outer folds.
Corrected repeated-CV t-test (`topoqspr.stats.corrected_ttest`, two-sided, 95%). Benjamini–
Hochberg over the 4 endpoints at q = 0.05. "TI adds information beyond composition for endpoint
X (cluster split)" is written only if q ≤ 0.05 and ΔR² > 0.

**Equivalence (declared margin):** ±0.05 ΔR², the frozen plan's practical threshold. "Practically
equivalent to composition" is written only if the 90% corrected-t interval lies inside
[−0.05, +0.05] (TOST at 5%). Otherwise, if no gain: "inconclusive / no detectable gain", never
"no effect". This applies to every endpoint × split for the primary contrast. Random and
scaffold are reported with this label but are not part of H1 (already observed in E1).

**Secondary (descriptive, no multiplicity claim), all splits:** ECFP − composition,
ECFP − TIa_composition, rdkit2d_asinh − TIa_composition, TIa_composition − size, and ΔMAE for each.
Fold-mean R² and MAE for every set. Worst-fold R².

## 6. Checks, smoke test, stopping rule

* Determinism check: on the random and scaffold splits, `size`, `composition`,
  `TIa_composition`, `TI_composition` and `rdkit2d_asinh` must reproduce E1's fold-mean R² to 1e-9
  (E1 names: `TI_asinh_composition` = `TIa_composition`). Failure stops the run.
* Smoke test (before the real run): the cluster split for `logd74_az` with the target **randomly
  permuted** (`numpy.random.default_rng(12345)`), seed 0 only, all sets. Its only use is to
  catch code errors. Its numbers are discarded unread except for completion and runtime.
* One real run. If it crashes on a code bug unrelated to results, fix it and rerun everything
  from scratch, logging each attempt. No design change after any real-target output is seen.
* Failed folds stay NaN (no imputation); a set with a failed fold has a NaN summary, reported as such.

## 7. Outputs

`~/topoqspr_runs/audit_2026-10-12_E2/`: `summary.csv`, `paired_tests.csv`, `folds.csv`,
`predictions.csv.gz`, `cluster_labels.csv`, `ecfp_sha256.txt`, `run_info.json` (commit, input
hashes, versions, timings). Copies of `summary.csv` and `paired_tests.csv` go to
`docs/audit-2026-10-12/out/E2_*.csv`. Entrypoint:
`docs/audit-2026-10-12/scripts/e2_cluster_ecfp.py`. Compute: 16 cores, `n_jobs=14`, estimated
2–3 h.
