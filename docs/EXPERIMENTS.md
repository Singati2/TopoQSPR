# TopoQSPR v1.0 experiment runners (EXPERIMENT_PLAN §7–8)

Code: `src/topoqspr/experiments.py` (runners), `tools/run_experiment.py` (command line).
Tests: `tests/test_experiments.py` (offline; synthetic data with known answers, plus the
recorded pilot cohort). Fixture: `tests/fixtures/experiments/synthetic_study.json`.

**Status of the plan.** `docs/EXPERIMENT_PLAN.md` is a draft. These runners implement its
§8 experiments A, C, D and the §7 size-dependence analysis so that they can be debugged on the
quarantined pilot. They set no threshold, metric, target or panel. Nothing here reports a
model result for any main-panel disease or for the pooled approved-drug set on measured
targets; none was computed.

## 0. Common rules

* **Quarantine guard.** `experiments.check_modelling_allowed(mesh_ids, plan_frozen=False)`
  raises `QuarantineError` for any MeSH ID other than the pilot **D008881** (Migraine
  Disorders). An **empty or blank** list is refused as "disease unknown" (only `None`, i.e.
  no cohort at all, passes). Where it applies:
  * `exp-c` and `size-dependence` with `--mesh`: checked **before** any ChEMBL access (tested
    with an empty offline cache);
  * `--manifest`: the MeSH IDs of the included rows must all be D008881. A manifest with an
    empty `mesh_ids` column (what `cohort.build_atc_cohort` writes for the pooled
    approved-drug set) or without the column is refused before the freeze, so the pooled
    set cannot reach exp-c or size-dependence through a manifest (with a measured
    `--target-csv` or otherwise);
  * `experiment_a` / `exp-a`: the controlled cohort's MeSH IDs (`mesh_ids` argument united
    with `cohort_mesh_ids(controlled)`: `parameters['mesh_ids']` and the manifest column)
    must be exactly D008881, and the naive cohort's (the headings its free text matched,
    from its manifest) must lie in `PILOT_NAIVE_MESH_IDS` = D008881 plus its two MeSH
    children with ChEMBL rows, D020325 (Migraine with Aura) and D020326 (Migraine without
    Aura). Otherwise exp-a runs only in **counts-only** mode (`counts_only=True`,
    `--counts-only`): cohort sizes, overlaps, membership and size distributions, but no
    targets, correlations, top-k or models (EXPERIMENT_PLAN §3.3 step 2: "counts only; no
    properties, no models"). The CLI refuses a non-pilot `--mesh` without `--counts-only`
    before any retrieval; the naive cohort is checked once it is retrieved.

  `--plan-frozen` lifts the guard; use it only after the plan header records the frozen
  commit. A SMILES CSV carries no disease and is not checked: the guard is a seat belt, not
  a lock.
* **Underpowered flag.** `--min-n` / `--min-n-test` (API `min_n`, `min_n_test`; default 50 and
  10, EVALUATION_PROTOCOL §4) only set the `underpowered` flag of nested rows (n below
  `min_n`, or a smallest outer test fold below `min_n_test`); they never block a run. At
  replication scale (Experiment D, n ≈ 10–40) every row is flagged at the defaults: that is
  the intended reading (fold-mean R² from 5-row test folds is very noisy), not an error.
  Pass smaller values only to silence the flag deliberately, and report the values used
  (`summary.json` → `settings`).
* **Structures.** `prepare_structures` re-implements the pipeline's standardise stage
  (docs/PIPELINE.md §1; the stage function there is private): `graph.standardize`,
  rejections kept with their code, full-InChIKey duplicates rejected as
  `duplicate_inchikey:<kept id>`, and parents sharing the InChIKey first block labelled with
  a `connectivity_group`. Every CV keeps a connectivity group in one fold (random split:
  `balanced_group_kfold` when such a group exists, else `KFold(shuffle=True)`; scaffold split:
  members take their first member's Bemis–Murcko group).
* **Missing values.** NaN + reason everywhere; never 0. Rows with a missing target are dropped
  for that target and counted. NaN features are median-imputed inside the training folds
  (EVALUATION_PROTOCOL §2), except in the v0.x protocol, which drops such rows as v0.x did.
* **Determinism.** Fixed seeds (outer 0–4, inner `1000·seed + fold`, estimator 0, control
  draws `random_seed`); rows in input/registry order; CSVs with `\n` endings and shortest
  round-trip floats. Running a runner twice gives byte-identical CSVs (tested for D).
* **Statistics.** Fold-mean R² CIs use the corrected repeated k-fold t interval
  (`stats.corrected_ttest`; Nadeau & Bengio 2003). It is an interval for a mean of fold R²
  values and can exceed 1. BH via `stats.benjamini_hochberg`.
* **Reading the CSVs.** Empty cells are missing values; read with
  `pd.read_csv(path, keep_default_na=False, na_values=[""])`.

### Command line

```
PYTHONPATH=src python tools/run_experiment.py exp-a --mesh D008881 --naive-text migraine --out DIR [--model | --counts-only]
PYTHONPATH=src python tools/run_experiment.py exp-c (--mesh ID.. | --manifest F | --smiles-csv F)
        (--target NAME [--provider rdkit] | --target-csv F --target-column C) [--feature-set TI] --out DIR
PYTHONPATH=src python tools/run_experiment.py exp-d --spec study.json --out DIR
PYTHONPATH=src python tools/run_experiment.py size-dependence (--mesh ID.. | --manifest F | --smiles-csv F) --out DIR
```

Shared options: `--families` (inner-loop candidates; `poly2_ridge` names the degree-2 ridge
family below), `--seeds`, `--outer-splits`, `--inner-splits`, `--alpha`, `--n-jobs` (outer
folds in parallel; results identical), `--min-n`, `--min-n-test`, `--plan-frozen`,
`--cache-dir`, `--offline`,
`--release`, `--timestamp`, `--overwrite`. Defaults: `ridge pls knn` (A, C), `ridge pls linear`
(D), `linear poly2_ridge knn` (size dependence), seeds 0–4, 5 outer × 3 inner folds.

Every run writes the experiment's tables, `summary.json` (NaN → `null`) and
`run_info.json` (argv, time stamp, `provenance.git_info()`, package version, input paths and
SHA-256, ChEMBL release and licence, wall time). Files are staged in `.<out>.partial` and
moved into `--out` only on success; a non-empty `--out` is refused without `--overwrite`.
Exit codes: 0 success, 2 input error or quarantine refusal, 3 internal error.

## A. Cohort sensitivity (`experiment_a`, `exp-a`)

Input: the naive cohort (`cohort.naive_cohort(text)`, v0.x retrieval) and the controlled
cohort (`cohort.build_cohort(mesh_ids)`) of one disease. The naive cohort is used per
molecule (salt forms are separate rows, as in v0.x); the controlled cohort per curated parent.
**Both** are standardised with the v1.0 standardiser, so the naive salt forms collapse onto
their parents (`duplicate_inchikey`) and the comparison isolates *which drugs* each
definition retrieves. (The effect of feeding salts to the index code, the other v0.x defect,
is documented in the audit and not repeated here.)

| output | content |
|---|---|
| `summary.json` | `mode` (`full` or `counts_only`), `quarantine` (why counts-only was required, else empty), `mesh_ids` of both cohorts; `counts` per cohort (manifest rows included, distinct parent ChEMBL IDs, standardised, rejections by code, connectivity groups); `jaccard`: `parent_chembl_id`, `structures_used` (the ChEMBL `standard_inchi_key` of the structures each retrieval returns: naive salts vs controlled parents; COHORT_DEFINITION §5's InChIKey column), `standardised_parent_inchikey` (after v1.0 standardisation); `intersection_standardised`; `controlled_in_naive_share` |
| `membership.csv` | one row per standardised parent InChIKey: naive id, controlled id, membership flags |
| `distributions.csv` | per variable (`n_heavy`, `n_bonds`, `n_rings`, each target): n, mean, SD, median, quartiles, range per cohort; two-sample KS statistic and p (`scipy.stats.ks_2samp`), BH `ks_q` across the variables. The two samples overlap, so the p values are descriptive |
| `correlations.csv` | per cohort × target × index: Spearman ρ, Fisher-z CI, p, n, BH `q` across the indices of that cohort × target, reason if undefined |
| `top_overlap.csv` | per target: the top-k indices by \|ρ\| in each cohort (`all`: every defined ρ; `bh_significant`: only q ≤ α, so possibly fewer than k), the overlap count and its Jaccard. In the strict columns ties in \|ρ\| are broken by registry order. Many indices are rank-identical on small cohorts, so the `*_with_ties` columns repeat the comparison with every index tied with the k-th \|ρ\| included (realised set sizes `n_naive_with_ties`, `n_controlled_with_ties`, which can exceed k); read the tie-inclusive overlap when the two disagree |
| `model.csv` (only `model=True`) | per cohort × target: nested `compare_feature_sets({size, TI})` with identical folds and connectivity groups: R² of each set, ΔR²(TI − size) with the corrected-t CI and p, `underpowered` |
| `molecules_naive.csv`, `molecules_controlled.csv` | the standardisation tables |
| `manifest_naive.csv`, `manifest_controlled.csv` (CLI) | both cohort manifests (ChEMBL data, CC BY-SA 3.0) |

Targets come from any target provider (`--provider`, `--targets`; default the RDKit computed
targets). Before the plan is frozen, only the counts-only steps (counts, overlaps, membership,
size distributions) run for a non-pilot disease; targets, correlations, top-k and `--model`
are quarantined (§0). In counts-only mode `correlations.csv`, `top_overlap.csv` and
`model.csv` are not written and `distributions.csv` holds the size variables only.

On the pilot, exp-a needs the naive `migraine` query (`drug_indication` with
`mesh_heading__icontains=migraine` and its molecule records), which the test fixture cache
(`tests/fixtures/pipeline/chembl`) does not contain: run it online once with a writable
`--cache-dir` (it can be seeded with a copy of the fixture cache); later runs can be
`--offline` against that cache. Offline against the fixture cache alone it exits 2 with a
`CacheMiss`.

## C. Optimism chain (`experiment_c`, `exp-c`)

Input: a feature table X, a target y, and optional `connectivity` labels and `scaffolds`
(or `mols`). The CLI builds X with `experiments.feature_set` (`size`, `size_rings`,
`composition`, `TI`, `TI_<family>`, `rdkit2d` with the target's exclusion rule) from a
cohort and a target.

| step | estimate |
|---|---|
| `train` | `evaluate.legacy_protocol`: best v0.x family (by its own CV), refitted on all rows, R² on the same rows (`r2_train`) |
| `train_selected` | secondary: v1.0 inner-CV selection among `families` on all rows (3 folds, seed 0), refitted on all rows, in-sample R² |
| `cv_legacy` | `legacy_protocol` fold-mean R² of the best family (tuned and chosen on all rows; `r2_cv`) |
| `nested_random` | `evaluate.nested_cv`, random outer folds, connectivity groups kept together (`r2_nested`) |
| `nested_scaffold` | `nested_cv` with Bemis–Murcko groups (acyclic molecules singletons; `r2_scaffold`) |

B_train = r2_train − r2_nested, B_CV = r2_cv − r2_nested, B_scaffold = r2_nested − r2_scaffold
(fold-mean R²; `pooled` column: same with the mean pooled R²). CIs: B_train and B_CV treat
r2_train and r2_cv as fixed numbers (no sampling distribution is defined for them); their CI is
the nested corrected-t interval reflected. B_scaffold combines the two corrected standard
errors as if the two CV estimates were independent (df = the smaller df); they share the data,
so this is approximate (likely conservative). Undefined steps are NaN with a reason (e.g. fewer
scaffold groups than folds).

**Same rows.** The v0.x protocol drops every row with a NaN feature (as v0.x did), nested CV
median-imputes them. B_train and B_CV must subtract R² values estimated on the same molecules,
so when the v0.x protocol drops rows, the nested CV is re-run on its complete-case rows (chain
step `nested_random_legacy_rows`) and B_train/B_CV use it; `B_rows` in `summary.json` and the
`rows` column of `optimism.csv` say which rows were used, and `B_train_all_rows` /
`B_CV_all_rows` keep the all-rows contrast for comparison. B_scaffold compares two all-rows
nested runs. If fewer than 10 complete-case rows remain, the v0.x protocol does not run:
`r2_train`, `r2_cv`, B_train and B_CV are NaN with `legacy_reason`.

| output | content |
|---|---|
| `chain.csv` | the five R² values with CIs where defined (six when `nested_random_legacy_rows` is needed) |
| `optimism.csv` | B_train, B_CV, B_scaffold, CI, pooled variant, CI method, rows used |
| `nested.csv` | both nested runs: n, p, fold-mean R² + CI, pooled R², MAE, failed/valid folds, smallest test fold, underpowered, chosen families |
| `legacy.csv` | the v0.x per-family table |
| `folds_random.csv`, `folds_scaffold.csv` | outer folds (`NestedCVResult.folds`) |
| `summary.json` | all of the above plus `scaffold` (`n_groups`, `largest_group`, `largest_group_share`, `n_singletons`), the split used, connectivity groups, `legacy_n_used`/`legacy_dropped` |

`--legacy-families` restricts the v0.x families (default all ten, including the MLP).

## D. Replication harness (`experiment_d`, `exp-d`)

### Study spec (JSON; YAML if PyYAML is installed)

```jsonc
{
  "name": "…", "reference": "doi:…", "notes": "…",
  "compounds": [ {"name": "drug A", "smiles": "…", "properties": {"BP": 301.2, "MR": 80.1}},
                 {"name": "drug B", "properties": {…}} ],          // no SMILES -> reported as unresolved
  "properties": ["BP", "MR"],
  "indices": ["first Zagreb", "Randic", "ABC_4", {"label": "M_1"}],  // the study's labels
  "model_forms": ["linear", "quadratic", "cubic"],                    // polynomial in one index
  "reported_fits": [{"property": "BP", "index": "first Zagreb", "form": "linear", "r2": 0.91}],
                                                                      // or "r" (squared)
  "reported_index_values": {"drug A": {"first Zagreb": 106}}         // printed values to check
}
```

Unknown keys raise. Compound names are used as `str(name).strip()` everywhere (in
`compounds` and as keys of `reported_index_values`), so numbered compounds (`"name": 1`) work;
names must be unique after this normalisation. Model forms are case-insensitive and must be
`linear`, `quadratic` or `cubic` (in `model_forms` and in every reported fit). Drug names are
**not** resolved automatically: each compound needs a
SMILES written into the spec with its source (a compound without one is reported, not
guessed). Index labels are mapped by `map_index_name`: exact registry name; literature name
recorded in the registry (`ABC_4`, `GA_5`, `ND_1`–`ND_5`, …); the alias table
`INDEX_ALIASES` (case, spaces, accents and the word "index" ignored: "first Zagreb" → M1,
"Randić" → R, "harmonic" → H, "Wiener" → W, …); a unique case-insensitive registry match.
Anything else is `unmappable` with the reason.

### Steps (per property; compounds with a missing value dropped and listed)

1. **Indices** recomputed with TopoQSPR (`index_values.csv`); printed values compared
   (`index_check.csv`: agree if |diff| ≤ 0.01 or relative diff ≤ 10⁻³).
2. **In-sample fit** of the property on each mapped index for each model form
   (least squares; `fits.csv`: n, R², adjusted R², coefficients in ascending powers of the
   index, reported R², difference). Each index is fitted on the compounds where it is
   defined (n per row). Undefined when n < degree + 2 or fewer than degree + 1
   distinct index values. Every reported fit is listed in `reported_fits.csv` with
   `compared` and, if it could not be compared, the reason (property not analysed, index
   label unmappable or not among the study's indices, form not in `model_forms`, no r2/r,
   duplicate); `summary.json` → `reported_fits_unmatched` repeats the uncompared ones.
3. **LOOCV** of the same fit (`loo_q2` = 1 − PRESS/SST, `loo_rmse`, `loo_q2_f3`, `loo_ccc`).
4. **Nested CV** (`nested.csv`), outer `seeds` × `outer_splits`, connectivity groups kept
   together, no applicability domain (n is genre-sized; at the default `min_n` = 50 every
   row is flagged `underpowered`, see §0):
   * learner `genre`: the genre's own procedure made honest: index **and** form are chosen by
     inner CV (`single_index_family`: `SingleIndexPolynomial` over every column × form) and the
     choice is scored on the outer fold. Applied to `TI_study` (the mapped indices), and as
     baselines to `size` (heavy atoms, bonds), `composition` and `null`;
   * learner `families`: the baseline ladder (`null`, `size`, `composition`, `TI_study`) with
     the multivariable `families` (default ridge, PLS, OLS);
   * `d_r2_vs_size` + corrected-t CI and p: per learner, each set against `size` on identical folds;
   * scaffold-grouped nested CV of the `genre` learner on `TI_study` (`split = scaffold`, with
     the group statistics); NaN + reason if there are fewer scaffold groups than folds;
   * `q2_f3` and `ccc` of the out-of-fold predictions (mean over repeats).
5. **Chance controls** for the best in-sample R² over the study's k indices × forms
   (`controls.csv`; p = (1 + #{null ≥ observed})/(1 + B), Phipson & Smyth 2010). The
   observed statistic and every null use the same rule as `fits.csv`: each index is fitted
   on the rows where it is finite. `observed_best` / `observed_best_r2` / `observed_best_n`
   come from one computation; `observed_matches_fits` checks them against the best R² in
   `fits.csv` (`observed_reason` if not); `rows_per_index` gives the rows each study index
   used and `excluded_indices` any index with no defined fit, with the reason:
   * `pseudo_*`: best R² among k i.i.d. N(0, 1) pseudo-descriptors, pseudo-descriptor j
     missing on the same rows as study index j, `n_random` draws (descriptor-pool
     randomisation, Rücker et al. 2007);
   * `registry_*`: best R² among k indices drawn without replacement from the full TopoQSPR
     registry (every index with a defined fit on its finite rows; `registry_pool`, of which
     `registry_pool_partial_nan` are undefined for some compounds), `n_random` draws. This is
     not a chance control (registry indices correlate with size); it asks whether the study's
     index choice beats an arbitrary choice of the same number of indices;
   * `perm_*`: y-permutation, best R² over the study's own indices and forms,
     `n_permutations` draws.
6. **Q²_F3 and CCC** (secondary, for comparability with the genre). Q²_F3 (Consonni et al.
   2009) is computed as 1 − [Σ(y − ŷ)²/n_out] / [mean over training sets of Σ(y − ȳ_train)²/n_train]
   on pooled out-of-fold predictions (one external set gives the published formula). CCC is
   Lin's (1989) concordance correlation with population moments.

The fixture `synthetic_study.json` (20 small molecules, BP an exact quadratic of M1, MV an exact
linear function of R; one compound without a SMILES, one unmappable index label, one wrong
printed index value) checks that the harness reproduces an exact fit (R² = 1, coefficients
12.5, 3.2, −0.015), maps and reports labels, and that both controls give p = 1/(B + 1).

## Size dependence (`size_dependence`, `size-dependence`)

For every descriptor (default: every registry index on the cohort) and for the predictor
sets `size` (`n_heavy`, `n_bonds`) and `size_rings` (+ `n_rings`): repeated nested CV R² of
the descriptor regressed on the set, family chosen in the inner loop among `families`
(default OLS, `poly2_ridge` = ridge on all degree-2 monomials of the scaled size features,
α ∈ {10⁻⁴, 10⁻², 1, 100}, and kNN), connectivity groups kept together. Rows where the
descriptor is NaN are dropped for that descriptor (`n_missing`); a constant descriptor is NaN
with the reason.

| output | content |
|---|---|
| `per_descriptor.csv` | descriptor, family, set, n, n_missing, fold-mean R² + CI, pooled R², chosen families, underpowered, reason |
| `per_family.csv` | per index family × set: number of indices, number defined, median/mean/min/max R², counts with R² ≥ 0.9 and ≥ 0.99 |
| `molecules.csv`, `index_failures.csv` (CLI) | standardisation table; undefined index values with reasons |

For a single-fragment molecule n_rings = n_bonds − n_heavy + 1, so `size_rings` can add
information only through non-linear learners (EVALUATION_PROTOCOL §1). Test: `n_heavy` itself
gives R² = 1 on both sets.

## Pilot smoke run (Migraine Disorders, D008881; runtime only)

Wall times only (no performance number of these runs is reported or was inspected; only the
cohort counts below were read), measured 2026-10-03
on the shared lab machine under heavy load (1-min load average 18–34 on the runs below),
`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`, commit e897101. Pilot cohort:
33 controlled parents; naive `migraine` query 137 included molecules (headings D008881,
D020325, D020326). Expect large variation with machine load.

| run | settings | wall time |
|---|---|---|
| `exp-a --mesh D008881 --naive-text migraine --counts-only` | online, cache seeded from the fixture cache (fetches the naive query) | 32 s |
| `exp-a --mesh D008881 --naive-text migraine --model --offline` | default targets and all 102 indices, `ridge pls knn`, seeds 0–4, 5 × 3 folds, `--n-jobs 4` | 372 s |
| `exp-c --mesh D008881 --target MolLogP --feature-set TI --offline` | `ridge pls knn`, all ten v0.x families, seeds 0–4, 5 × 3 folds, `--n-jobs 4` | 203 s (65.6 s in an earlier, less loaded run) |
| `size-dependence --mesh D008881 --offline` (reduced smoke) | `--families linear knn --seeds 0 --outer-splits 3 --inner-splits 2 --n-jobs 4` | 280 s |
| `size-dependence --mesh D008881 --offline` (defaults) | `linear poly2_ridge knn`, seeds 0–4, 5 × 3 folds, `--n-jobs 4` | 2525 s (≈ 42 min; run on the pre-fix commit 7dab41d, whose size-dependence code is the same) |

The default size-dependence run fits 102 indices × 2 predictor sets × 25 outer folds, each
with an inner grid search, on 33 molecules; most of the time is per-fit and per-fold
overhead, not data size. Use the reduced command above as the smoke test; the default
settings are for the real panel runs.

## API notes

* No public API of `cohort`, `pipeline`, `evaluate` or `stats` was changed. Two pieces were
  re-implemented here because the pipeline's versions are private: the standardise stage
  (`prepare_structures`) and the CV-group rule (`cv_groups`). If the pipeline's rule changes,
  both must change together.
* `evaluate.optimism` uses the same B_train/B_CV definitions; `experiment_c` adds the scaffold
  step, the CIs and the group statistics.

## References

* Rücker C, Rücker G, Meringer M (2007) y-Randomization and its variants in QSPR/QSAR. *J Chem Inf Model* 47:2345–2357. doi:10.1021/ci700157b
* Consonni V, Ballabio D, Todeschini R (2009) Comments on the definition of the Q² parameter for QSAR validation. *J Chem Inf Model* 49:1669–1678. doi:10.1021/ci900115y
* Lin LI-K (1989) A concordance correlation coefficient to evaluate reproducibility. *Biometrics* 45:255–268. doi:10.2307/2532051
* Nadeau C, Bengio Y (2003) Inference for the generalization error. *Machine Learning* 52:239–281. doi:10.1023/A:1024068626366
* Phipson B, Smyth GK (2010) Permutation p-values should never be zero. *Stat Appl Genet Mol Biol* 9:39. doi:10.2202/1544-6115.1585
* Bemis GW, Murcko MA (1996) The properties of known drugs. 1. Molecular frameworks. *J Med Chem* 39:2887–2893. doi:10.1021/jm9602928
