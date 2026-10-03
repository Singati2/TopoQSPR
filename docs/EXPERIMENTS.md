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
  Disorders). It is called by `experiment_a(model=True)` and, in the command line, by
  `exp-a --model`, `exp-c` and `size-dependence` on a `--mesh` or `--manifest` input,
  **before** any ChEMBL access (tested with an empty offline cache). `--plan-frozen` lifts
  it; use it only after the plan header records the frozen commit. A SMILES CSV carries no
  disease and is not checked: the guard is a seat belt, not a lock.
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
PYTHONPATH=src python tools/run_experiment.py exp-a --mesh D008881 --naive-text migraine --out DIR [--model]
PYTHONPATH=src python tools/run_experiment.py exp-c (--mesh ID.. | --manifest F | --smiles-csv F)
        (--target NAME [--provider rdkit] | --target-csv F --target-column C) [--feature-set TI] --out DIR
PYTHONPATH=src python tools/run_experiment.py exp-d --spec study.json --out DIR
PYTHONPATH=src python tools/run_experiment.py size-dependence (--mesh ID.. | --manifest F | --smiles-csv F) --out DIR
```

Shared options: `--families` (inner-loop candidates; `poly2_ridge` names the degree-2 ridge
family below), `--seeds`, `--outer-splits`, `--inner-splits`, `--alpha`, `--n-jobs` (outer
folds in parallel; results identical), `--plan-frozen`, `--cache-dir`, `--offline`,
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
| `summary.json` | `counts` per cohort (manifest rows included, distinct parent ChEMBL IDs, standardised, rejections by code, connectivity groups); `jaccard`: `parent_chembl_id`, `structures_used` (the ChEMBL `standard_inchi_key` of the structures each retrieval returns: naive salts vs controlled parents; COHORT_DEFINITION §5's InChIKey column), `standardised_parent_inchikey` (after v1.0 standardisation); `intersection_standardised`; `controlled_in_naive_share` |
| `membership.csv` | one row per standardised parent InChIKey: naive id, controlled id, membership flags |
| `distributions.csv` | per variable (`n_heavy`, `n_bonds`, `n_rings`, each target): n, mean, SD, median, quartiles, range per cohort; two-sample KS statistic and p (`scipy.stats.ks_2samp`), BH `ks_q` across the variables. The two samples overlap, so the p values are descriptive |
| `correlations.csv` | per cohort × target × index: Spearman ρ, Fisher-z CI, p, n, BH `q` across the indices of that cohort × target, reason if undefined |
| `top_overlap.csv` | per target: the top-k indices by \|ρ\| in each cohort (`all`: every defined ρ; `bh_significant`: only q ≤ α, so possibly fewer than k), the overlap count and its Jaccard. Ties in \|ρ\| (many indices are rank-identical on small cohorts) are broken by registry order |
| `model.csv` (only `model=True`) | per cohort × target: nested `compare_feature_sets({size, TI})` with identical folds and connectivity groups: R² of each set, ΔR²(TI − size) with the corrected-t CI and p, `underpowered` |
| `molecules_naive.csv`, `molecules_controlled.csv` | the standardisation tables |
| `manifest_naive.csv`, `manifest_controlled.csv` (CLI) | both cohort manifests (ChEMBL data, CC BY-SA 3.0) |

Targets come from any target provider (`--provider`, `--targets`; default the RDKit computed
targets). Descriptive steps (counts, overlap, distributions, correlations) run for any disease;
the modelling step is quarantined (§0).

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

| output | content |
|---|---|
| `chain.csv` | the five R² values with CIs where defined |
| `optimism.csv` | B_train, B_CV, B_scaffold, CI, pooled variant, CI method |
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

Unknown keys raise. Drug names are **not** resolved automatically: each compound needs a
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
   index, reported R², difference). Undefined when n < degree + 2 or fewer than degree + 1
   distinct index values.
3. **LOOCV** of the same fit (`loo_q2` = 1 − PRESS/SST, `loo_rmse`, `loo_q2_f3`, `loo_ccc`).
4. **Nested CV** (`nested.csv`), outer `seeds` × `outer_splits`, connectivity groups kept
   together, no applicability domain (n is genre-sized):
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
   (`controls.csv`; p = (1 + #{null ≥ observed})/(1 + B), Phipson & Smyth 2010):
   * `pseudo_*`: best R² among k i.i.d. N(0, 1) pseudo-descriptors, `n_random` draws
     (descriptor-pool randomisation, Rücker et al. 2007);
   * `registry_*`: best R² among k indices drawn without replacement from the full TopoQSPR
     registry (indices defined and non-constant on these compounds), `n_random` draws. This is
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

PILOT_RUNTIME_PLACEHOLDER

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
