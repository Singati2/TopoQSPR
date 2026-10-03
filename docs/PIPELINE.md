# TopoQSPR v1.0 pipeline, CLI and provenance

Code: `src/topoqspr/pipeline.py` (`run_study`, `rerun`), `src/topoqspr/provenance.py`
(run record), `src/topoqspr/target_providers.py` (target interface),
`src/topoqspr/cli.py` (command line). Tests: `tests/test_pipeline.py` (offline; the
only real data is the recorded pilot cohort in `tests/fixtures/pipeline/chembl`).
Roadmap item M6.

## 1. Stages

| # | Stage | What happens | Module |
|---|---|---|---|
| 1 | cohort | Controlled cohort from confirmed MeSH IDs (ChEMBL, cached per release), or the pooled approved-drug set (`expdata.approved_drug_pool`), or a cohort manifest written earlier, or a SMILES CSV. Writes `manifest.csv` (+ `manifest.json`). | `cohort.build_cohort`, `expdata.approved_drug_pool` |
| 2 | standardize | `graph.standardize` on every included structure. Rejections stay in `molecules.csv` with their code. A structure whose standardised parent has the same full InChIKey as an earlier one is excluded as `duplicate_inchikey:<kept id>` (it would otherwise sit in a training and a test fold at once). Parents that share only the InChIKey **first block** (stereoisomers, racemate and enantiomer, isotopologues) are all kept but get a common `connectivity_group` (see *Connectivity groups* below). | `graph` |
| 3 | indices | All `indices.REGISTRY` entries (or `indices.names`). An undefined index is NaN with the registry's reason in `failures.csv`; an unexpected exception makes that molecule's indices NaN with the exception text. Never 0. | `indices` |
| 4 | features | `size`, `size_rings`, `composition`, `rdkit2d` (all of `descList`, subset per target below), `TI` and one `TI_<family>` set per index family. | `features` |
| 5 | targets | Each configured provider returns values (NaN + reason when missing). Rows with a missing value are dropped **for that target only**; the count is `n_target_missing` in `summary.csv`. | `target_providers` |
| 6 | evaluate | For every target and every learner: `evaluate.compare_feature_sets` on all feature sets with identical outer folds (EVALUATION_PROTOCOL.md §3–5). | `evaluate` |
| 7 | protocol | AD coverage and in-domain metrics for every result (`ad.csv`, `ad_summary.csv`). With `config.protocol`: the Y-randomised control, the whole-procedure paired bootstrap and the confirmatory family table (§9; EVALUATION_PROTOCOL.md §11). | `pipeline` |
| 8 | write | Result CSVs, `study.json`, then `run_record.json` last. Everything is written to a staging directory `.<out>.partial` next to `--out` and moved into `--out` only when the run has succeeded: a failed run leaves no partial results directory, and with `--overwrite` the old results are replaced only by a finished run. | `provenance` |

**Learners.** `selected` is the protocol's default: family and hyperparameters are
chosen in the inner loop among all configured families. Every family name listed in
`evaluation.learners` is also run on its own (inner loop tunes only its
hyperparameters), which gives the baseline ladder per learner family. A family that
cannot be fitted for a set (OLS when n_inner_train < p + 2, i.e. p ≥ n_inner_train − 1) produces failed folds and NaN
summaries with the reason, as specified in the protocol; nothing is dropped.

**rdkit2d per target.** The `rdkit2d` set of a target is `descList` minus
`features.TARGET_EXCLUSIONS[spec.rdkit2d_exclusion]`; the column lists actually used are
in `study.json` (`feature_sets.rdkit2d.<target>`).

**Groups.** `evaluation.split = "scaffold"` passes `evaluate.scaffold_groups` (Bemis–Murcko,
acyclic molecules as singletons) as `groups`, so outer and inner folds are scaffold-grouped.

**Connectivity groups.** COHORT_DEFINITION.md de-duplicates on the full InChIKey because
stereoisomers are different compounds with possibly different measured properties, so
pairs such as citalopram/escitalopram, cetirizine/levocetirizine or omeprazole/esomeprazole
stay in the cohort. They have the same molecular graph, so every `size`, `composition` and
`ti.*` feature, all class 1 and 2 targets, and every `rdkit2d` column except
`NumUnspecifiedAtomStereoCenters` are identical within a pair. Splitting such a pair
between a training and a test fold would hand the learner its test answer. The pipeline
therefore keeps both members and labels each group of standardised parents that share the
InChIKey first block with the ID of its first member (`molecules.csv`,
`connectivity_group`; empty for molecules without such a partner), and **every group is
kept in one fold**, both outer and inner:

* `split = "random"` with at least one connectivity group: `balanced_group_kfold` with one
  group per connectivity group and every other molecule as a singleton. Without any
  connectivity group the split stays `KFold(shuffle=True)`, so cohorts without such pairs
  (the pilot) are unchanged;
* `split = "scaffold"`: members take the scaffold group of their first member (they share a
  scaffold anyway, except acyclic molecules, which would otherwise be singletons).

The run record gives the number of groups and molecules (`counts.connectivity_groups`) and
the splitter actually used (`seeds.outer_splitter`). Whether to keep one member per group
instead (for class 1 and 2 targets) is left to the frozen plan; grouping is the safe default
that keeps every compound and never leaks.

## 2. Study configuration (JSON; YAML only if PyYAML is installed)

Every key is optional except `cohort`. Unknown keys raise `ValueError`, so a typo cannot
fall back to a default. `resolve_config` fills every default; the **resolved** config is
what is stored and hashed.

```jsonc
{
  "schema": 1,
  "name": "study",
  "cohort": {                       // one of four sources
    "source": "chembl",             // controlled cohort (COHORT_DEFINITION.md)
    "mesh_ids": ["D008881"],        // confirmed descriptor IDs (required)
    "min_phase": 4.0,               // null = no phase rule
    "evidence_sources": null,       // e.g. ["FDA", "EMA", "ATC"]
    "exclude_sole_evidence": null,  // e.g. ["DailyMed"]
    "exclude_diagnostics": false,
    "include_descendants": false,   // not supported by the pipeline (list IDs explicitly)
    "disease_label": null,
    "chembl_release": null          // e.g. "ChEMBL_37": fail if another release is served
  },
  // or {"source": "pooled", "chembl_release": null}   pooled approved drugs, any indication
  //    (expdata.approved_drug_pool: the union of the phase-4 cohorts of every MeSH heading)
  // or {"source": "manifest", "path": "cohort/manifest.csv"}          (decision == include rows)
  // or {"source": "smiles_csv", "path": "x.csv", "id_column": "id", "smiles_column": "smiles"}
  //    relative paths are relative to the config file's directory (CLI) / base_dir (API)
  "standardize": {"min_heavy_atoms": 3},
  "indices": {"names": null},       // null = the full registry
  "feature_sets": ["null", "size", "size_rings", "composition", "TI", "TI_family", "rdkit2d"],
  "targets": [{"provider": "rdkit", "names": ["MolWt", "MolLogP", "TPSA", "MolMR"], "options": {}}],
  "evaluation": {
    "learners": null,               // null = ["selected", *families]
    "families": ["ridge", "elasticnet", "pls", "linear", "knn", "svr", "rf", "extratrees", "gbr", "xgb"],
    "outer_splits": 5, "seeds": [0, 1, 2, 3, 4], "inner_splits": 3,
    "random_state": 0, "ad": "knn", "split": "random",
    "reference": "size", "alpha": 0.05, "min_n": 50, "min_n_test": 10
  },
  "protocol": null,                 // null = no protocol stage; see §9 for the section
  "output": {"features_format": "csv"},   // "parquet" additionally writes features.parquet
  "runtime": {"n_jobs": 1, "chembl_cache_dir": null, "offline": false,
              "expdata_cache_dir": null}   // null = <TOPOQSPR_CACHE or ~/.cache/topoqspr>/expdata
}
```

**Baselines.** `resolve_config` requires `null`, `size` and `composition` in
`feature_sets` and `evaluation.reference = "size"`. A config may opt out only with
`evaluation.allow_custom_baselines: true`, which is kept in the resolved config (and so
changes its hash) and recorded under `baselines` in the run record; compliant configs
hash exactly as before.

`protocol: null` (the default) leaves the key out of the resolved config, so configs
written before the protocol stage existed keep their `config_sha256` (the golden hash of
`examples/pilot_migraine.json` is unchanged) and old run records still re-run.

**Validation happens before any work.** Besides unknown keys, `resolve_config` rejects an
empty or duplicated `evaluation.seeds` (identical repeats would be counted as independent
and narrow every CI and p value), `outer_splits` or `inner_splits` below 2, empty
`families`/`learners`, an empty `targets` list or an empty or duplicated `names` list
(`null` means all of a provider's targets), and `features_format: "parquet"` when neither
pyarrow nor fastparquet is installed (neither is in `requirements-lock.txt`). Every target
provider is imported and instantiated, and target names are checked for collisions, before
the cohort stage, so these errors cost seconds, not an evaluation.

`runtime` holds settings that cannot change a result (outer folds run in parallel give
identical frames; the cache location and offline mode do not change what is read). It is
recorded but **left out of the config hash**, so the same study on two machines has the
same `config_sha256`.

`examples/pilot_migraine.json` is the quarantined pilot (EXPERIMENT_PLAN §2): Migraine
Disorders, MeSH **D008881** (looked up in the ChEMBL_37 heading index:
`topoqspr headings migraine` lists D008881 *Migraine Disorders*, 143 indication rows, plus
D020325 *Migraine with Aura* and D020326 *Migraine without Aura*, which are not used),
pinned to ChEMBL_37, all ten families plus `selected`.

## 3. Targets and the provider interface

`target_providers.TargetProvider` (a `typing.Protocol`):

| member | contract |
|---|---|
| `name` | provider name (recorded) |
| `specs() -> list[TargetSpec]` | one spec per target column: `name`, `target_class` (EXPERIMENT_PLAN §4: 1 structural, 2 computed physicochemical, 3 external experimental), `rdkit2d_exclusion` (key of `features.TARGET_EXCLUSIONS`, or `None` = exclude nothing; no default, so every target has an explicit leakage decision), `units`, `source`, `description` |
| `compute(molecules, mols) -> TargetValues` | `molecules`: DataFrame indexed by molecule ID with `input_smiles`, `smiles` (standardised parent) and `inchikey` (full standard InChIKey of the standardised parent); `mols`: the RDKit parents in the same order. Returns `values` (same index, one float column per spec, NaN = missing) and `failures` (`{"row", "feature", "reason"}` for each NaN) |
| `provenance() -> dict` | what the run record must say: source, version, licence, file hashes, matching rule (full key or first block) |

The shipped provider `rdkit` computes, on `Chem.RemoveHs` of the standardised parent:

| target | class | rdkit2d exclusion | units |
|---|---|---|---|
| `MolWt` | 1 (sanity check only) | `MolWt` | g/mol |
| `MolLogP` (Wildman–Crippen) | 2 | `MolLogP` | log units |
| `TPSA` (Ertl) | 2 | `TPSA` | Å² |
| `MolMR` (Wildman–Crippen) | 2 | `MolMR` | cm³/mol |

The shipped provider `expdata` (`target_providers.ExpdataTargetProvider`) serves the
measured targets of `expdata.get_target` (EXPERIMENTAL_DATA.md): `mp_bradley`,
`mp_opera`, `logp_opera`, `logd74_az`, `logs_aqsoldb`, `logs_opera`, `logs74_az`,
`logs68_biogen` (`names: null` = all eight). Every target is class 3, has
`rdkit2d_exclusion: null` (a measured value is computed from no RDKit descriptor; by the
leakage rule of EVALUATION_PROTOCOL.md §1 alternative estimators such as Crippen `MolLogP`
for measured logP are kept, even though Crippen's atom contributions were themselves fitted
to measured logP data) and carries the source in `source` (`expdata:<source key>`).

* **Join:** the full standard InChIKey of the standardised parent (the same
  `graph.standardize` that `expdata` applies to the source records). There is no
  first-block fallback, so a stereoisomer never inherits its partner's measurement.
* **Missing values** are NaN with one of the reasons `no_inchikey`, `no_measurement`,
  `no_measurement_stereo_mismatch` (only another stereoisomer or isotopologue has a value;
  not used) or `replicates_disagree` (replicate SD above the target's `disagreement_sd`;
  the provider loads with `drop_disagreeing=False` and applies the rule itself so that the
  reason is recorded).
* **Options** (hashed): `drop_disagreeing` (default true), `keep_salt_forms` (default
  false; sensitivity analysis), `chembl_release` (for the AstraZeneca targets; default:
  the release of the cohort's ChEMBL client, else ChEMBL_37).
* **Runtime** (not hashed): the pipeline calls `set_runtime(cache_dir=runtime.expdata_cache_dir,
  chembl_cache_dir=..., offline=..., client=<cohort client>)` on any provider that defines
  it, so the AstraZeneca (ChEMBL-hosted) targets come from the cohort's release and cache.
* **Provenance** (`study.json` and `run_record.json` → `target_providers`): matching rule,
  options, `expdata.PROCESSING_VERSION`, and per target the source key, title, citation,
  DOI, version, **licence**, redistribution flag, URL or ChEMBL assay, the raw-data
  SHA-256 and processing fingerprint, the number of source parents, full-key matches,
  parents dropped for disagreeing replicates, and values used.

A target with fewer molecules (or connectivity groups) than the CV needs is not evaluated;
it is listed with the reason in `run_record.json` → `counts.skipped_targets` (a pilot-size
cohort covers only a handful of drugs for some measured endpoints).

Experimental (class 3) providers are added by `register_target_provider(name, factory)`
or named in the config as `"package.module:factory"`; the factory is called as
`factory(names=<list or None>, **options)`. The pipeline checks that the returned columns
and index match the specs and that target names are unique across providers, and it
**enforces the NaN + reason rule**: a provider that returns a NaN without a matching
`{"row", "feature", "reason"}` failure record, or any ±inf, stops the run with a
`ValueError` naming the cells (before the evaluate stage).

## 4. Outputs

All CSVs: UTF-8, `\n` line endings, no index column, floats in Python's shortest
round-trip form, rows in a fixed order (manifest order; then target, learner, set,
repeat, fold/position). Missing values are empty cells.

**Reading the CSVs.** Only empty cells are missing values. A default `pandas.read_csv`
also turns the strings `null`, `NA`, `None`, `nan`, ... into NaN, which silently drops the
`null` baseline set (`set == "null"`) from every result table. Read with
`keep_default_na=False, na_values=[""]`, or use the helper
`topoqspr.pipeline.read_results(out_dir, "summary")`, which does exactly that.

**Licence.** When the cohort comes from ChEMBL, `manifest.csv`, `molecules.csv` and every
file derived from them contain ChEMBL data (CC BY-SA 3.0, EMBL-EBI). The sidecar
`manifest.json` and `run_record.json` (`chembl.data_licence`) name the licence;
redistributed results must keep the attribution and the licence (README, *Data*).

| file | content |
|---|---|
| `manifest.csv`, `manifest.json` | the cohort manifest and its sidecar (COHORT_DEFINITION.md §3, plus `data_licence`); for a SMILES input, `id,smiles` and a sidecar with the input's SHA-256. With `cohort.source = "manifest"` the input CSV must hash to its sidecar's `manifest_sha256`, otherwise the run stops (an edited manifest would otherwise be labelled with the ChEMBL release); without a sidecar the run proceeds and the record says `chembl.manifest_verified: false` |
| `molecules.csv` | one row per candidate: `id, input_smiles, status (ok/rejected), reason, detail, smiles, inchikey, parent_note, n_heavy, n_edges, connectivity_group` |
| `features.csv` | one row per standardised molecule; columns prefixed by block: `size.*`, `composition.*`, `ti.*`, `rdkit2d.*` |
| `targets.csv` | `id` + one column per target |
| `failures.csv` | `block, row, feature, reason` for every NaN feature or target value |
| `folds.csv` | per (target, learner, set, repeat, fold): split sizes, chosen family and params, inner MSE, R², MAE, RMSE, AD coverage, status, reason (`NestedCVResult.folds`) |
| `predictions.csv` | out-of-fold predictions per (target, learner, set, repeat, row) with the AD flag |
| `paired.csv` | per-fold differences d_R², d_MAE, d_RMSE against the reference set |
| `inner.csv` | the inner-loop selection table of every outer fold |
| `summary.csv` | the baseline ladder per target × learner: `target, target_class, learner, set, reference, n, n_target_missing, p`, fold-mean and pooled R², MAE, failed/valid fold counts, `underpowered` (+ reason), `d_r2_mean/ci_low/ci_high/p`, `d_mae_*`, `d_r2_q` (BH across the sets of that target × learner), `family_counts`, `min_n_test` |
| `ad.csv` | per (target, learner, set, repeat, fold): n_test, n_in_domain, n_outside, share_outside, R²/MAE of the in-domain and of the outside test rows, ad_reason (written whenever `evaluation.ad` is set) |
| `ad_summary.csv` | per (target, learner, set): mean and max share outside over folds, pooled share, and mean over repeats of the pooled R²/MAE on all, in-domain and outside rows; `reason` when the AD was not fitted (e.g. the zero-column `null` set) |
| `permutation.csv` | protocol: one row per (target, set, permutation): the permuted-target nested-CV statistics and the RNG used |
| `permutation_summary.csv` | protocol: per (target, set): observed mean pooled R², K, empirical p (+1), null mean/SD/95th percentile/max, valid and failed permutations, `uninformative` |
| `bootstrap.csv` | protocol: per (target, set, reference, resample): the RNG, rows and units drawn, d_R², d_MAE and both sets' fold-mean R²/MAE; `resample = -1` is the original data |
| `confirmatory.csv` | protocol: the confirmatory family table (EVALUATION_PROTOCOL.md §11) |
| `study.json` | feature-set column lists (rdkit2d per target), target specs, provider provenance |
| `run_record.json` | §5 |

`d_r2_q` in `summary.csv` corrects across the sets of one target × learner only. The
confirmatory BH/BY correction is in `confirmatory.csv`, over the family that the run's
`protocol.family` defines; a family that spans several runs (e.g. diseases run separately)
is corrected by the analysis that collects `d_r2_p` across runs.

## 5. Run record

`run_record.json` (schema 1) is written last and contains:

| key | content |
|---|---|
| `created_utc` | time stamp **passed in** by the caller (`--timestamp`, default: now) |
| `topoqspr_version` | `topoqspr.__version__` |
| `code` | `commit`, `dirty` (tracked files differ from HEAD), `untracked` (count), `describe`, `root` of the git tree the package was imported from, and `source = "git"`. An installed copy outside any git tree (the Docker image) uses the commit and dirty flag written into `topoqspr/_build_info.py` at build time (`source = "build_info"`, §7). Otherwise all `None` plus `reason` |
| `environment` | Python version/implementation/executable, platform, machine, `packages` (every installed distribution and version, `importlib.metadata`), and `env` (`PYTHONHASHSEED`, `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `OMP_PROC_BIND`) |
| `chembl` | `release`, `release_date` (from the client), source, offline flag, `data_licence`, and for a manifest input `manifest_verified`; `null` for SMILES input |
| `config`, `config_sha256` | the full resolved config and the SHA-256 of the canonical JSON (sorted keys, `(",", ":")` separators, ASCII, NaN forbidden) of everything except `runtime` |
| `base_dir` | directory that relative input paths refer to |
| `seeds` | outer seeds, the inner-seed rule (`1000 * outer_seed + fold`), the estimator `random_state`, the outer splitter actually used (§1 *Connectivity groups*) |
| `inputs` | SHA-256 of every input file (SMILES CSV / manifest + sidecar) |
| `manifest_sha256` | SHA-256 of the written `manifest.csv` |
| `counts` | candidates, standardised, rejections by code, non-missing values per target, feature-set sizes, number of failure records, `connectivity_groups` (groups and molecules), `skipped_targets` (target → reason) |
| `targets` | every target spec: name, class (EXPERIMENT_PLAN §4) and label, rdkit2d exclusion, units, source, description |
| `target_providers` | each provider's `provenance()` (for `expdata`: sources, licences, hashes, match counts; §3) |
| `protocol` | `null` without `config.protocol`; otherwise the expanded confirmatory **family** (target, set, reference), the family rule, the confirmatory learner and its families, the **practical threshold**, q, alpha, the multiplicity, **family-membership** (`family_membership`, `min_n`, `min_n_test`), informative and practical rules, and for `permutation` and `bootstrap` `enabled` (with the reason when off), K/B, seed, RNG rule, permutation scheme (row-wise or block, per target), learner or families (a fast bootstrap subset is labelled as such), statistic and interval rules |
| `outputs` | SHA-256 and size of every file written before the record |
| `result_csvs` | the files a re-run must reproduce byte for byte |
| `timings_s` | wall time per stage, per target × learner for `evaluate`, and per target (permutation) or family row (bootstrap) for `protocol` |

**Re-run.** `topoqspr rerun RUN_RECORD --out DIR` (`pipeline.rerun`) checks the record's
config hash, checks every input file against its recorded SHA-256, builds an **offline**
ChEMBL client pinned to the recorded release (`--cache-dir` points at the cache that holds
it; a missing entry raises `CacheMiss`, it is never fetched), runs the study and compares
each result CSV with the recorded SHA-256. It writes `rerun_report.json` (`identical`,
`different`, `missing`, package versions that differ from the record, both commits) and
exits 1 if any result CSV differs. The time stamp is not part of any result CSV, so a
re-run on another day is byte-identical; `manifest.json` and `run_record.json` differ only
in their time stamps and environment fields. Tested in `tests/test_pipeline.py`
(`test_rerun_byte_identical`, on the recorded pilot cohort, and the CLI round trip).

Determinism rests on: seeded splitters and estimators, `n_jobs = 1` inside every
estimator (outer folds may run in parallel; results are identical), sorted rows, sorted
JSON keys, no clock or hash-seed dependence in any CSV (the config hash is tested under
two `PYTHONHASHSEED` values). Different package versions (numpy, scikit-learn, xgboost,
RDKit) can change floating-point results; the record lists them, and the re-run report
names any that differ.

**Python version (fixed; historical note).** Before commit `069133a`, `topoqspr indices` on
`tests/fixtures/pipeline/tiny.csv` in the Docker image (Python 3.12, same locked packages)
and in the venv (Python 3.10.12) differed in 6 of 3,564 cells, all Balaban J, by at most
2.0·10⁻¹⁶ relative (last digit of the shortest float repr): Python 3.12's built-in `sum()`
of floats is compensated (Neumaier), 3.10's is not, and `indices._balaban_j` summed with
`sum()`. `_balaban_j` now sums with `math.fsum` (correctly rounded on every version).
Re-checked on 2026-10-03 at `v1.0-dev` (e35dc4a): an image built from that commit
(Python 3.12.15) and the venv (Python 3.10.12) write byte-identical `indices` output and
failure files for that fixture. Cross-version identity therefore now depends on the package
versions only; `environment.python_version` stays in the record.

## 6. Command line

```
topoqspr headings TEXT [--limit 10] [--cache-dir D] [--offline] [--release R]
topoqspr cohort --mesh D######  [--mesh ...] --out DIR [--min-phase 4 | --all-phases]
                [--evidence-sources T]... [--exclude-sole-evidence T]... [--exclude-diagnostics]
                [--label L] [--timestamp ISO] [--cache-dir D] [--offline] [--release R]
topoqspr indices --smiles-csv F --out G [--id-column id] [--smiles-column smiles]
topoqspr run --config study.json --out DIR [--timestamp ISO] [--n-jobs N] [--cache-dir D]
             [--offline] [--overwrite]
topoqspr rerun RUN_RECORD --out DIR [--timestamp ISO] [--cache-dir D] [--base-dir D]
               [--n-jobs N] [--online]
```

`indices` writes the standardisation status and every index per molecule (exact
integers, e.g. Hosoya Z beyond 2⁵³, are kept) and `<out stem>.failures.csv` with the
reason for every undefined value or rejected structure. `python -m topoqspr` is the
same CLI.

| exit code | meaning |
|---|---|
| 0 | success |
| 1 | `rerun`: a result CSV differs from the record (only this) |
| 2 | input or configuration error: invalid or non-mapping config, unknown key, missing file, non-empty `--out`, offline cache miss, release mismatch, missing optional dependency (parquet engine), target provider that cannot be imported or breaks its contract |
| 3 | unexpected internal error (a bug); the traceback is printed to stderr |

A run that fails leaves no results directory behind (§1, stage 8), so it can be repeated
without `--overwrite`.

## 7. Packaging

* `openpyxl` is needed to build the `mp_bradley` target from its raw Excel file
  (`expdata`; it is in the `expdata` extra of `pyproject.toml`). Check that the lock or
  image you install includes it before running a study with that target
  (EXPERIMENTAL_DATA.md §2).
* `requirements-lock.txt`: exact versions of the distributions `src/topoqspr` imports and
  their runtime requirements, from the project venv (Python 3.10.12); `xgboost` is pinned
  as `xgboost-cpu` (no CUDA wheels), so the package itself is installed with `--no-deps`.
  `requirements-lock-test.txt` adds pytest and the cross-check libraries.
* `Dockerfile`: `python:3.12-slim`, installs the lock, then the package; entry point
  `topoqspr`; BLAS/OpenMP threads set to 1. Built and checked locally on 2026-10-03
  (`topoqspr --version`, `indices` and an offline `cohort` from the fixture cache).
  In the image the package sits in site-packages, outside any git tree, so git cannot
  name the commit at run time (mounting the source tree does not help). The build
  therefore **requires** the commit and writes it into `topoqspr/_build_info.py`:

  ```bash
  docker build -t topoqspr \
      --build-arg GIT_COMMIT="$(git rev-parse HEAD)" \
      --build-arg GIT_DIRTY="$(git status --porcelain --untracked-files=no | grep -q . && echo 1 || echo 0)" .
  ```

  A build without these arguments fails. Run records from the image then carry
  `code.commit`, `code.dirty` and `code.source = "build_info"`; the image is also labelled
  `org.opencontainers.image.revision`. `_build_info.py` is git-ignored and never used inside
  a git checkout, where the live git state wins.
* CI (`.github/workflows/tests.yml`) installs `requirements-lock-test.txt`, then the
  package with `--no-deps`, and runs `pytest -q` (`pyproject.toml` adds
  `-m "not network"`).

## 8. Pilot smoke run (Migraine Disorders only)

Purpose: check that the full pipeline runs end to end on real data with every learner
family and measure its cost before the plan is frozen (EVALUATION_PROTOCOL §9). Only the
quarantined pilot disease was used; **no performance number was read or is reported
here**, and none may set a threshold, metric, target or panel (EXPERIMENT_PLAN §2).

Command (2026-10-03, mpcrlab12, 16 cores, load average 35–57 from other jobs; each of the
13 workers got about 25 % of a core):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  topoqspr -v run --config examples/pilot_migraine.json --out <scratch>/pilot_smoke --n-jobs 13
```

Run record: commit `ec7ae2c` (clean tree), Python 3.10.12, ChEMBL_37 (2026-05-01, live
API, default cache), `config_sha256` `f9478a2c…e2cdea`, `manifest_sha256` `7d3e104d…ccb62a`.
The results directory stays in scratch space and is not committed.

**Cohort.** Migraine Disorders (D008881): 152 manifest rows: 33 included parents,
30 salt forms collapsed onto their parents, 84 parents below phase 4, 5 not small
molecules. All 33 passed standardisation (no rejections, no InChIKey duplicates) and
all four targets were defined for all 33 (n = 33). No feature value was undefined
(`failures.csv` is empty; every parent has maximum degree ≥ 3, so the R₃ indices
exist). Every summary row is labelled `underpowered` (n = 33 < 50; smallest outer test
fold 6 < 10), as the protocol requires.

**Scale.** 4 targets × 11 learners (`selected` + 10 families) × 12 feature sets
(null, size, size_rings, composition, TI, 6 TI families, rdkit2d) × 25 outer folds =
13,200 outer folds.

**Wall time per stage**

| stage | seconds |
|---|---|
| cohort (live ChEMBL incl. the 61-page heading index; cached afterwards) | 32.9 |
| standardize | 2.0 |
| features (102 indices + rdkit2d + baselines) | 3.5 |
| targets | 0.2 |
| evaluate | 23,050.6 (6.4 h) |
| write | 8.9 |

**Evaluate, by learner** (seconds summed over the 4 targets, 12 sets × 25 folds each):

| selected | rf | extratrees | gbr | xgb | elasticnet | knn | ridge | svr | pls | linear |
|---|---|---|---|---|---|---|---|---|---|---|
| 11,637 | 4,083 | 2,905 | 2,076 | 961 | 717 | 183 | 160 | 141 | 130 | 40 |

By target: MolWt 7,281 s (ran first, under the heaviest load), MolLogP 5,538 s, TPSA 5,271 s,
MolMR 4,941 s. The four tree ensembles take 88 % of the single-family time. `selected`
repeats the inner fits of all ten families, so it costs about as much as the ten
single-family runs together (11,637 s vs 11,396 s). Taking `selected` from the
single-family inner tables (same splits and seeds, so the same choice) would halve
`evaluate`; not implemented.

**Failures (expected by the protocol, reported, not dropped).** 400 of 13,200 outer folds
failed, all with reason "no family could be fitted": the `linear` (OLS) learner on `TI`,
`TI_degree`, `TI_reverse_degree` and `rdkit2d`, where p ≥ n_inner_train − 1 (e.g. "needs
n_train >= 104 for p=102, got 17"), for all 25 folds and 4 targets. Those 16 summary rows
have NaN fold means and NaN CIs with the reason. In the `selected` learner OLS is skipped
on those sets in the inner loop (800 `skipped` rows in `inner.csv`) and another family is
chosen. No other fold failed or had an undefined metric; the AD was fitted in every fold
except the 1,100 folds of the zero-column `null` set ("no feature columns").

**Output sizes**: predictions.csv 6.8 MB, inner.csv 2.6 MB, folds.csv 2.3 MB,
paired.csv 1.1 MB, summary.csv 0.19 MB (528 rows), features.csv 0.12 MB.

**Budget implication** (cost only): at this load, one disease × 4 targets with all
learners took 6.4 h for n = 33. Tree-ensemble cost grows with n_train, and the protocol's
own measurement (§9) is about 340 s per outer fold at n_train = 123 for the ten families,
so main-panel cohorts (n ≥ 50) with all eleven learners will need a dedicated machine or
a reduced learner list, decided before the plan is frozen.

## 9. Protocol stage (confirmatory analysis)

Enabled by a `protocol` section in the study config (EXPERIMENT_PLAN §7; statistics in
EVALUATION_PROTOCOL.md §11). Defaults:

```jsonc
"protocol": {
  "family": null,                 // "auto" = every class-3 target x (TI vs evaluation.reference);
                                  // or [{"target": "mp_bradley", "set": "TI", "reference": "size"}, ...]
                                  // or the shorthand ["mp_bradley:TI-size", "logp_opera:TI_degree"]
  "learner": "selected",          // whose results form the table; must be in evaluation.learners
  "practical_threshold": 0.05,    // minimal meaningful Delta R2 (plan: to be confirmed)
  "q": 0.05,                      // FDR level for BH and BY
  "permutation": {"enabled": true, "n_permutations": 200, "seed": 0,
                  "learner": null,          // null = protocol.learner (must be in evaluation.learners)
                  "feature_sets": null,     // null = the sets and references of the family
                  "quantile": 0.95},
  "bootstrap": {"enabled": null, "n_boot": 200, "seed": 0,  // null = on for pooled cohorts only
                "families": null}           // null = protocol.learner's families; or a fast subset
}
```

**Bootstrap scope.** The plan (EXPERIMENT_PLAN §7) uses the whole-procedure bootstrap on
the **pooled** set only; at per-disease n the resampled test folds are too small for R² to
be defined (§10). `bootstrap.enabled` defaults to `null`: on for `cohort.source = "pooled"`, off for every
other cohort unless set to true explicitly; the rule is recorded in the run record.

Validation (before any work): the family's sets must be evaluated feature sets (a
`TI_<family>` set needs `TI_family`), the reference must be `evaluation.reference`,
entries must be unique, the learners must be in `evaluation.learners`, K ≥ 1, B ≥ 2,
bootstrap families must exist, and every family target must be supplied by a provider
(checked before the cohort stage).

Outputs: `confirmatory.csv`, `permutation.csv`, `permutation_summary.csv`,
`bootstrap.csv` (§4), and the run record's `protocol` entry (§5), which holds the expanded
family, the threshold and every seed and rule. All are in `result_csvs`, so `rerun`
checks them byte for byte.

**Cost.** Each permutation reruns one nested CV per set; each bootstrap resample reruns
two. With the plan's K = B = 200 and a family of four targets × (TI vs size), that is
4 × 2 × 200 + 4 × 2 × 201 = 3,208 nested CVs on top of the main analysis, so the main run
needs a fast learner list (CHECKPOINT decision 4) and, for the bootstrap, possibly a
fast family subset (`protocol.bootstrap.families`; labelled in the run record).
Permutations and resamples run in parallel with `runtime.n_jobs` (identical results).

## 10. Protocol smoke run (Migraine Disorders only)

Purpose: run the protocol stage end to end on real data, with measured targets, and
measure its cost. Only the quarantined pilot was used (EXPERIMENT_PLAN §2); **no
performance number, CI, p or q value was read or is reported here**. Only counts, status
fields and timings were inspected (with a script that prints nothing else).

Command (2026-10-03, mpcrlab12, load average 25–30 from other jobs):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  topoqspr -v run --config pilot_protocol_smoke.json --out <scratch>/pilot_smoke --n-jobs 4 --offline
```

Config (kept in scratch space, not committed): the D008881 cohort pinned to ChEMBL_37;
feature sets `null, size, composition, TI, rdkit2d`; targets `rdkit` (MolLogP, TPSA) and
`expdata` (mp_bradley, logp_opera, logd74_az, logs_aqsoldb); learner `selected` over the
fast subset ridge, PLS, kNN; 5 × 5 outer, 3 inner; `protocol.family = "auto"` (the four
measured targets × TI vs size), K = 19 permutations, B = 19 bootstrap resamples with the
bootstrap family subset `["ridge"]`.

Run record: commit `7a39c1c` (clean tree), Python 3.10.12, ChEMBL_37 (2026-05-01, offline
cache), `config_sha256` `9d84370a…128d24`, `manifest_sha256` `7d3e104d…ccb62a` (the same
cohort as §8).

**Counts.** 33 parents, all standardised. Measured values matched on the full InChIKey:
mp_bradley 15, logs_aqsoldb 12, logp_opera 10, logd74_az 8 (the coverage counts of
EXPERIMENTAL_DATA.md). Missing values carry their reasons: `no_measurement` 17–22 per
target, `no_measurement_stereo_mismatch` (only another stereoisomer has a value) 3 for
logd74_az and 1 for logp_opera, `replicates_disagree` mp_bradley 1, logp_opera 2,
logd74_az 0, logs_aqsoldb 1. Every confirmatory row is labelled underpowered (n < 50).
(This run predates the family-membership rule of EVALUATION_PROTOCOL §11.2: it still gave
these rows BH/BY q values; under the current code every one of them is left out of the
multiplicity family, `family_m` = 0, with NaN q and no informative flag.) One target (logd74_az, n = 8) has outer test folds with a single
row, so 10 of 25 fold R² are undefined and its corrected-t p is NaN with the reason, as
specified. The bootstrap intervals are NaN with "k of 19 resamples undefined" for all four
rows: at n = 8–15 a resampled test fold often holds copies of one or two molecules
(R² undefined), and the protocol refuses to compute an interval from the surviving
resamples. All 152 permutation runs had a defined statistic. `null` is the only set
without an AD (no columns).

**Wall time per stage** (seconds; n_jobs = 4 over outer folds in the main analysis and
over permutations / resamples in the protocol stage)

| stage | s |
|---|---|
| cohort (offline cache) | 0.4 |
| standardize | 1.7 |
| features | 3.5 |
| targets (`expdata.get_target` for 4 targets from the processed cache, plus the full-key join) | 170.2 |
| evaluate (6 targets × 5 sets × 25 folds, `selected` over 3 families) | 529.8 |
| protocol: AD report | 3.2 |
| protocol: permutations (4 targets × 2 sets × 19 = 152 nested CVs, `selected`) | 1,906.3 (415–550 per target) |
| protocol: bootstrap (4 rows × 20 runs × 2 sets = 160 nested CVs, ridge only) | 885.9 (178–246 per row) |
| write | 1.0 |

So at pilot size and this load one nested CV costs about 12.5 s wall for `selected`
over ridge/PLS/kNN and 5.5 s for ridge alone (4 parallel workers). **Budget implication**
(cost only): the plan's K = B = 200 is about 10.5 times this run, i.e. about 5.6 h for the
permutations and 2.5 h for the bootstrap per four-target family at n ≈ 10–15, before any
growth with n; on the pooled set (hundreds of drugs per target) with tree ensembles in the
learner list it needs a dedicated machine, a reduced learner list, or a smaller K/B
fixed in the frozen plan (CHECKPOINT decision 4). The bootstrap also needs a cohort large
enough that resampled outer test folds have defined R²; at pilot size it cannot produce an
interval by design.
