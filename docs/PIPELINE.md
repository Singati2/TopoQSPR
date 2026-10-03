# TopoQSPR v1.0 pipeline, CLI and provenance

Code: `src/topoqspr/pipeline.py` (`run_study`, `rerun`), `src/topoqspr/provenance.py`
(run record), `src/topoqspr/target_providers.py` (target interface),
`src/topoqspr/cli.py` (command line). Tests: `tests/test_pipeline.py` (offline; the
only real data is the recorded pilot cohort in `tests/fixtures/pipeline/chembl`).
Roadmap item M6.

## 1. Stages

| # | Stage | What happens | Module |
|---|---|---|---|
| 1 | cohort | Controlled cohort from confirmed MeSH IDs (ChEMBL, cached per release), or a cohort manifest written earlier, or a SMILES CSV. Writes `manifest.csv` (+ `manifest.json`). | `cohort.build_cohort` |
| 2 | standardize | `graph.standardize` on every included structure. Rejections stay in `molecules.csv` with their code. A structure whose standardised parent has the same full InChIKey as an earlier one is excluded as `duplicate_inchikey:<kept id>` (it would otherwise sit in a training and a test fold at once). | `graph` |
| 3 | indices | All `indices.REGISTRY` entries (or `indices.names`). An undefined index is NaN with the registry's reason in `failures.csv`; an unexpected exception makes that molecule's indices NaN with the exception text. Never 0. | `indices` |
| 4 | features | `size`, `size_rings`, `composition`, `rdkit2d` (all of `descList`, subset per target below), `TI` and one `TI_<family>` set per index family. | `features` |
| 5 | targets | Each configured provider returns values (NaN + reason when missing). Rows with a missing value are dropped **for that target only**; the count is `n_target_missing` in `summary.csv`. | `target_providers` |
| 6 | evaluate | For every target and every learner: `evaluate.compare_feature_sets` on all feature sets with identical outer folds (EVALUATION_PROTOCOL.md §3–5). | `evaluate` |
| 7 | write | Result CSVs, `study.json`, then `run_record.json` last. | `provenance` |

**Learners.** `selected` is the protocol's default: family and hyperparameters are
chosen in the inner loop among all configured families. Every family name listed in
`evaluation.learners` is also run on its own (inner loop tunes only its
hyperparameters), which gives the baseline ladder per learner family. A family that
cannot be fitted for a set (OLS when p ≥ n_inner_train) produces failed folds and NaN
summaries with the reason, as specified in the protocol; nothing is dropped.

**rdkit2d per target.** The `rdkit2d` set of a target is `descList` minus
`features.TARGET_EXCLUSIONS[spec.rdkit2d_exclusion]`; the column lists actually used are
in `study.json` (`feature_sets.rdkit2d.<target>`).

**Groups.** `evaluation.split = "scaffold"` passes `evaluate.scaffold_groups` (Bemis–Murcko,
acyclic molecules as singletons) as `groups`, so outer and inner folds are scaffold-grouped.

## 2. Study configuration (JSON; YAML only if PyYAML is installed)

Every key is optional except `cohort`. Unknown keys raise `ValueError`, so a typo cannot
fall back to a default. `resolve_config` fills every default; the **resolved** config is
what is stored and hashed.

```jsonc
{
  "schema": 1,
  "name": "study",
  "cohort": {                       // one of three sources
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
  "output": {"features_format": "csv"},   // "parquet" additionally writes features.parquet (needs pyarrow)
  "runtime": {"n_jobs": 1, "chembl_cache_dir": null, "offline": false}
}
```

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

Experimental (class 3) providers are added by `register_target_provider(name, factory)`
or named in the config as `"package.module:factory"`; the factory is called as
`factory(names=<list or None>, **options)`. The pipeline checks that the returned columns
and index match the specs and that target names are unique across providers.

## 4. Outputs

All CSVs: UTF-8, `\n` line endings, no index column, floats in Python's shortest
round-trip form, rows in a fixed order (manifest order; then target, learner, set,
repeat, fold/position). Missing values are empty cells.

| file | content |
|---|---|
| `manifest.csv`, `manifest.json` | the cohort manifest and its sidecar (COHORT_DEFINITION.md §3); for a SMILES input, `id,smiles` and a sidecar with the input's SHA-256 |
| `molecules.csv` | one row per candidate: `id, input_smiles, status (ok/rejected), reason, detail, smiles, inchikey, parent_note, n_heavy, n_edges` |
| `features.csv` | one row per standardised molecule; columns prefixed by block: `size.*`, `composition.*`, `ti.*`, `rdkit2d.*` |
| `targets.csv` | `id` + one column per target |
| `failures.csv` | `block, row, feature, reason` for every NaN feature or target value |
| `folds.csv` | per (target, learner, set, repeat, fold): split sizes, chosen family and params, inner MSE, R², MAE, RMSE, AD coverage, status, reason (`NestedCVResult.folds`) |
| `predictions.csv` | out-of-fold predictions per (target, learner, set, repeat, row) with the AD flag |
| `paired.csv` | per-fold differences d_R², d_MAE, d_RMSE against the reference set |
| `inner.csv` | the inner-loop selection table of every outer fold |
| `summary.csv` | the baseline ladder per target × learner: `target, target_class, learner, set, reference, n, n_target_missing, p`, fold-mean and pooled R², MAE, failed/valid fold counts, `underpowered` (+ reason), `d_r2_mean/ci_low/ci_high/p`, `d_mae_*`, `d_r2_q` (BH across the sets of that target × learner), `family_counts`, `min_n_test` |
| `study.json` | feature-set column lists (rdkit2d per target), target specs, provider provenance |
| `run_record.json` | §5 |

The study-level BH correction over disease × target pairs (EXPERIMENT_PLAN §7) is applied
by the analysis that collects `d_r2_p` across runs, not by a single run.

## 5. Run record

`run_record.json` (schema 1) is written last and contains:

| key | content |
|---|---|
| `created_utc` | time stamp **passed in** by the caller (`--timestamp`, default: now) |
| `topoqspr_version` | `topoqspr.__version__` |
| `code` | `commit`, `dirty` (tracked files differ from HEAD), `untracked` (count), `describe`, `root` of the git tree the package was imported from; all `None` plus `reason` outside a checkout |
| `environment` | Python version/implementation/executable, platform, machine, `packages` (every installed distribution and version, `importlib.metadata`), and `env` (`PYTHONHASHSEED`, `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `OMP_PROC_BIND`) |
| `chembl` | `release`, `release_date` (from the client), source, offline flag; `null` for SMILES input |
| `config`, `config_sha256` | the full resolved config and the SHA-256 of the canonical JSON (sorted keys, `(",", ":")` separators, ASCII, NaN forbidden) of everything except `runtime` |
| `base_dir` | directory that relative input paths refer to |
| `seeds` | outer seeds, the inner-seed rule (`1000 * outer_seed + fold`), the estimator `random_state`, the outer splitter |
| `inputs` | SHA-256 of every input file (SMILES CSV / manifest + sidecar) |
| `manifest_sha256` | SHA-256 of the written `manifest.csv` |
| `counts` | candidates, standardised, rejections by code, non-missing values per target, feature-set sizes, number of failure records |
| `outputs` | SHA-256 and size of every file written before the record |
| `result_csvs` | the files a re-run must reproduce byte for byte |
| `timings_s` | wall time per stage, and per target × learner for `evaluate` |

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

**The Python version matters too.** `topoqspr indices` on `tests/fixtures/pipeline/tiny.csv`
in the Docker image (Python 3.12, same locked packages) and in the venv (Python 3.10.12)
differs in 6 of 3,564 cells, all Balaban J, by at most 2.0·10⁻¹⁶ relative (last digit
of the shortest float repr). Python 3.12's built-in `sum()` of floats is compensated
(Neumaier), 3.10's is not, and `indices._balaban_j` sums with `sum()`. Byte identity is
therefore guaranteed for the same Python minor version and package versions;
`environment.python_version` is in the record. Using `math.fsum` in the index sums would
remove this dependence (indices module; not changed here).

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
same CLI. Exit codes: 0 success, 1 re-run mismatch, 2 input or configuration error.

## 7. Packaging

* `requirements-lock.txt`: exact versions of the distributions `src/topoqspr` imports and
  their runtime requirements, from the project venv (Python 3.10.12); `xgboost` is pinned
  as `xgboost-cpu` (no CUDA wheels), so the package itself is installed with `--no-deps`.
  `requirements-lock-test.txt` adds pytest and the cross-check libraries.
* `Dockerfile`: `python:3.12-slim`, installs the lock, then the package; entry point
  `topoqspr`; BLAS/OpenMP threads set to 1. Built and checked locally on 2026-10-03
  (`topoqspr --version`, `indices` and an offline `cohort` from the fixture cache).
* CI (`.github/workflows/tests.yml`) installs `requirements-lock-test.txt`, then the
  package with `--no-deps`, and runs `pytest -q` (`pyproject.toml` adds
  `-m "not network"`).

## 8. Pilot smoke run (Migraine Disorders only)

SMOKE_RUN_PLACEHOLDER
