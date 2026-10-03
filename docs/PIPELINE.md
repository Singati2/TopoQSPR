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
| 2 | standardize | `graph.standardize` on every included structure. Rejections stay in `molecules.csv` with their code. A structure whose standardised parent has the same full InChIKey as an earlier one is excluded as `duplicate_inchikey:<kept id>` (it would otherwise sit in a training and a test fold at once). Parents that share only the InChIKey **first block** (stereoisomers, racemate and enantiomer, isotopologues) are all kept but get a common `connectivity_group` (see *Connectivity groups* below). | `graph` |
| 3 | indices | All `indices.REGISTRY` entries (or `indices.names`). An undefined index is NaN with the registry's reason in `failures.csv`; an unexpected exception makes that molecule's indices NaN with the exception text. Never 0. | `indices` |
| 4 | features | `size`, `size_rings`, `composition`, `rdkit2d` (all of `descList`, subset per target below), `TI` and one `TI_<family>` set per index family. | `features` |
| 5 | targets | Each configured provider returns values (NaN + reason when missing). Rows with a missing value are dropped **for that target only**; the count is `n_target_missing` in `summary.csv`. | `target_providers` |
| 6 | evaluate | For every target and every learner: `evaluate.compare_feature_sets` on all feature sets with identical outer folds (EVALUATION_PROTOCOL.md §3–5). | `evaluate` |
| 7 | write | Result CSVs, `study.json`, then `run_record.json` last. Everything is written to a staging directory `.<out>.partial` next to `--out` and moved into `--out` only when the run has succeeded: a failed run leaves no partial results directory, and with `--overwrite` the old results are replaced only by a finished run. | `provenance` |

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
  "output": {"features_format": "csv"},   // "parquet" additionally writes features.parquet
  "runtime": {"n_jobs": 1, "chembl_cache_dir": null, "offline": false}
}
```

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
| `code` | `commit`, `dirty` (tracked files differ from HEAD), `untracked` (count), `describe`, `root` of the git tree the package was imported from, and `source = "git"`. An installed copy outside any git tree (the Docker image) uses the commit and dirty flag written into `topoqspr/_build_info.py` at build time (`source = "build_info"`, §7). Otherwise all `None` plus `reason` |
| `environment` | Python version/implementation/executable, platform, machine, `packages` (every installed distribution and version, `importlib.metadata`), and `env` (`PYTHONHASHSEED`, `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `OMP_PROC_BIND`) |
| `chembl` | `release`, `release_date` (from the client), source, offline flag, `data_licence`, and for a manifest input `manifest_verified`; `null` for SMILES input |
| `config`, `config_sha256` | the full resolved config and the SHA-256 of the canonical JSON (sorted keys, `(",", ":")` separators, ASCII, NaN forbidden) of everything except `runtime` |
| `base_dir` | directory that relative input paths refer to |
| `seeds` | outer seeds, the inner-seed rule (`1000 * outer_seed + fold`), the estimator `random_state`, the outer splitter actually used (§1 *Connectivity groups*) |
| `inputs` | SHA-256 of every input file (SMILES CSV / manifest + sidecar) |
| `manifest_sha256` | SHA-256 of the written `manifest.csv` |
| `counts` | candidates, standardised, rejections by code, non-missing values per target, feature-set sizes, number of failure records, `connectivity_groups` (groups and molecules) |
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
same CLI.

| exit code | meaning |
|---|---|
| 0 | success |
| 1 | `rerun`: a result CSV differs from the record (only this) |
| 2 | input or configuration error: invalid or non-mapping config, unknown key, missing file, non-empty `--out`, offline cache miss, release mismatch, missing optional dependency (parquet engine), target provider that cannot be imported or breaks its contract |
| 3 | unexpected internal error (a bug); the traceback is printed to stderr |

A run that fails leaves no results directory behind (§1, stage 7), so it can be repeated
without `--overwrite`.

## 7. Packaging

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
