# Changelog

## Unreleased (v1.0.0, branch `v1.0-dev`)

### Added
- `src/topoqspr/` package: controlled disease-to-drug cohorts with versioned manifests
  (`cohort.py`, `chembl.py`), 102 verified topological indices (`graph.py`, `indices.py`),
  leakage-free repeated nested CV with a baseline ladder (`features.py`, `models.py`,
  `evaluate.py`), statistics (`stats.py`: corrected repeated k-fold t-test, BH and BY FDR,
  Spearman CIs, y-scrambling, applicability domain).
- Study pipeline and command line: `pipeline.py` (`run_study`, `rerun`; stages cohort,
  standardize, indices, features, targets, evaluate, protocol, write; atomic results
  directory), `cli.py` (the `topoqspr` entry point: `headings`, `cohort`, `indices`, `run`,
  `rerun`), `provenance.py` (run record: commit, dirty flag, package versions, config and
  manifest hashes, seeds; byte-identical re-run check), `target_providers.py` (target
  interface with an explicit target class and leakage exclusions; `rdkit` and `expdata`
  providers). See `docs/PIPELINE.md`.
- Measured (class-3) targets: `expdata.py` with licence-checked sources (AqSolDB, Bradley,
  OPERA, AstraZeneca ChEMBL deposition, Biogen), curation rules, full-InChIKey matching, and
  the pooled approved-drug cohort (`expdata.approved_drug_pool`). See `docs/EXPERIMENTAL_DATA.md`.
- Confirmatory protocol stage: applicability-domain report, Y-randomised control,
  whole-procedure paired bootstrap, confirmatory family table (BH/BY with fixed m,
  underpowered members excluded). See `docs/EVALUATION_PROTOCOL.md` §11.
- Experiment runners A (cohort sensitivity), C (optimism chain), D (replication harness) and
  size dependence: `experiments.py`, `tools/run_experiment.py`, with a quarantine guard
  until the plan is frozen. See `docs/EXPERIMENTS.md`.
- Reproducible environment: `requirements-lock.txt`, `requirements-lock-test.txt`,
  `Dockerfile` (commit and dirty flag written in at build time), GitHub Actions CI running
  the offline tests.
- Audit, verification and planning documents under `docs/`; reference, literature,
  cross-check and replication tools under `tools/`.

### Fixed (relative to the v0.x prototype on `main`)
- Hosoya Z (wrong for > 18 bonds), four wrong index numerators, duplicate indices, Balaban J for
  salts, failed indices filled with 0, non-nested CV, shared model instance across targets,
  non-deterministic drug selection, salt/parent duplicates, missing phase filter.

### Removed
- Invented Caco-2 formula and unit-mixed bioactivity columns (not part of the v1.0 package).
- Not ported from v0.x to the v1.0 package: SHAP explanations (attribution deferred,
  `docs/EXPERIMENT_PLAN.md` §9), PubChem/UniChem annotations and Lipinski/Veber/PAINS
  filters.

## 0.x (`main`)
- Original MolecuLens web prototype.
