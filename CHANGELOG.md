# Changelog

## Unreleased (v1.0.0, branch `v1.0-dev`)

### Added
- `src/topoqspr/` package: controlled disease-to-drug cohorts with versioned manifests
  (`cohort.py`, `chembl.py`), 102 verified topological indices (`graph.py`, `indices.py`),
  leakage-free repeated nested CV with a baseline ladder (`features.py`, `models.py`,
  `evaluate.py`), statistics (`stats.py`: corrected repeated k-fold t-test, BH and BY FDR,
  Spearman CIs, y-scrambling, applicability domain).
- Audit, verification and planning documents under `docs/`; reference and literature tools under `tools/`.

### Fixed (relative to the v0.x prototype on `main`)
- Hosoya Z (wrong for > 18 bonds), four wrong index numerators, duplicate indices, Balaban J for
  salts, failed indices filled with 0, non-nested CV, shared model instance across targets,
  non-deterministic drug selection, salt/parent duplicates, missing phase filter.

### Removed
- Invented Caco-2 formula and unit-mixed bioactivity columns (not part of the v1.0 package).

## 0.x (`main`)
- Original MolecuLens web prototype.
