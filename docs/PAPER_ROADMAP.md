# TopoQSPR — roadmap to the paper

**Goal:** a SoftwareX software article describing TopoQSPR v1.0.
The manuscript lives in [`paper/`](../paper) (Overleaf-ready: `main.tex`, `refs.bib`).
Every code change below exists to produce a section, table or figure of that paper;
no number goes into `paper/main.tex` unless it was produced by a tagged release.

Target venue: SoftwareX (Journal of Cheminformatics only if the replication result is strong).
External review and its verification: [`external-review-2026-10.md`](external-review-2026-10.md).

## Milestones

| # | Milestone | Feeds paper section | Status |
|---|---|---|---|
| M0 | Repository consolidated under `Singati2/TopoQSPR`, MIT license, README, paper draft | Code metadata (C1–C9) | done |
| M1 | Expert audit of code, paper and 2026 landscape ([`AUDIT_TOPOQSPR.md`](AUDIT_TOPOQSPR.md), [`LITERATURE_GAP_2026.md`](LITERATURE_GAP_2026.md)) | All; decides scope | done |
| M2 | **Drug-set curation**: disease → controlled MeSH/EFO term confirmed by the user; clinical-phase threshold (default approved); small molecules only; parent/salt normalisation; dedup by InChIKey; deterministic ordering; downloadable drug manifest with ChEMBL release | §2.2, Table 1, Fig. 2 | module done (`src/topoqspr/cohort.py`, `docs/COHORT_DEFINITION.md`); default evidence filter still to decide |
| M3 | **Verified indices**: unit tests against closed forms and literature values; cross-check against ChemGraphX; fix or drop Hosoya Z; fix DS1≡DS2 duplicate; reconcile reverse-degree definition; failures recorded as missing (never 0); cite the source of every family | §2.3, Table 2, Suppl. S1 | module done (`indices.py`, 102 verified indices, `docs/INDEX_VERIFICATION.md`); ChemGraphX cross-check done (41/41 shared indices agree); hybrid indices resolved from the primary paper; open: original v0.x SS_/reverse scripts from the Arockiaraj group (optional), source of BMA/TMA |
| M4 | **Honest evaluation**: true (repeated) nested CV with model-family selection inside the outer fold; null, size and RDKit-2D baselines; ΔR²(TI − size) with CI as headline; collapse duplicate targets; Benjamini–Hochberg FDR for correlations | §2.4, Fig. 3 | module done (`evaluate.py`, `stats.py`); open: full-family smoke run on real data before freezing |
| M5 | **Clean annotations**: remove invented Caco-2 formula and unit-mixed bioactivity; fix model count claim in UI | §2.5 | done in v1.0 package (no Caco-2 heuristic, no bioactivity); legacy web app still shows them |
| M6 | **Official release engineering**: `pyproject.toml` package, CLI, pytest + GitHub Actions CI, pinned environment + Dockerfile, run-record (provenance) file, docs, v1.0.0 tag, Zenodo DOI, website renamed | Code metadata, §2.1 | package, CLI (`topoqspr`), run records, lock file, Dockerfile, CI done (`docs/PIPELINE.md`); open: tag v1.0.0, Zenodo DOI, docs site, web demo hardening |
| M6b | **GO / NO-GO checkpoint** (`docs/CHECKPOINT.md`, per RECOMMENDATIONS.md); pilot on quarantined disease; freeze experiment plan | gate for M7 | pilot smoke run done (Migraine, n = 33); `docs/CHECKPOINT.md` drafted: NO-GO pending author decisions and experiment runners (Phase 2 building) |
| M7 | **Validation study**: 8–12 diseases above a pre-set minimum approved-cohort size, current ChEMBL release; per-disease manifests and results archived | §3 Illustrative examples | runners implemented (`experiments.py`, protocol stage); blocked on freeze (CHECKPOINT decisions) |
| M7b | **Replication**: re-analyse 3–5 published disease-specific topological-index QSPR studies (chosen with all co-authors): published in-sample R² vs nested-CV R² vs size baseline | §3 | harness implemented; 10 studies extracted and verified (`REPLICATION_DATA.md`); co-authors to confirm rows 1-5 |
| M8 | Write-up: results, impact, conclusions, figures; co-author review; submission | whole paper | todo |

## Rules for the paper

- Results come only from a tagged release; the run records and manifests are archived (Zenodo) and cited.
- Placeholders in `main.tex` are marked `\todo{}` and are never filled by hand.
- Authorship, author order and CRediT statement are agreed with all co-authors before submission.
