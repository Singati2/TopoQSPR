# TopoQSPR — roadmap to the paper

**Goal:** a SoftwareX software article describing TopoQSPR v1.0.
The manuscript lives in [`paper/`](../paper) (Overleaf-ready: `main.tex`, `refs.bib`).
Every code change below exists to produce a section, table or figure of that paper;
no number goes into `paper/main.tex` unless it was produced by a tagged release.

Target venue (provisional, being re-checked in the October 2026 audit): SoftwareX.

## Milestones

| # | Milestone | Feeds paper section | Status |
|---|---|---|---|
| M0 | Repository consolidated under `Singati2/TopoQSPR`, MIT license, README, paper draft | Code metadata (C1–C9) | done |
| M1 | Expert audit of code, paper and 2026 landscape (`docs/audit-2026-10.md`) | All; decides scope | in progress |
| M2 | **Drug-set curation**: disease → controlled MeSH/EFO term confirmed by the user; clinical-phase threshold (default approved); small molecules only; parent/salt normalisation; dedup by InChIKey; deterministic ordering; downloadable drug manifest with ChEMBL release | §2.2, Table 1, Fig. 2 | todo |
| M3 | **Verified indices**: unit tests against closed forms and literature values; fix or drop Hosoya Z; remove duplicated/ill-defined indices; cite the source of every family | §2.3, Table 2, Suppl. S1 | todo |
| M4 | **Honest evaluation**: true (repeated) nested CV; null and molecular-size baselines; ΔR² over size as headline; collapse duplicate targets; multiple-testing control for correlations | §2.4, Fig. 3 | todo |
| M5 | **Clean annotations**: remove invented Caco-2 formula and unit-mixed bioactivity; fix model count claim in UI | §2.5 | todo |
| M6 | **Official release engineering**: `pyproject.toml` package, CLI, pytest + GitHub Actions CI, pinned environment + Dockerfile, run-record (provenance) file, docs, v1.0.0 tag, Zenodo DOI, website renamed | Code metadata, §2.1 | todo |
| M7 | **Validation study**: fixed disease panel, approved drugs, current ChEMBL release; per-disease manifests and results archived | §3 Illustrative examples | todo |
| M8 | Write-up: results, impact, conclusions, figures; co-author review; submission | whole paper | todo |

## Rules for the paper

- Results come only from a tagged release; the run records and manifests are archived (Zenodo) and cited.
- Placeholders in `main.tex` are marked `\todo{}` and are never filled by hand.
- Authorship, author order and CRediT statement are agreed with all co-authors before submission.
