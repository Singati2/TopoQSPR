# Publication readiness (2026-10-13)

**Status: not ready to submit.** Blocking items are wording and provenance (cheap), the
deployment and the public release (need accounts and authorisation). No new experiment is
required for a defensible paper.

## Claim scope the evidence supports

1. **Software:** a versioned disease → MeSH → ChEMBL approved-drug manifest; 104 indices verified
   against independent sources; leakage-safe repeated nested CV with mandatory null, size and
   composition baselines; byte-identical rerun from run records.
2. **Frozen confirmatory result** (102 indices, random grouped folds): indices beat size for logP
   only after BH (+0.15 [0.08, 0.23]). With the registry held fixed, this shrinks to +0.10
   [−0.05, +0.25] on scaffold folds.
3. **Exploratory amendment E2** (locked before running): no detectable gain of the indices beyond
   an 11-count element composition on any endpoint (similarity-cluster split, BH q ≥ 0.45). The
   same holds on random and scaffold splits. Equivalence within ±0.05 R² is not shown. ECFP and RDKit
   2D are stronger baselines (RDKit 2D significantly better than TI+composition for MP, logP and logD).
4. **Methodological findings:** retrieval rule changes the drug set 3.6–16×. Unbounded descriptors
   (Hosoya Z, RDKit Ipc) make linear learners fail silently on a few large molecules. A fit-free
   asinh input transform removes the failure for composition-augmented sets but not for TI alone
   on logP. Published single-index studies reproduce numerically where graphs match; only
   one has the precision to bound its gain over size.
5. **Representation:** all index collisions in these data are graph isomorphisms (110 classes).
   They are real (ephedrine/pseudoephedrine, MP 34 vs 117.5 °C) but too rare to cap accuracy
   (R² ceilings ≤ 0.94–1.0, i.e. upper bounds well above the achieved 0.2–0.6).

## Must not claim

"Preregistered" · "1821 drugs with measured endpoints" · "indices add no information" or
"equivalent" · "the replicated studies are equivalent to size" · "the selected-model test is
validated by the bootstrap" · "the web interface is deployed" · disease-specific predictive
validity (no disease × endpoint cell has ≥ 50 drugs; maximum 47).

## Journal strategy

Recommended: **Journal of Cheminformatics, Research article**, evaluation-led ("what topological
indices add, measured properly"). Fallback: Digital Discovery (Full paper). SoftwareX only as a
cut-down tool paper with Experiments A–E moved to a methods paper. Reasons for and against are
in PUBLICATION_ROADMAP.md. The official guideline pages (word limits, APCs, article types) were
**not readable by the automated fetcher**; they must be checked in a browser. Acceptance
probability is not estimated.

## Unresolved data and licence items

* OPERA and Biogen data recorded as "MIT (repository LICENSE)". The dataset terms need checking (F14).
* Epirubicin MP 344.53 °C (Bradley, DrugBank-sourced) is suspect (F10). Keep it; flag it in the data card.
* ChEMBL CC BY-SA 3.0: manifests and result files that contain ChEMBL data must carry the licence.
* `related_papers/` holds copyrighted page images. They must leave the history before the repository is
  public. That needs a history rewrite, which needs explicit authorisation, and should happen
  once, before the release tag.
* ChemGraphX co-author overlap (M. Arockiaraj) must be disclosed to the editor.

## Reproducibility checklist

| item | state |
|---|---|
| Full test suite | 1111 passed, 3 skipped (local, `90ca304` tree) |
| Paper and supplement build | 34 + 14 pp, 0 errors, 0 undefined references (at `fecb800`; text unchanged since) |
| Archived results regenerate | yes, to 1e-9 from archived features (E1, E2 asserts) |
| Run-record commit IDs resolvable on GitHub | **no**: publish `sha_map_2026-10-12.tsv` or a Zenodo archive (F0) |
| Amended SVR control exact tree | **unrecoverable** (dirty, no diff). Package source unchanged by circumstantial evidence. A clean rerun (~6 h) was not done |
| Future dirty runs | now archive `code_diff.patch` (`90ca304`) |
| Permutation fallback | 0 rows in all archived targets (replayed); future runs withhold p if it fires |
| End-to-end raw → table reconstruction | not done in this audit (`topoqspr rerun`, ~5 h per study) |
| Deployed v1.0 web app | **no** (v0.x live) |
| Public repo, release tag, DOI | **no** |
| Audit branch pushed | **no** (local `audit-2026-10-12`) |
