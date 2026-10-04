# Paper format: following ChemGraphX (Advik, 2026-10-03)

Advik asked that the paper follow the format and writing style of ChemGraphX (Jacob, Clement,
Arockiaraj, Jindal, *J. Cheminform.* 18:39, 2026, PMC13020044) as closely as appropriate, without
copying text, and be written after the website works reliably and the code is verified.

## ChemGraphX structure (from the open-access full text)

Abstract with a "Scientific contribution" statement · Introduction · ChemGraphX (Description; Input
files; Work flow) · Comparative analysis · Validation and impact · Concluding remarks · Code metadata ·
Availability and requirements · Author contributions · Funding · Data availability · Code availability ·
Declarations (competing interests) · References · Supplementary information.
About 5,700 words of body text, 10 figures, 3 tables.

## Mapping for TopoQSPR

| ChemGraphX section | TopoQSPR content (current `paper/main.tex`) |
|---|---|
| Abstract + Scientific contribution | Abstract; contribution = reproducible disease-to-drug cohorts, verified indices, size/composition baselines, leakage-free evaluation, replication of published studies |
| Introduction | §1 Motivation and significance (genre growth, weaknesses, prior art on size and element blindness) |
| TopoQSPR: Description | §2 Architecture, run record |
| TopoQSPR: Input / disease resolution | §2.2 Disease resolution and drug-set curation (MeSH headings, manifest) |
| TopoQSPR: Work flow | CLI and web workflow figure (cohort preview → confirm → run → results) |
| Comparative analysis | Related-tools table (Mordred, ChemGraphX, DOPtools, QSPRpred, QSAR in the Browser, PaDEL, alvaDesc, PyTIE, RANQSAR) |
| Validation and impact | Index verification (closed forms, brute force, networkx, Mordred, ChemGraphX 41/41, printed values incl. errata, independent literature check); Experiments A–D results; impact |
| Concluding remarks | Conclusions and limitations |
| Code metadata; Availability and requirements | Code metadata table (C1–C9) and an "Availability and requirements" block (project name, home page, OS, language, requirements, licence) |
| Back matter | CRediT, funding, data and code availability (ChEMBL CC BY-SA, Zenodo DOI), declarations, generative-AI declaration |

## Venue implication

ChemGraphX's format and length are those of a *Journal of Cheminformatics* software article. SoftwareX
articles are much shorter (about 3,000–4,000 words, at most 6 figures; limit still unconfirmed because the
guide returns HTTP 403). Recommendation: **target J. Cheminform. (software article)**, which matches
Advik's requested format, the manuscript's current length and the study component; SoftwareX remains the
fallback after cutting. To decide with all co-authors.

## Order of work (Advik)

1. Website reliable on the v1.0 code (Phase 3, in progress).
2. Code accuracy and index definitions proof-checked against further independent papers (Phase 3).
3. Freeze the plan (CHECKPOINT decisions 3–5), run the experiments from a tagged release.
4. Write the paper in the structure above.
