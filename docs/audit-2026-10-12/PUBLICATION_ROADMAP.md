# Publication roadmap (audit 2026-10-12)

Journal facts come from LITERATURE_2026.md §4. The official guideline pages of all three
journals returned 403 or a bot challenge to the automated fetcher, so the word limits and
APCs below are **secondhand** and must be checked in a browser before choosing.

## What the paper can claim today (survives this audit)

1. Software: a disease → MeSH → ChEMBL approved-drug manifest with explicit, versioned rules.
   No tool doing this was found (search scope in LITERATURE_2026.md). 104 indices with
   multi-source verification. Leakage-safe repeated nested CV with mandatory size and
   composition baselines. Run records with byte-identical rerun.
2. Empirical: on 351–711 approved drugs per measured endpoint, the indices beat a 2-feature
   size baseline only for logP after multiplicity correction (+0.15 [0.08, 0.23], random
   split, 102 indices). That gain shrinks on scaffold folds and its interval includes zero. In
   no endpoint do the indices reach the 11-count composition baseline. RDKit 2D descriptors
   are best throughout.
3. Methodological: the free-text retrieval rule changes the drug set 4–16×. Unbounded
   descriptors (Hosoya Z, RDKit Ipc) make linear learners fail catastrophically and silently
   on a few large molecules. Published single-index studies reproduce numerically where the
   graphs match, but their samples are too small to bound the gain over size, except Godlin's.

## What it cannot claim

* "Preregistered" (F4). "1821 drugs with measured endpoints" (F3). That the rdkit2d failure
  is due to Z (F1). That log P's scaffold gain is "the same size" (F5). That the replicated
  studies are *equivalent* to size (F6). That the web interface is deployed (F7).
* Disease-specific *predictive* validity. The modelling is on the pooled set; no disease ×
  endpoint cell has ≥ 50 measured drugs. Disease conditioning is a retrieval feature, not a
  validated predictor.

## Route A (minimum acceptable): SoftwareX

* Fit: tool paper. The official March 2026 template needs a public GitHub repo (forked by the
  journal), README.md, and a licence file named `Licence.txt`.
* Case **against**: secondhand word limit 3,000–4,000. The current body is about 4,700 words
  plus 34 pages with tables, so Experiments A–D would have to go. The scientific findings (the
  main value) would be lost or moved to a second paper. It is also a weaker venue for a
  methodological claim.
* Work: 2–3 weeks (cut, fix F0–F8, public repo, Zenodo DOI, deploy).

## Route B (recommended): Journal of Cheminformatics, Research/Methodology article

* Fit: evaluation-led paper ("what do topological indices add, measured properly"), with the
  software as the vehicle. J Cheminform Software articles must show a "significant advance
  over previously published software" by direct comparison. A Research article framed on
  the evaluation sidesteps that test. The 3-sentence Scientific Contribution statement is
  mandatory, and the draft already has one.
* Case **against**: ChemGraphX (same journal, Feb 2026) shares co-author Arockiaraj. Editors
  may see salami-slicing or overlap unless the cover letter is explicit. Reviewers will ask for
  an ECFP baseline and labelled-graph indices. APC is about USD 2,390 (secondhand).
* Required before submission:
  1. Fix F0–F8 (wording and provenance), keeping frozen numbers.
  2. Add E1 (composition-conditional, asinh and identical-registry scaffold; done in this audit
     as exploratory) as a dated, labelled sensitivity section.
  3. Add E4 (ECFP + ridge/extra-trees), declared as an amendment *before* running.
  4. Cite and contrast Ash et al. 2025 (JCIM 65:9398): they recommend repeated-measures
     ANOVA + Tukey HSD; TopoQSPR uses a corrected t-test + BH. Say why (pre-specified, two-set
     paired contrasts) in one sentence.
  5. Discuss the unlabeled-graph limitation explicitly, with the in-data collisions of F9
     (ephedrine/pseudoephedrine and others) as the illustration.
  6. Public repo, release tag, Zenodo DOI, deployed v1.0 site, and `related_papers/` purged
     *before* the tag (one history rewrite only).
* Time: 4–6 weeks.

## Route C (fallback): Digital Discovery, Full paper

* Fit: no page limit, a data reviewer checks code and data, and DOPtools is a precedent. APC is
  £2,200 (secondhand).
* Case **against**: an RSC readership oriented to ML/automation may find unlabeled TIs dated.
  The negative headline needs the E3 isomer study or the ECFP comparison to read as a
  contribution rather than a null.

## Ordered implementation plan

1. (Days) Manuscript corrections F1, F3–F6. Cite the new commit ID and the SHA map (F0). Fix the README (F8).
2. (Days) Disclose the ChemGraphX overlap. Decide the venue after checking the guideline pages by hand.
3. (1 week) Amendment entry for E1 (exploratory, already run) and E4 (declared before run). Run E4.
4. (1 week) Data card E6, including the F10 suspect value and the licence check F14.
5. (Account needed) Deploy v1.0 (`render.yaml`). Make the v0.x `/health` honest.
6. (Once, at release) Purge `related_papers/`, make the repo public, tag `v1.0.0`, Zenodo DOI. Fill `\todo`s.
7. (Optional, second paper) E3, the isomer-controlled study plus the Proposition 1 floor; Lean only if a toolchain is set up.

## What not to do

* Do not re-run the frozen analysis with transformed features and present it as the headline.
* Do not register on OSF now and call it preregistration.
* Do not add GNN/foundation-model baselines to the main paper. The question is TI vs size and
  composition (plus ECFP as the conventional floor).
* Do not drop the extrapolating molecules from headline metrics.
* Do not rewrite git history again after the release tag.
