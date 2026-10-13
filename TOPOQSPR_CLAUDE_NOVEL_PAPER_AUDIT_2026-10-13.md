# CLAUDE CODE MASTER PROMPT — TOPOQSPR
## Adversarial verification, mathematical novelty, chemical validity, reproducible experiments and paper decision
**Prepared:** 2026-10-13  
**Repository:** https://github.com/Singati2/TopoQSPR  
**Status:** Instructions to investigate; NOT proof of novelty or a request to manufacture a positive result.

### Your role and objective

Act as an independent team of specialists: chemical graph theorist, cheminformatician, statistician, rigorous software auditor, research-reproducibility engineer, and skeptical journal referee. Challenge both the repository AND this prompt. Your goal is **the strongest scientifically defensible and potentially publishable original paper**, whether the result is positive, mixed or negative.

Do **not** accept that the current project is novel, that any audit finding is correct, or that a positive result is obtainable. Do **not** claim that experiments were run or errors fixed unless you actually executed and documented them. Never optimize the study to obtain p < 0.05, hide failures, reinterpret exploratory evidence as prospective confirmation, or misrepresent older code as the current release. Avoid changing original frozen experiments, original outputs, tags, or provenance.

You may work in an isolated branch, create tests, run computations within the authorized environment, and prepare a reviewable PR. **Do not merge, overwrite frozen artifacts, force-push, or publish manuscripts without approval.** Preserve negative results as first-class scientific evidence.

### 0 — Establish the actual state before anything else

1. `git fetch --all --tags`; record UTC timestamp, all relevant branch heads, commit SHA, `git status --porcelain`, upstream divergence, dependency locks and available datasets. Check `main`, `v1.0-dev`, any local `audit-2026-10-12` branch, and the status of PR #2. If uncommitted work exists, protect it; do not silently discard, rebase or overwrite.
2. Read `CLAUDE.md`; `README.md`; `docs/EXPERIMENT_PLAN.md` INCLUDING amendments; `docs/CHECKPOINT.md`; `docs/FEYNMAN_REVIEW.md`; `docs/INDEPENDENT_ADVERSARIAL_AUDIT_2026-10-12.md`; `docs/POSITIVE_RESULT_FALSIFICATION_PLAN_2026-10-13.md`; `docs/AUDIT_TOPOQSPR.md`; `docs/INDEX_VERIFICATION.md`; `docs/EVALUATION_PROTOCOL.md`; `docs/EXPERIMENTAL_DATA.md`; `docs/PIPELINE.md`; `docs/LITERATURE_GAP_2026.md`; `docs/REPLICATION_DATA.md`; `docs/PAPER_ROADMAP.md`; paper sources and bibliography.
3. Inventory all active code, tests, data, raw-source hashes, archived results, run records, experimental configurations, scripts, changes since frozen tag `plan-freeze-v1`, and claims in the manuscript. Do not assume `main` root `ml.py` is the v1.0 implementation: the README says it is a flawed v0.x prototype; find the actual `src/topoqspr/` package.
4. Check whether October 13 E2 work is pushed. Do not treat a local report described by a proposal as an independently accessible result. Verify every reported result against actual archived predictions.
5. Write `docs/research/00_STATE_AND_CLAIM_LEDGER.md` with: claim ID, exact textual claim, file:line, underlying output, source SHA, evidence status [CONFIRMED / REFUTED / NOT VERIFIED], possible counterexample and next action. Reconcile every paper number to a run record.

### 1 — Research question and novelty: reject weak contributions

Default *candidate* question, to be accepted only if prior-art checks and feasibility survive:

> Under rigorous leakage controls, when does unlabeled molecular topology supply predictive information beyond molecular size and elemental composition, and what are the mathematical and empirical limits imposed by chemically unlabeled graph representations?

Do not frame the work as "first platform to use 104 indices", "first disease-driven QSPR", "first to combine descriptors and ML", "new atom-labeled indices", or "topological indices improve all properties": investigate whether each is already in the literature. **Explicitly separate** (A) methodological benchmark/replication contribution, (B) proven mathematical representation limit, (C) experimentally verified incremental value, and (D) reusable validated research software. These are different types of novelty and require different evidence.

Candidate manuscript title ONLY IF supported: **"What Does Molecular Topology Add? Representation Limits and Leakage-Resistant Benchmarks for QSPR"**.

Possible routes, in order:
- Route I: a rigorous research paper quantifying the representation ceiling, conditional incremental value, failures under scaffold/similarity transfer, and replication of published claims across independent measured-property datasets.
- Route II: a reproducibility/methodology paper documenting where conventional topological-index QSPR results collapse under proper molecular identity, composition baselines, nested CV, and external validation. Scientifically valid negative results are welcome.
- Route III: a software article on the independently verified, reproducible disease-cohort-to-QSPR workflow, if Routes I/II cannot meet novelty requirements. Must demonstrate features distinct from ChemGraphX, QSPRpred, DOPtools and related tools.
- Stop a claimed *novel descriptor family* if substantially equivalent constructions already exist.

### 2 — Systematic, source-grounded literature audit (through October 2026)

Search Crossref, PubMed/Europe PMC, OpenAlex, Semantic Scholar, Google Scholar, arXiv, ACS, RSC, Springer/Nature, and mathematics databases such as MathSciNet/zbMATH **where legally/accessibly available**. Record exact query, date, source, inclusion criteria, DOI, peer-reviewed/preprint status, and duplicate records. Read original methods and supplementary material—not only abstracts. Never invent search results or claim an exhaustive review of inaccessible databases.

Mandatory starting papers/tools, then expand citation graphs forward/backward:
- ChemGraphX (Jacob et al., *Journal of Cheminformatics*, 2026), DOI 10.1186/s13321-025-01141-x.
- QSPRpred (van den Maagdenberg et al., *Journal of Cheminformatics*, 2024), DOI 10.1186/s13321-024-00908-y.
- DOPtools (Byadi et al., *Digital Discovery*, 2025), DOI 10.1039/D4DD00399C.
- QSAR in the Browser (Masud et al., *JCIM*, 2026), DOI 10.1021/acs.jcim.6c01010.
- MoleculeNet (Wu et al., *Chemical Science*, 2018), DOI 10.1039/C7SC02664A.
- SIMPD (Landrum et al., *J. Cheminformatics*, 2023), DOI 10.1186/s13321-023-00787-9.
- "Scaffold Splits Overestimate Virtual Screening Performance" (Guo et al., 2024 preprint): check peer-reviewed successor and applicability to continuous endpoints rather than assuming generality.
- Ash et al., "Practically Significant Method Comparison Protocols for Machine Learning in Small Molecule Drug Discovery", *JCIM*, 2025, DOI 10.1021/acs.jcim.5c01609.
- Ramakrishnan et al., QM9, *Scientific Data*, 2014, DOI 10.1038/sdata.2014.22.
- Ivanciuc, "QSAR Comparative Study of Wiener Descriptors for Weighted Molecular Graphs" (2000), DOI 10.1021/ci000068y; and 2000 weighting-scheme paper DOI 10.1021/ci990129w.
- "QSPR analysis of some novel neighbourhood degree-based topological descriptors" (*Complex & Intelligent Systems*, DOI 10.1007/s40747-020-00262-0): pre-existing isomer-degeneracy studies.
- "From graph theory to chemoinformatics: modified bond-based indices and a hypothesis-driven multi-task QSAR/QSPR benchmark" (*Scientific Reports*, published 2026-02-21), DOI 10.1038/s41598-026-40969-7: especially strong overlap with a naive "new bond index + multi-task QSPR" idea.
- Foundational molecular descriptor work (Wiener, Randić, Zagreb, Sombor, Hosoya, Balaban, spectral/weighted/atom-labeled indices), original Arockiaraj-group papers, PADEL/Mordred/alvaDesc, modern ECFP/Morgan, graph kernels and molecular GNN baselines.

Produce a **prior-art claim-by-claim matrix**: paper, descriptor math, chemical labels/bond orders, data, endpoints and whether measured, identity handling, split strategy, nested selection, baselines, uncertainty, external validation, code and actual differentiator. Cite source locations. Search for exact mathematical equivalence under different notation, not merely matching titles. Decide NOVEL / INCREMENTAL / ALREADY KNOWN / UNRESOLVED for every proposed contribution, with references.

### 3 — Mathematical and graph-theoretic correctness

Verify the exact definition, hypotheses, domain, units/dimensional interpretation, asymptotic computational complexity, and counterexamples for every active index. Check published primary definitions, adjacency/degree conventions, hydrogen handling, distances, bond orders, salts, disconnected/zero-degree graphs, reverse degrees, neighborhood sums, and modified reverse-degree indices. Confirm whether the current active registry is 102 or 104 in every experiment; do not retroactively mix registries. Cross-check with independent NetworkX/Mordred/ChemGraphX implementations, brute force on graphs with <= 8 vertices and analytical values for paths, stars, cycles, complete graphs, regular graphs, trees, disjoint unions and edge cases. Use property-based tests and exact arithmetic when possible.

Seek: accidental duplicated/affinely dependent indices; identities implied by n/m/degree sequence; degeneracy and collisions; mistaken denominators; NaN silently converted to zero; overflow (especially Hosoya matching counts); performance and cache errors; inconsistent chemical graph standardization; wrong bond multiplicities. Try to falsify existing cross-check claims and classify any disagreement [CODE BUG / DEFINITION DIFFERENCE / INPUT NORMALIZATION / SOURCE TYPO / FLOATING POINT / UNRESOLVED].

Investigate and, where justified, formally prove the following **basic information-limit facts**; do not pretend they are necessarily new mathematics:
1. If two molecules have isomorphic *unlabeled hydrogen-suppressed heavy-atom graphs*, every exact invariant of that graph agrees, irrespective of atom labels and bond-order arrangements. Provide a graph-isomorphism proof.
2. A model whose *entire feature vector* is identical on two molecules must output the same prediction for both under a deterministic fixed fitted predictor. This does NOT mean all constitutional isomers have identical graph invariants.
3. Under squared loss, the optimal population prediction using representation T has irreducible Bayes risk E[Var(Y | T)]. Derive the conditional-variance identity and explain what collision groups can and cannot establish from finite/noisy samples.
4. In a connected graph, cyclomatic number = m-n+1. Check whether each alleged "ring count" feature adds genuinely independent information when size is already supplied.
5. Distinguish descriptive redundancy, exact algebraic dependence, training-set multicollinearity, and genuinely independent predictive information. Report numerical rank only with tolerance, scaling and sample-size disclosures.

**Reproducible counterexample:** ethanol `CCO` and dimethyl ether `COC` have formula C2H6O and the same unlabeled heavy-atom path P3 but different atom placements. Their RDKit-computed Crippen MolLogP may differ; verify with the installed RDKit and cite version. This illustrates *representation collision*, not evidence that new indices predict experimentally measured endpoints better. Include counterexamples where same molecular formula has different unlabeled topology, and where more descriptor discrimination *worsens* generalization. Do not confuse constitutional isomers with stereoisomers.

Optional novel mathematical work: propose a precise and nontrivial theorem about descriptor collisions or bounds that is demonstrably absent from prior art, independently prove or counterexample it, and provide Lean 4 verification only if it adds meaningful rigor. Do not generate decorative theorems.

### 4 — Data provenance and chemical identity audit

Inspect `src/topoqspr/{chembl,cohort,graph,expdata,target_providers,pipeline}.py`, raw-source manifests, hashes, and selection rules. Independently re-download/check pinned sources when licensing/connectivity permits. Check:
- ChEMBL version, exact MeSH disease ID, indication-specific versus molecule-wide phase, DailyMed-only exclusion, approved-small-molecule definition, evidence provenance and possible false positives.
- Parent/salt stripping, mixtures, solvates, stereochemistry, tautomers, protonation, isotopes, full InChIKey versus connectivity first block, molecule-group and assay-source duplicates.
- Stable ordering and deterministic cohort assembly; whether full source forms must be retained for melting point and solid-state solubility. Test with salts, stereoisomers, metformin formulations, and structure/assay counterexamples.
- What OPERA/PHYSPROP, Bradley, AqSolDB and AstraZeneca records *actually measure*. OPERA is also a predictive modeling platform: independently establish whether each specific file contains observed measurements or predictions/curated external values. Verify source IDs, units, temperatures, pH and transformations, censoring, replicate spread, inter-lab heterogeneity and data licences.
- Missingness related to chemical size, class, or property; sample-selection bias from restricting to approved drugs and matched target coverage; overlap between training target sources and descriptor-estimator training sources.
- Cross-dataset chemical identity leakage and near duplicates. Do not join measured values merely by the first InChIKey connectivity block without a justified chemical identity assumption.

Create `docs/research/01_DATA_PROVENANCE_AUDIT.md` with an ID-level discrepancy list and a machine-readable audit manifest. Do not replace experimental measurements with calculated or predicted labels without a clear separate category.

### 5 — Code audit and destructive tests

Focus on `src/topoqspr/{indices,features,models,evaluate,stats,pipeline,experiments,expdata}.py`, CLI, website, tests and CI. Trace each data-dependent operation: filtering, scaling, imputation, feature selection, duplicate removal, target transformations, applicability domain, model-family selection, grid search, permutation, bootstrap and graph clustering.

MANDATORY adversarial issues:
- A1: original six-family `selected` model versus Ridge-only and amended SVR-only permutation/bootstrap controls. Check whether inference statements incorrectly treat family-specific controls as validation of selected-family selection.
- A2: genuinely paired `TI+composition` vs `composition` incremental-value contrast, properly labeled exploratory; do not confuse with `TI` vs `size`.
- A3: group-preserving permutations: a vector permutation does NOT establish exchangeability under within-group outcome correlation. Simulate clustered null data with unequal group sizes and test Type I error and alternative valid randomization schemes; specify assumptions.
- A4: `stats.corrected_ttest` currently implements zero observed variance and nonzero mean as p=0 with a point-width interval. This is a **confirmed code convention**, not proof of certainty. Reproduce with `[0.1,0.1,0.1]`, discuss statistical contract and report warnings or undefined inference where justified; do not silently change published statistics.
- A5: 102-index frozen experiment versus 104-index later control, and dirty working-tree status. Trace exactly which results use which registry and whether comparison is genuinely matched.
- A6: do not call a Git tag created before the run an externally preregistered OSF protocol unless verifiable dated registration exists.
- Ridge numerical explosion for extreme large cardiac glycosides; check shape, solver convergence, regularization, extrapolation and whether predictions are bounded under amended SVR. Preserve original failure metrics.
- Oversized/scaffold clusters, tiny test folds and variance-undefined R²; improper treatment of dependent repeated folds as independent n; zero/nonfinite targets; silent `NaN` omissions; p-value multiplicity, y-randomization, sign conventions.
- Re-run targeted attack tests on v0.x: Hosoya >18 bonds, Balaban salts, wrong index numerators, shared estimator instance, full-data hyperparameter selection, ordered/unordered drug slicing, duplicate formulas. Confirm none persist in active v1.0 rather than assuming an old audit applies.
- Property-based/mutation tests: shuffle y; randomly permute IDs; reorder input molecules; duplicate stereoisomers; inject missing targets/features; remove an expected index; deliberately alter an index formula; hold out a target to test preprocessing integrity; compare multi-CPU runs and thread counts.
- API and website: actual v1.0 startup and run end to end, sanitization, timeout/jobs, concurrent sessions, insecure paths, resource exhaustion, filesystem containment, error diagnostics, clear separation from legacy v0.x.

Every bug report must include file:line, root cause, minimal failing fixture, expected versus actual, scientific consequence, unit/regression test, and before/after metrics. Mark a suspicion FALSE when evidence refutes it. Do not run Codex and say "Codex verified" unless actually invoked.

### 6 — Statistical review from raw out-of-fold predictions

Independently recompute every original confirmatory endpoint and experimental table: n, exact folds, R² definition, MAE/RMSE, paired delta sign, fold and pooled metrics, mean, corrected repeated-CV standard error, intervals, p-values, BH/BY q-values, adjusted family size and practicality threshold (ΔR² >= 0.05). Reproduce one by hand on a small fixture. Compare source hashes with the manuscript macros. Investigate coverage of corrected-t inference and bootstrap under correlated/repeated resampling, unequal group folds and multiple endpoints; report inference as approximate where necessary.

Check independent sample unit (molecule versus chemical series versus formula versus assay), repeated CV reuse, group leakage, uncertainty across datasets, target-selection bias, outer seed uniqueness, applicability-domain selection, calibration, outlier sensitivity and influential observations. Try a principled cluster bootstrap/hierarchical analysis if justified. Distinguish *failure to detect benefit* from *evidence of equivalence* (requires a prespecified margin and an adequately powered equivalence procedure). Never substitute a family-specific SVR interval for the selected-family model's interval.

Original archive to independently reconcile:
`results/expB/main_confirmatory/confirmatory.csv` currently reports `TI-size`, selected learner:
- mp_bradley n=711, ΔR²=+0.071005, BH q=0.075828
- logp_opera n=505, ΔR²=+0.152927, BH q=0.001553
- logd74_az n=369, ΔR²=+0.038623, BH q=0.348587
- logs_aqsoldb n=351, ΔR²≈−0.012533, BH q≈0.673940

These four statistics are **not** direct proof that TI adds information beyond full chemical composition. A later October-13 proposal reports exploratory E2 results but says they came from local, then-unpushed work; locate and verify actual files before quoting them as independently reproducible.

### 7 — Candidate new experiment: representation and isomer stress test (NOT yet approved)

Do **target-blind feasibility before model fitting**:
- For each of the four measured endpoints, census formula groups with at least 2/3/5 standardized distinct constitutional structures, number of observations, source/assay consistency, number of chemical identity duplicates, target variability and chemistry class.
- If current approved-drug subsets are too sparse, say NO-GO. Consider well-documented external *measured* continuous endpoints with compatible units and chemical space; record source and licensing.
- QM9 C7H10O2 (6,095 constitutional isomers) is a candidate **computed-quantum-property mechanistic benchmark** only. It is not a measured-drug logP, solubility or melting-point benchmark and cannot independently confirm such claims.

Design three separate estimands, not an accidental mixture:
**E3a — known formula, unseen constitutional structure:** assess ability to discriminate structural variation within represented molecular formulas, with leakage-free molecular grouping and exact source identity.
**E3b — unseen molecular formula:** group-disjoint transfer to genuinely novel formula groups; do not use formula-level target means from the test set.
**E3c — chemical-dissimilarity transfer:** compare scaffold grouping with verified similarity-cluster split; quantify train/test Tanimoto distribution. Neither scaffold nor clustering guarantees external validity.

Feature ladder evaluated on identical rows/folds/search budgets:
(0) train-mean/null;
(1) heavy atoms + bonds;
(2) full elemental composition including hydrogen;
(3) existing version-locked unlabeled TI;
(4) composition + unlabeled TI;
(5) chemical labels/bond-aware TI if mathematics and prior art justify;
(6) RDKit 2D excluding direct target leakage;
(7) Morgan ECFP fingerprints and a capacity-matched compact fingerprint;
(8) RDKit+TI and ECFP+TI incremental tests;
(9) optional efficient GNN only with enough data/compute and a meaningful controlled comparison.

For representation collisions report collision rates, conditional target spread, per-formula MAE/RMSE, pairwise errors with dependence-aware uncertainty and non-equivalence caveats. A representation distinguishing more molecules is not by itself evidence of predictive improvement. Within a fixed-formula group, formula-only features are constant; choose metrics well-defined when outcome variance is tiny. Predefine a primary benchmark, primary contrast and smallest meaningful improvement before looking at performance. Use formula/series-level resampling where appropriate. Archive held-out predictions and cluster assignments.

Test independent-source replication and negative controls; report where adding TI worsens accuracy. No retrospective "confirmatory" relabeling of work chosen after prior results were seen.

### 8 — Optional new chemical-label-aware indices: mathematical work before programming

Only pursue if literature audit discovers a specific unresolved gap and data feasibility passes. Mathematically define a minimal family of atom/bond-labeled graph invariants with exact reduction to existing unlabeled indices in the unlabeled limit. Demonstrate relabeling invariance under chemically label-preserving graph isomorphism, domain, numerical stability, computational complexity, degeneracy and informative non-isomorphic witnesses. Compare against **old weighted Wiener and atom-weighted graph indices**, known atom-pair descriptors, Morgan ECFP and RDKit descriptors. Require proof that the construction is not an algebraic reparameterization already published. Use at least two independent calculation implementations.

A potential positive conclusion must be **incremental improvement beyond composition AND competitive fingerprints on independent measured data**, with appropriate uncertainty; if it fails, publish the limitation and collision/benchmark result instead. Do not claim the 2000 weighted-index literature has been superseded without evidence.

### 9 — Reliability and reproducibility

Use one clean, pinned environment; include RDKit/NumPy/scikit-learn and system versions, deterministic sorting, feature registry hash, raw-source SHA-256, standardized-structure IDs, dataset/manifest hash, outer/inner fold assignments, random seeds, config hash, commit, dirty flag, runtime and exact reproduction command. Keep ordinary tests fast; expensive experiments behind intentional flags and compute budget estimates based on measured smoke timing. Run `pytest`/CI, external independent numerical checks and at least one clean-environment end-to-end reproduction. Declare unexecuted checks rather than implying they passed. Public artifacts must respect data redistribution licenses and consent for original source materials.

### 10 — Manuscript development: conditional, not predetermined

Before writing results, produce a journal-positioning table:
- **Journal of Cheminformatics** only if methodological/representation contribution and competitive baselines are strong and clearly distinct from ChemGraphX and QSPRpred.
- **Digital Discovery** if the chemical ML benchmark and mechanistic insight are sufficiently substantial.
- **SoftwareX** or a suitable software venue if primary deliverable is a reliable, differentiated platform and not a demonstrated new scientific prediction claim.
- Any other venue only after verifying 2026 author guidelines and appropriate article type.

Suggested structure: (1) falsifiable research question; (2) precise mathematical representation and existing-work boundary; (3) source and identity provenance; (4) prespecified statistics/splits; (5) benchmark and independent replication; (6) collisions, negative results, mechanism and limitations; (7) reproducible software/workflow; (8) practical implications.

Every statement in abstract, figures, conclusion and supplement must map to a proof or artifact. Do not plagiarize ChemGraphX structure/text; emulate appropriate academic organization without copying.

### 11 — Deliverables, gates and final verdict

GATE A **state and literature**: no expensive modeling. Deliver `00_STATE_AND_CLAIM_LEDGER.md` and `02_NOVELTY_PRIOR_ART_MATRIX.md`. Mark each recommendation from THIS prompt and the repository's Oct-13 plan VALID / VALID WITH CHANGES / INVALID / UNRESOLVED, with counterexample and citations.

GATE B **correctness and data**: add `01_DATA_PROVENANCE_AUDIT.md`, `03_MATH_COUNTEREXAMPLES.md`, `04_CODE_STATISTICS_AUDIT.md` with failing/passing tests and severity [CRITICAL/HIGH/MEDIUM/LOW]. Stop if primary sources are not valid measurements or if chemical identity cannot be trusted.

GATE C **feasibility**: create `05_TARGET_BLIND_DATA_FEASIBILITY.md` with isomer counts and a defensible GO/NO-GO decision; reject infeasible positive-result hypotheses.

GATE D **protocol**: write `06_PROSPECTIVE_STUDY_AMENDMENT.md` before viewing new test outcomes: independent units, hypotheses, endpoints, contrasts, split hashes, models, exact alpha/multiplicity plan, practical effect, stopping conditions, computation, data rights. Existing E2/E3 exploratory work cannot become retroactively confirmatory.

GATE E **run and independent verification**: only then run justified studies; archive `07_RESULT_RECONCILIATION.md` plus immutable data/prediction manifests, negative controls, all failures, tests, intervals and comparisons. If no runs occur, report "NOT RUN" explicitly.

GATE F **publication**: create `08_PAPER_NOVELTY_AND_SUBMISSION_DECISION.md`, a reproducible journal-appropriate outline, and `09_REVIEWER_TWO_ADVERSARIAL_REPORT.md` written from the perspective of a hostile but fair reviewer. Include strongest evidence **against** novelty and strongest counterexample to each major claim.

Each gate must conclude with a concise independent verdict and concrete citations, exact file locations, test commands and observed output. Do not run later gates automatically if an earlier gate fails. Do not rewrite prior results or change paper statements without traced justification.

### Final response to owner — REQUIRED format

1. Current verified commit/branches and unpushed changes.
2. Top 10 actual blockers ranked, clearly separated as CONFIRMED BUG / UNVERIFIED RISK / STATISTICAL INFERENCE LIMIT / PRIOR ART / DATA FEASIBILITY.
3. Which previous recommendations were wrong, with counterexamples.
4. Exactly which papers preempt the naive novelty claims, and what *specific* unfilled niche remains.
5. Corrected mathematical statements, proofs/counterexamples, exact index count per experiment.
6. Reconciliation of original confirmatory and all exploratory/amended results, without suppressing negative cells.
7. Which tests and experiments you ACTUALLY executed, with commit, run IDs, runtime and hashes.
8. Recommended scientific contribution, whether it currently meets that bar, and an honest [GO / CONDITIONAL GO / NO-GO].
9. Minimal, sequenced changes needed to make a publishable submission.
10. PR link (if created), no destructive merge, and a short paper abstract **only if** evidence warrants it.

**Controlling principle:** "Novel" means not already in prior work; "positive" means supported by a fair out-of-sample comparison, not a selected p-value. A reproducible, insightful negative result can be more publishable than a misleading positive one.