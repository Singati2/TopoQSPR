# TopoQSPR: independent adversarial audit and coding-agent worklist

**Review date:** 2026-10-12  
**Repository:** https://github.com/Singati2/TopoQSPR  
**GitHub HEAD reviewed:** fecb800d2b10fb8802d889b6d0b4c3a192eef090 (2026-10-06)  
**Target:** v1.0 implementation, frozen Experiment B, amended control, Experiments A–D, mathematical/stats code, manuscript.  
**Review limitations:** GitHub source and archived results inspected via connected GitHub API; **not** an execution of Codex, the entire pytest suite, complete data rebuild, or independently repeated nested CV. Do not claim an issue is independently reproduced until its failing test has been run.

## Read this first (for Claude Code or Codex)

Work adversarially, not as a defender of the current manuscript. Begin with docs/EXPERIMENT_PLAN.md (including 2026-10-05/06 amendments), docs/CHECKPOINT.md, docs/FEYNMAN_REVIEW.md, docs/EVALUATION_PROTOCOL.md, docs/INDEX_VERIFICATION.md, the two expB confirmatory directories, and src/topoqspr/{evaluate,stats,pipeline,features,indices}.py. **Keep the original frozen results and the amended SVR controls separate. Do not rewrite numerical conclusions to improve publication prospects.**

For each item below: inspect implementation; try to falsify the concern; write the smallest failing counterexample/test if one exists; identify whether it is a **confirmed defect**, **inference limitation**, **documentation mismatch**, or **additional experiment**; implement the narrowest justified fix; rerun focused tests plus pytest; publish before/after numeric deltas and exact commit SHA. If no defect exists, say why with an explicit test or proof. Preserve ChEMBL 37 data provenance, frozen protocol, and negative results. Do not claim an independent Codex validation from this file.

## Grounded snapshot: what the archived experiment says

Original selected-learner confirmatory comparison (TI vs two-feature size baseline, 5 repeats x 5 folds; four-member BH family):

| Measured endpoint | n | mean paired Delta R2 | corrected 95% CI | BH q | BH positive? |
|---|---:|---:|---|---:|---|
| mp_bradley | 711 | +0.071005 | [+0.004302, +0.137708] | 0.075828 | No |
| logp_opera | 505 | +0.152927 | [+0.076329, +0.229526] | 0.001553 | Yes |
| logd74_az | 369 | +0.038623 | [-0.030688, +0.107935] | 0.348587 | No |
| logs_aqsoldb | 351 | -0.012533 | [-0.073262, +0.048195] | 0.673940 | No |

Data: results/expB/main_confirmatory/confirmatory.csv. A positive Delta R2 vs size does not imply improvement beyond full element composition. The repo's Feynman review calls composition's better absolute performance **descriptive**, because a paired inferential contrast against composition was not run.

### A1 [HIGH: inference scope] Selected-learner hypothesis vs single-family controls

- Evidence: results/expB/main_confirmatory/run_record.json has evaluation.families = ridge, elasticnet, pls, knn, svr, extratrees, protocol.learner = selected, but protocol.permutation.learner = ridge and protocol.bootstrap.families = [ridge]. See src/topoqspr/pipeline.py lines 1145-1172 and 1479 onward.
- **Known issue, not a newly found bug.** docs/EXPERIMENT_PLAN.md explicitly documents the post-run amendment. Original Ridge-only bootstrap produced extreme negative intervals (for logP approximately [-1,662,029, +0.190]); original Ridge-only permutation called some TI sets uninformative. The later SVR-only run is archived in results/expB/main_confirmatory_protocol_svr and reports logP bootstrap Delta R2 interval approximately [+0.046, +0.189].
- Follow-up: establish whether full six-family selection controls are computationally feasible (or use a transparent reduced-cost, independent confirmatory replication). If not, explicitly label both Ridge and SVR as *family-specific sensitivity analyses*, not whole-selected-procedure confirmation. Do not alter or discard either result.
- Acceptance: manuscript table, abstract, captions and source-macro generation all distinguish the selected-model main test from Ridge and SVR controls; automated test verifies correct learner names and feature-registry versions.

### A2 [HIGH: missing scientific contrast, not pre-freeze bug] Incremental graph information beyond composition

- Evidence: src/topoqspr/pipeline.py lines 775-777 includes TI_composition; docs/EXPERIMENT_PLAN.md lines 160-162 prespecifies the secondary contrast; confirmatory family tests TI vs size.
- Add **exploratory** paired TI_composition vs composition and TI vs composition on identical outer folds with consistent inner family selection, corrected-CV intervals and multiplicity explicitly identified. Assess sensitivity to standard RDKit 2D baseline. Evaluate whether a positive result survives scaffold or chemically dissimilar holdout.
- Counterexample conceptual check: different chemical structures can share heavy-atom/bond and elemental counts, so an index might add useful connectivity; conversely the reported TI advantage over size may be entirely composition-related. Either outcome must remain reportable.
- Acceptance: archived scripts, fold-level paired differences, effect sizes, CIs, corrections, train/test split hashes; conclusions explicitly separate 'beyond two size counts' from 'beyond composition'.

### A3 [MEDIUM-HIGH: permutation exchangeability needs proof/test]

- Evidence: src/topoqspr/pipeline.py lines 1012-1048. For groups with a unique size, the algorithm mixes their rows with singleton rows; for repeated equal-size groups it swaps whole target blocks. Code comments say a permutation of y 'stays valid'.
- That last implication is not generally true **without an exchangeability assumption**. A permutation vector is not by itself a valid null randomization: for strongly correlated outcomes within a cluster, row-wise reassignment can destroy within-cluster dependence.
- Construct a synthetic clustered null with group sizes [3,2,2,1,1,...], correlated responses within groups and X independent of y; compare null Type I calibration of current scheme with a principled group-respecting permutation, keeping the same statistical target. Assess whether the published grouped runs are affected before proposing changes.
- Acceptance: explicit null/exchangeability assumptions in docs; test of permutation indices and simulation-based calibration; when assumptions fail, flag rather than silently report a valid p-value.

### A4 [MEDIUM: degeneracy / numerical-inference convention]

- Evidence: src/topoqspr/stats.py lines 99-112; tests/test_stats.py lines 47-51 intentionally require a nonzero constant vector of fold differences to give p=0 and zero-width CI.
- Minimal analytic counterexample: corrected_ttest([0.1,0.1,0.1], n_train=8, n_test=2) yields p=0, CI=[0.1,0.1]. This is a *code-confirmed convention*, not evidence of generalization certainty. On overlapping CV folds, zero empirical dispersion need not justify a population point interval.
- Investigate whether to retain the mathematical limit with an explicit 'degenerate observed fold variance' warning, or produce undefined inferential p/CI rather than p=0. Decide based on statistical contract and document effects; do not silently change published values.
- Acceptance: deterministic test for this condition, defined behavior in stats documentation, clear downstream inferential labeling.

### A5 [MEDIUM: release/provenance and mixed descriptor versions]

- Evidence: original main_confirmatory run_record reports TI=102, code commit 659fac38, dirty=false; amended SVR run_record reports TI=104, commit fe92837, dirty=true. docs/EXPERIMENT_PLAN.md explicitly explains two added indices for Experiment D and the amended control's 104-index run. Manuscript abstract currently says '104 indices'.
- This is **documented**, not an undisclosed index-computation bug. Verify every table and figure identifies which registry version produced it; compare 102-vs-104 main-fold selected-model effects where appropriate, without redefining frozen confirmatory inference.
- Acceptance: machine-readable run-to-table provenance ledger, clean pinned versioned rerun for release where feasible, CI checks that the manuscript's '104' refers to the software and not falsely to the original 102-index confirmatory analysis.

### A6 [HIGH: integrity of manuscript preregistration terminology]

- Evidence: docs/EXPERIMENT_PLAN.md heading: plan frozen on a timestamped git tag on 2026-10-05, but 'OSF registration is pending'; paper/main.tex describes 'preregistered statistics'. External OSF registration was not shown in files reviewed.
- Fix terminology until a public prospective registration with verifiable timing exists. An internally timestamped, pre-result git plan is *pre-specified/frozen*, not automatically an externally registered preregistration. Any OSF record submitted **after** results should be called retrospective registration rather than prospective preregistration.
- Acceptance: manuscript, abstract, README, and registration status all agree with actual public records; include verifiable URL and creation date if one exists. Preserve original freeze provenance.

### A7 [MEDIUM: robustness of chemistry and splits]

- Evidence: docs/CHECKPOINT.md confirms stereo-connectivity grouping; src/topoqspr/evaluate.py lines 696-717 groups by Bemis-Murcko scaffold, assigning acyclic molecules separate groups. Experiment C/scaffold results exist.
- Add, as exploratory external-validity work, a fingerprint similarity or Butina-cluster holdout in addition to (not as replacement for) existing random and scaffold group experiments. Verify no near-identical stereoisomer/tautomer/standardized-parent leakage; check groups before any target imputation or tuning.
- Acceptance: manifests with group identifiers and largest-group proportions; verified no group split across training/test; estimates and uncertainty side-by-side without cherry-picking favourable folds.

### A8 [MEDIUM: outlier/extrapolation diagnosis]

- Evidence: docs/FEYNMAN_REVIEW.md corrects original 'three cardiac glycosides' statement to four named outliers: lanatoside C (logP), digoxin and digitoxin (solubility), atracurium (melting point). See tools/check_extrapolation.py and original/bootstrap CSV.
- Reproduce residual >10 target-SD rule on actual endpoint-specific training folds and report indexwise train-range ratios, predictions, and breakdown by learner. Do not remove molecules from pre-specified all-molecule primary metrics. Compare robust or bounded prediction models only as explicit sensitivity.
- Acceptance: one archived molecule-by-fold witness table and regression tests for the outlier-identification script, preserving the corrected names.

### A9 [MEDIUM: data/cohort and publication audit]

- Independently verify a stratified sample of ChEMBL 37 MeSH indication evidence, DailyMed-only exclusion, phase=4 semantics, structure parents, InChIKey exact matches, external MP/logP/logD/solubility source licences and units. The existing 1,821-parent cohort and 711/505/369/351 endpoint coverages come from the repo; this audit has *not* rebuilt them.
- Evaluate 3-study replication sample selection against predeclared eligibility, and label very broad confidence intervals for studies with n around 10-25 as inconclusive rather than proof of no added value.
- Acceptance: evidence table of source identifiers, retrieval timestamps, versioned molecule hashes, and independently checked random sample; all inclusion/exclusion counts reproducible.

## Priority implementation sequence

1. **First pass (no number changes):** reconcile protocol-versus-controls claims, clarify preregistration terminology, add run-to-manuscript provenance table, add degenerate-test explanation.
2. **Adversarial tests:** group permutation exchangeability toy/calibration, constant-difference corrected_ttest, feature-version consistency, 4 outlier witnesses and bootstrap provenance.
3. **New exploratory experiments:** composition-adjusted paired comparisons and stronger chemical-disjoint validation; preserve original pre-specified family and no post hoc rebranding.
4. **Release engineering:** full pytest/CI, locked environment clean tree, deterministic regeneration, code/figure/statement matrix, correct licences/DOIs, coauthor review.

## Required agent completion report

Return a Markdown table: issue ID | claim tested | exact source path:line | failing counterexample/test | verdict (confirmed / non-bug / limitation / new experiment) | fix commit | tests run | numbers changed? Include output of pytest, git diff --stat, git status, run-record SHA and manifest hashes. Never imply tests ran if they were merely read. Do not push runtime-heavy changes or alter paper's frozen results without recording an amendment, results-diff ledger, and justification.