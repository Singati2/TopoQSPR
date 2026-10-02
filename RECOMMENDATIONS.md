# TopoQSPR — Publication Recommendations

## Objective

Continue TopoQSPR, but do **not** add new methods or broaden the project yet. The next phase should focus on making the current scientific claim publication-grade.

The central question should remain:

> **When disease-specific drug cohorts are rigorously defined and evaluation is leakage-free, how much predictive information do graph-theoretic topological indices carry beyond molecular size?**

The study should be optimized for a result that survives skeptical peer review, not for a favorable conclusion.

## 1. Finish the audit and freeze the experiment protocol

Finish the code/methodology audit before inspecting any main-panel model results. Debug only on the quarantined pilot disease. Once the pipeline is correct:

- freeze `docs/EXPERIMENT_PLAN.md`;
- record the exact commit SHA and date;
- define all primary/secondary endpoints before the main run;
- record every later change in the Amendments section with date and reason;
- do not change disease selection, targets, baselines, metrics, or statistical tests after seeing main-panel results unless the change is transparently documented as post hoc.

## 2. Independently verify every topological index

Before modelling, resolve all mathematical/implementation issues, especially:

- Hosoya Z;
- DS1/DS2 duplication;
- reverse-degree definitions;
- missing-value behavior;
- any discrepancy between mathematical definitions and implementation;
- exact duplicate or algebraically redundant descriptors.

Add unit tests using analytically tractable graphs/molecules and published reference values. Cross-check overlapping indices against an independent implementation such as ChemGraphX. Failed index calculations must be represented as missing/errors and must never be silently replaced by zero.

Maintain a table containing, for every retained index: mathematical definition, literature citation, implementation location, unit test/reference case, and verification status.

## 3. Make disease cohorts deterministic and defensible

Use controlled disease terms and a reproducible disease-to-drug pipeline:

- exact controlled disease heading rather than substring matching;
- indication-specific clinical-phase filtering;
- small molecules only for the primary analysis;
- parent/salt normalization;
- molecular standardization;
- InChIKey deduplication;
- deterministic ordering;
- fixed and recorded ChEMBL release;
- complete exported cohort manifest for every disease.

The exact drug cohort used in every analysis must be reconstructable from the archived manifest.

## 4. Treat molecular size as the critical confounder

The central scientific test is not whether TIs produce a high raw R². It is whether they contain useful predictive information **beyond molecular size**.

Compare at minimum:

1. null model;
2. molecular-size baseline (heavy-atom count and bond count);
3. size+ sensitivity baseline (e.g. ring count added);
4. verified topological indices;
5. individual TI families;
6. standard RDKit 2D descriptors, excluding descriptors that directly define the prediction target.

Make paired **ΔR²(TI − size)** the primary effect. Also report uncertainty intervals and paired ΔMAE. Raw R² should remain descriptive rather than being the headline result.

Quantify how much variance in individual TIs and TI families can itself be explained by molecular size.

## 5. Use leakage-free evaluation

All data-dependent decisions must occur inside the appropriate training fold, including:

- preprocessing;
- scaling;
- missing-value handling;
- constant/near-constant feature removal;
- feature selection;
- hyperparameter tuning;
- model-family selection.

Use the pre-specified repeated nested CV protocol. In addition, run scaffold-grouped outer validation using Bemis–Murcko scaffolds to measure analog/scaffold leakage.

Explicitly report the change:

`ordinary/random-split performance -> nested-CV performance -> scaffold-CV performance`.

## 6. Run the multi-disease panel exactly according to the frozen rule

Do not select diseases because their model results are favorable. Apply the pre-specified cohort-size rule before modelling and retain all eligible diseases.

Negative findings, failed generalization, small effect sizes, and diseases where TIs do not outperform size are results and should be retained.

## 7. Prioritize the published-study replication experiment

This may become the strongest scientific contribution.

Select 3–5 published disease-specific topological-index QSPR studies using pre-specified eligibility criteria. The selected studies should have sufficiently recoverable compound lists/property values and clear index definitions.

For each study, reconstruct as faithfully as possible:

`published R²`

→ `reproduced in-sample R²`

→ `ordinary CV R²`

→ `nested-CV R²`

→ `scaffold-CV R²`

→ `molecular-size baseline`

Report uncertainty intervals and quantify the optimism gap. Determine whether the reported TI advantage survives rigorous validation and whether any apparent performance is largely recoverable from molecular size alone.

Do not choose replication papers based on whether they produce the desired conclusion.

## 8. Strengthen the statistical analysis

Report paired fold-level differences, effect sizes, uncertainty intervals, and multiplicity-adjusted inference where appropriate.

For the main disease × target comparisons:

- use the pre-specified paired repeated-CV analysis;
- report 95% confidence intervals;
- apply Benjamini–Hochberg/FDR control to the confirmatory family of comparisons;
- distinguish confirmatory analyses from exploratory analyses;
- report sample size for every endpoint;
- characterize TI collinearity and family structure;
- avoid interpreting statistical significance without effect magnitude.

Where assumptions of a proposed inferential procedure are questionable, document them and include an appropriate sensitivity analysis.

## 9. Keep SHAP secondary

Do not use SHAP to make a weak predictive model appear scientifically important.

First establish that the model generalizes and that the TI feature set improves meaningfully over the molecular-size baseline. Only then interpret TI families.

Compute interpretation only within properly trained outer-fold models and report stability of attribution across folds/repeats. If attribution is unstable, report that instability rather than selecting a favorable fold.

## 10. Make reproducibility publication-grade

Every final number in the paper should be generated from a tagged release and traceable to:

- Git commit SHA/tag;
- package/software version;
- ChEMBL release;
- disease cohort manifest;
- configuration file;
- random seeds;
- dependency lock/environment;
- model/CV specification;
- output artifact.

Add/complete:

- `pytest` coverage for critical calculations;
- GitHub Actions CI;
- pinned environment;
- package/CLI interface;
- one-command or clearly scripted reproduction workflow;
- provenance/run-record files;
- archived release artifacts (e.g. Zenodo) for the submission version.

Do not manually copy model-performance numbers into the manuscript when they can be generated programmatically.

## 11. Keep the novelty claim narrow and defensible

Do **not** claim novelty merely because TopoQSPR computes many topological indices or combines topological descriptors with machine learning.

The stronger potential contribution is the combination of:

1. rigorous disease-anchored cohort construction;
2. independently verified TI computation;
3. explicit measurement of predictive information beyond molecular size;
4. leakage-free nested and scaffold-aware evaluation;
5. quantification of evaluation optimism;
6. systematic replication/audit of conventional disease-specific TI-QSPR conclusions;
7. a reproducible open-source implementation of the complete workflow.

Position TopoQSPR relative to existing descriptor/index software as complementary rather than claiming functionality that already exists elsewhere.

## 12. Let the evidence determine the final paper

Do not force a positive TI story.

If TIs strongly and reproducibly outperform size across diseases, targets, repeated nested CV and scaffold splits, characterize **when and why** that happens.

If the TI advantage largely disappears after controlling for molecular size or using rigorous validation, make that finding central. Quantify how much apparent conventional TI-QSPR performance is attributable to:

- molecular-size confounding;
- cohort construction;
- model-selection/evaluation optimism;
- analog/scaffold leakage.

A rigorous negative result can be more scientifically valuable than another small-data QSPR paper reporting a very high in-sample R².

# Publication-readiness checkpoint

Before launching the expensive main-panel experiments, produce a checkpoint report containing:

- unresolved code errors;
- unresolved mathematical-definition issues;
- verified vs unverified TI list;
- independent cross-check status;
- disease-cohort validation status;
- leakage audit;
- frozen hypotheses and endpoints;
- exact disease-selection rule;
- final baselines;
- final CV/statistical protocol;
- candidate replication papers and eligibility rationale;
- remaining threats to internal validity;
- remaining threats to external validity;
- reproducibility status;
- **GO / NO-GO decision for the main experiment**.

Do not start the main experiment until the checkpoint is GO.

# Recommended order of work

**Audit → mathematical/index verification → cohort pipeline → pilot/debugging → freeze protocol → main multi-disease experiment → scaffold sensitivity → published-study replication → statistical synthesis → manuscript.**

Avoid adding extra deep-learning models, graph neural networks, quantum methods, or additional descriptor families unless the frozen core experiment reveals a specific scientific reason for them. Additional complexity should not replace the central test of whether TIs provide information beyond simple molecular size under rigorous validation.
