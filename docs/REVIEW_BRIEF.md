# Claude Prompt — TopoQSPR End-to-End Audit, Validation, and Paper-Strengthening Plan

Repository: `https://github.com/Singati2/TopoQSPR`

## Role

Act as a skeptical senior researcher in:

- chemical graph theory,
- cheminformatics,
- QSPR/QSAR,
- statistics and machine learning,
- scientific software engineering,
- reproducibility and research methodology.

Do **not** assume that the existing code, manuscript, README, roadmap, equations, or claimed novelty are correct. Verify them independently.

The goal is **not** to make the project look publishable. The goal is to determine what is actually correct, identify what is wrong or weak, fix the scientifically important problems, and turn TopoQSPR into the strongest defensible paper possible.

---

# 1. Read the Entire Repository First

Audit the current repository end-to-end, including at minimum:

- `README.md`
- `docs/PAPER_ROADMAP.md`
- `paper/main.tex`
- `paper/refs.bib`
- `paper/audit_chembl_indications.py`
- `config.py`
- `fetcher.py`
- `disease_validator.py`
- `ml.py`
- `main.py`
- `drug_pipeline_expanded.py`
- `website.py`
- tests, scripts, data-processing utilities, deployment files, and any recently added files
- git history and the most recent commits

Do not rely only on comments or README descriptions. Trace the actual code execution.

Create an audit table with these categories:

| Item | Status | Evidence | Severity | Required fix |
|---|---|---|---|---|
| Correct | verified | file/function/line | — | none |
| Bug | verified error | file/function/line | critical/high/medium/low | fix |
| Methodological weakness | statistically/scientifically weak | evidence | severity | fix |
| Unverified claim | insufficient evidence | evidence missing | severity | validation needed |
| Manuscript/code mismatch | inconsistent definition or behavior | both locations | severity | reconcile |

---

# 2. Verify the Mathematical Definitions of Every Topological Index

The repository currently implements a large family of degree-, reverse-degree-, neighborhood-degree-sum-, distance-, cut-based-, and entropy-related indices.

Do **not** assume that an implementation is correct just because the function runs.

For every index or index family:

1. Identify the exact mathematical definition used.
2. Find the primary or authoritative literature source.
3. Compare the published equation with the implementation.
4. Check edge cases.
5. Check whether two differently named indices are actually duplicates.
6. Check whether any implementation silently changes the mathematical definition.
7. Check whether normalization conventions differ across papers.
8. Flag any ambiguous or nonstandard definitions.

Produce a provenance table:

| Index | Mathematical definition | Literature source | Code implementation | Verified? | Test case |
|---|---|---|---|---|---|

Pay special attention to:

- Hosoya \(Z\)
- reverse-degree indices
- \(R_k\) reverse-degree variants
- degree-sum variants
- neighborhood-degree-sum variants
- Balaban \(J\)
- Szeged variants
- Mostar variants
- PI
- Schultz
- Gutman index
- graph entropy
- any index whose name is not standard in the literature

### Important reverse-degree issue

The current manuscript and code may use different reverse-degree conventions.

Do **not** simply force the code to match the manuscript.

Instead determine which definition is intended from the cited literature. For example, distinguish definitions such as

\[
\bar d(v)=n+1-d(v)
\]

from

\[
\bar d(v)=\Delta-d(v)+1.
\]

Document the exact convention for every reverse-degree family and cite the original source.

### Hosoya Z

Independently verify the Hosoya index implementation.

If the current algorithm becomes approximate, returns incorrect values beyond some graph size, or silently substitutes another formula, replace it with a correct algorithm or remove the descriptor from the release.

No approximate value may be presented as the exact Hosoya index.

---

# 3. Build Mathematical Unit Tests

Create rigorous automated tests for the descriptor library.

At minimum test:

- paths \(P_n\),
- cycles \(C_n\),
- stars \(S_n\),
- complete graphs \(K_n\),
- several small chemically valid molecules with published reference values.

Whenever a closed form is known, test the implementation against it.

Examples include:

\[
W(P_n)=\frac{n^3-n}{6}.
\]

Use exact arithmetic where possible.

Tests must fail loudly rather than replacing undefined or failed calculations by zero.

### Missing values

Audit any code that currently does something equivalent to:

```python
fillna(0.0)
```

for failed topological descriptors.

A computational failure, undefined descriptor, and true mathematical value \(0\) are not the same event.

Replace this behavior with explicit missing/error handling and provenance.

---

# 4. Audit Disease-to-Drug Cohort Construction

This is potentially one of TopoQSPR's main scientific contributions.

The current project must not define "drugs for disease X" through an uncontrolled substring search.

Design a reproducible pipeline:

\[
\text{free text}
\rightarrow
\text{controlled disease concept}
\rightarrow
\text{ChEMBL indication records}
\rightarrow
\text{clinical-phase filter}
\rightarrow
\text{small-molecule filter}
\rightarrow
\text{parent normalization}
\rightarrow
\text{salt/fragment handling}
\rightarrow
\text{deduplication}
\rightarrow
\text{versioned drug manifest}.
\]

Verify:

- MeSH mapping,
- EFO identifiers where available,
- indication-specific clinical phase,
- approved vs investigational drugs,
- molecule type,
- parent-compound relationships,
- salts,
- mixtures,
- disconnected structures,
- duplicate structures,
- InChIKey-based deduplication,
- deterministic ordering,
- ChEMBL version,
- caching,
- provenance.

Every final cohort should have a machine-readable manifest containing at least:

- disease label,
- controlled identifier,
- ChEMBL release,
- ChEMBL ID,
- preferred compound name,
- canonical structure identifier,
- indication-specific phase,
- parent molecule,
- inclusion/exclusion decision,
- reason for decision.

The exact same configuration should regenerate the same manifest.

---

# 5. Audit the Current Machine-Learning Evaluation

Trace `run_ml_qspr()` and all associated model-selection code.

Verify whether hyperparameters are tuned using any information from the outer test folds.

The evaluation must be **genuinely nested**.

For each outer split:

\[
D =
D_{\text{train}}^{(k)}
\cup
D_{\text{test}}^{(k)}
\]

hyperparameter selection and model-family selection must occur using only

\[
D_{\text{train}}^{(k)}.
\]

Then:

1. inner CV selects hyperparameters/model,
2. refit using the entire outer-training partition,
3. evaluate once on the untouched outer-test partition.

Do not:

- tune once on the full dataset,
- select the winning algorithm using the same estimates later reported as unbiased performance,
- use test-fold information for feature selection,
- perform preprocessing outside the training fold.

All preprocessing must live inside the appropriate pipeline.

---

# 6. Use Repeated Nested Cross-Validation

For datasets large enough to support it, use repeated nested CV.

Report distributions rather than only one mean.

For every model/target/disease report:

- \(R^2\),
- MAE,
- RMSE if useful,
- median,
- mean,
- standard deviation,
- confidence intervals.

If a disease cohort is too small for a reliable ML analysis, explicitly label it as underpowered instead of forcing a model.

---

# 7. Introduce Mandatory Baselines

Do not claim that topological indices are informative simply because a model obtains a large \(R^2\).

Topological indices are often strongly related to molecular size.

For each property compare at least:

### Null model

\[
M_0 : \hat y = \bar y_{\text{train}}
\]

### Molecular-size baseline

Use only very simple size descriptors such as:

- heavy-atom count,
- bond count,
- possibly carbon count,
- ring count where justified.

Avoid accidentally constructing a baseline that already contains the target.

### Topological-index model

\[
M_{\mathrm{TI}}:
y=f(T_1,\ldots,T_p).
\]

### Conventional molecular-descriptor baseline

Where scientifically appropriate include a standard RDKit descriptor model.

The headline quantity should include:

\[
\Delta R^2_{\text{TI-size}}
=
R^2_{\text{TI}}
-
R^2_{\text{size}}.
\]

Also consider:

\[
\Delta MAE_{\text{TI-size}}
=
MAE_{\text{size}}
-
MAE_{\text{TI}}.
\]

Use paired comparisons because predictions arise from the same outer folds.

---

# 8. Main Scientific Question

Do not frame the manuscript around:

> "TopoQSPR calculates 135 topological indices."

That is not a strong 2026 novelty claim.

Instead test the research question:

> **How sensitive are disease-specific topological-index QSPR conclusions to drug-set definition, molecular-size confounding, and evaluation methodology?**

This should become the central scientific story if the results support it.

---

# 9. Run Three Core Experiments

## Experiment A — Cohort-definition sensitivity

Compare:

### A1. Naive retrieval

The current broad/substr-based disease retrieval.

### A2. Controlled cohort

Controlled disease concept + approved/phase threshold + parent/salt normalization + deterministic deduplication.

Measure cohort overlap:

\[
J(A,B)=\frac{|A\cap B|}{|A\cup B|}.
\]

Also measure how much:

- sample size,
- property distributions,
- selected topological descriptors,
- correlation conclusions,
- model performance

change when cohort construction changes.

This can demonstrate that defining "drugs for disease X" is itself a nontrivial reproducibility problem.

---

## Experiment B — Molecular-size confounding

For every disease and target compare:

\[
R^2_{\mathrm{null}},
\qquad
R^2_{\mathrm{size}},
\qquad
R^2_{\mathrm{TI}},
\qquad
R^2_{\mathrm{standard}}.
\]

Then calculate:

\[
\Delta R^2_{\mathrm{TI-size}}.
\]

Determine which topological-index families still provide signal beyond molecular size.

Do not hide negative results.

If complex graph descriptors fail to outperform simple size variables, that is an important scientific finding.

---

## Experiment C — Evaluation optimism

For selected disease/target combinations compare:

\[
R^2_{\text{training}},
\]

\[
R^2_{\text{ordinary CV}},
\]

and

\[
R^2_{\text{nested CV}}.
\]

Quantify the optimism:

\[
B_{\text{train}}
=
R^2_{\text{training}}
-
R^2_{\text{nested}}
\]

and, if relevant,

\[
B_{\text{CV}}
=
R^2_{\text{ordinary CV}}
-
R^2_{\text{nested}}.
\]

This could become one of the strongest figures in the manuscript.

---

# 10. Reproduce Published Disease-Specific QSPR Studies

Find several representative published studies using graph-theoretic topological indices for small disease-specific drug sets.

Prefer studies for which:

- the compound list is available,
- formulas are clear,
- target properties are available,
- the analysis can be reconstructed.

For each study:

1. reproduce its dataset as faithfully as possible,
2. reproduce its reported regression,
3. verify whether the reported result is reproducible,
4. recompute using validated TopoQSPR descriptors,
5. evaluate using proper nested CV,
6. compare against the molecular-size baseline.

For example, compare:

\[
R^2_{\mathrm{reported}}
\]

against

\[
R^2_{\mathrm{nested}}
\]

and

\[
R^2_{\mathrm{size}}.
\]

The goal is not to attack previous authors.

The goal is to measure how methodological choices affect reported QSPR strength.

---

# 11. Multiple Testing and Correlation Analysis

The current project calculates many descriptor-target correlations.

Do not report raw \(p<0.05\) stars across hundreds of tests without correction.

Use a justified multiple-testing method, such as Benjamini-Hochberg FDR.

For correlations report:

- effect size,
- confidence interval,
- adjusted \(q\)-value,
- sample size.

Do not make the manuscript depend primarily on significance stars.

Also quantify descriptor collinearity.

Consider:

- correlation clustering,
- rank deficiency,
- VIF carefully,
- PCA only if justified,
- family-level summaries.

---

# 12. SHAP and Interpretability

Do not compute SHAP on a model fitted using all available data and then interpret it as if it were externally validated.

Separate:

- predictive evaluation,
- model interpretation.

If SHAP is retained:

- explain exactly which fitted model is being interpreted,
- avoid using the test data to select features,
- summarize strongly correlated descriptors by families,
- perform stability analysis across outer folds.

If SHAP is too unstable because of collinearity and small \(n\), say so.

---

# 13. Remove or Reclassify Unsupported ADMET Claims

Audit every "ADMET" field.

Any heuristic quantity such as a manually constructed Caco-2 formula must not be presented as an experimentally validated ADMET predictor.

Either:

- remove it,
- clearly relabel it as a heuristic,
- or replace it with a validated external model whose licensing and methodology permit use.

Likewise audit all bioactivity values for:

- target identity,
- assay type,
- units,
- relation operators,
- standardization,
- duplicate measurements.

Do not take the "first available bioactivity value."

---

# 14. Compare TopoQSPR Against Existing Tools

Perform an updated literature/software comparison, including at minimum:

- ChemGraphX,
- Mordred,
- RDKit descriptors,
- PaDEL-Descriptor if relevant,
- DOPtools,
- other current QSPR/QSAR workflow systems.

## ChemGraphX

Read the actual 2026 ChemGraphX paper.

Do not describe it inaccurately.

Determine precisely:

- what inputs it accepts,
- how many and what types of indices it computes,
- whether it performs QSPR modeling,
- whether it performs disease-to-drug cohort construction,
- whether it performs controlled ChEMBL indication mapping,
- whether it performs nested CV,
- whether it constructs provenance-tracked drug manifests.

The likely distinction to test is:

> ChemGraphX focuses on graph-descriptor computation, whereas TopoQSPR aims to construct an end-to-end disease-anchored, provenance-tracked, statistically defensible QSPR workflow.

Only make this claim if the literature supports it.

Cite ChemGraphX prominently.

Do not use "number of indices" as the principal novelty claim.

---

# 15. Decide the Appropriate Paper Type

After the audit and experiments, give three separate verdicts.

## A. SoftwareX-style software paper

Assess whether TopoQSPR is sufficiently:

- useful,
- documented,
- tested,
- reproducible,
- versioned,
- novel as scientific software.

## B. Journal of Cheminformatics-style methodological paper

Assess whether the cohort-sensitivity + size-confounding + evaluation-bias results are scientifically strong enough for a methodological cheminformatics paper.

## C. Lower-risk software venue

Give fallback options if the contribution is primarily software rather than methodology.

For each venue provide:

- fit,
- strengths,
- weaknesses,
- missing experiments,
- realistic acceptance risk.

Do not recommend a prestigious venue merely to be encouraging.

---

# 16. Software Engineering Required for v1.0

Before generating final paper numbers, implement or verify:

- proper Python package structure,
- `pyproject.toml`,
- pinned dependencies,
- deterministic seeds,
- pytest test suite,
- CI,
- descriptor unit tests,
- integration tests,
- deterministic manifests,
- explicit ChEMBL version,
- run configuration file,
- machine-readable provenance,
- Docker environment,
- documented CLI,
- documented API,
- semantic versioning,
- `v1.0.0` release,
- Zenodo archive/DOI if practical.

Every result used in the manuscript must be generated from a tagged commit/release.

Record:

- git commit,
- software version,
- database version,
- dependency versions,
- random seeds,
- configuration,
- cohort manifest.

---

# 17. Recommended Validation Disease Panel

Do not blindly use the disease list currently in the manuscript.

Select diseases based on:

- adequate number of approved small-molecule drugs,
- clean controlled ontology mapping,
- diverse chemical classes,
- different cohort sizes.

Candidate examples include:

- type 2 diabetes mellitus,
- hypertension,
- asthma,
- Alzheimer disease,
- breast neoplasms,
- epilepsy,

but verify whether each produces a statistically useful final cohort.

If a cohort is too small, replace it.

Record the selection rule before looking at model performance.

---

# 18. Avoid Target Leakage

Audit every target carefully.

Do not make headline claims for properties that are nearly deterministic functions of simple atom counts or that are computed from the same molecular representation in a way that makes prediction trivial.

Classify targets into:

1. nearly deterministic structural quantities,
2. physicochemical properties,
3. genuinely external experimental/biological targets.

Clearly distinguish them in the manuscript.

Predicting molecular weight from size-related graph invariants is mainly a sanity check, not a major scientific result.

---

# 19. Statistical Reporting

For all main comparisons:

- use effect sizes,
- confidence intervals,
- paired fold-level comparisons,
- corrected multiple testing when appropriate.

Where possible estimate uncertainty for

\[
\Delta R^2_{\mathrm{TI-size}}.
\]

Consider bootstrap or repeated-CV intervals, but ensure resampling respects the experimental unit.

Avoid pseudo-replication.

---

# 20. Figures the Final Paper Should Aim For

Produce publication-quality figures such as:

### Figure 1
TopoQSPR workflow:

\[
\text{Disease}
\rightarrow
\text{Controlled concept}
\rightarrow
\text{Drug manifest}
\rightarrow
\text{Molecular graph}
\rightarrow
\text{Verified indices}
\rightarrow
\text{Nested QSPR}
\rightarrow
\text{Reproducible report}.
\]

### Figure 2
Naive vs controlled disease cohort sizes / overlap.

### Figure 3
Heatmap of:

\[
\Delta R^2_{\mathrm{TI-size}}
\]

across disease × target.

### Figure 4
Training vs ordinary CV vs nested-CV performance.

### Figure 5
Descriptor-family contributions beyond molecular size.

### Figure 6
Published-study reproduction:
reported performance vs TopoQSPR nested-CV performance.

Only keep figures whose results survive the validation.

---

# 21. Manuscript Positioning

Do not claim:

> "TopoQSPR is novel because it computes 135 topological indices."

Instead, if supported by the experiments, build the paper around:

> **TopoQSPR provides a reproducible disease-anchored QSPR workflow that controls cohort definition, descriptor correctness, molecular-size confounding, model-selection bias, and computational provenance.**

A possible stronger research-paper title is:

> **TopoQSPR: Reproducible Disease-Anchored QSPR Reveals the Effects of Cohort Definition, Molecular Size, and Evaluation Bias on Graph-Topological Models**

A safer software-paper title is:

> **TopoQSPR: A Provenance-Tracked Platform for Disease-Anchored QSPR Using Verified Graph-Theoretic Descriptors**

Do not choose the final title until the experiments are complete.

---

# 22. Preserve Negative Results

If the data show that:

- topological indices rarely beat molecular-size baselines,
- many indices are redundant,
- published high \(R^2\) values collapse under nested CV,
- disease cohort definition changes conclusions materially,

do not regard that as failure.

Those may be the most important results of the paper.

The software should help reveal when graph-topological QSPR is informative and when it is not.

---

# 23. Required Final Deliverables

At the end of this work produce:

## A. `AUDIT_TOPOQSPR.md`

Include:

- verified bugs,
- mathematical errors,
- methodological weaknesses,
- reproducibility problems,
- manuscript inconsistencies,
- severity,
- fixes,
- exact file/function references.

## B. `INDEX_VERIFICATION.md`

For every descriptor:

- definition,
- citation,
- implementation mapping,
- unit tests,
- verification status.

## C. `EXPERIMENT_PLAN.md`

Freeze the main experiment design before looking at final results.

Include:

- diseases,
- targets,
- baselines,
- models,
- CV procedure,
- metrics,
- statistical tests,
- figures.

## D. `LITERATURE_GAP_2026.md`

Compare TopoQSPR with ChemGraphX and other relevant systems.

State clearly what is and is not novel.

## E. `PAPER_VERDICT.md`

Give a direct verdict:

- ready to submit,
- almost ready,
- major experiments still needed,
- or not publishable in current form.

Provide reasons.

## F. Code fixes

Implement all critical/high-priority corrections.

Add tests.

Do not silently change mathematical definitions.

## G. Updated manuscript

Only after the verified v1.0 experiments finish, revise `paper/main.tex`.

Every quantitative claim must trace to a reproducible result file generated from the tagged code release.

---

# 24. Working Style

Work skeptically.

For every important claim use this pattern:

> **Claim → evidence → independent verification → verdict.**

Never accept:

- README statements,
- code comments,
- manuscript equations,
- previous AI conclusions,

as proof.

Check them.

If an external source is required, open the original paper/documentation whenever possible.

Prefer primary literature over blogs and summaries.

When you discover that an earlier assumption was wrong, explicitly correct it rather than defending it.

---

# 25. Final Research Question

At the end, determine whether the evidence supports the following central question:

> **When disease-specific drug cohorts are rigorously defined and model evaluation is leakage-free, how much predictive information do graph-theoretic topological descriptors provide beyond molecular size, and how sensitive are conventional QSPR conclusions to cohort construction and evaluation methodology?**

If this question is well answered with robust multi-disease experiments, TopoQSPR should be treated not merely as a calculator but as a methodological and reproducibility study.

If the evidence does not support this framing, say so and recommend the narrower SoftwareX version.

Do not overclaim novelty, statistical evidence, or publication readiness.
