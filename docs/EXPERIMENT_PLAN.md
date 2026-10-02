# TopoQSPR v1.0 — experiment plan

**Status: DRAFT (not frozen).** It will be frozen, with the commit hash recorded below, after
(1) the October 2026 audit is merged and (2) the pilot run has debugged the pipeline.
After freezing, any change is logged in the *Amendments* section with date and reason.
No model performance from the main panel may be looked at before freezing.

Frozen at commit: `—`  Date: `—`

## 1. Question

When disease-specific drug cohorts are rigorously defined and evaluation is leakage-free,
how much predictive information do graph-theoretic topological indices (TIs) carry beyond
molecular size? How sensitive are conventional QSPR conclusions to cohort construction and
evaluation methodology?

Negative results (TIs not beating size; published R² collapsing) are reported as results.

## 2. Pilot (quarantined)

- One disease **outside** the candidate list is used only to debug the pipeline end to end.
  Proposed: *Migraine Disorders*; replace if its curated cohort is < 30.
- The pilot never sets thresholds, metrics, targets or the disease panel. Its results are not
  reported as findings.

## 3. Cohorts

### 3.1 Controlled cohort (v1.0 default)
Free text → user-confirmed MeSH heading (and EFO ID where present) → ChEMBL `drug_indication`
rows for that exact heading (no substring, no keyword fallback) → indication-specific
`max_phase_for_ind` = 4 → small molecules → parent compound → largest organic fragment →
dedup by full InChIKey of the standardised parent → sorted by ChEMBL ID → manifest.
The ChEMBL release is recorded; the study uses one fixed release (ChEMBL 37 unless a newer
release is out at freeze time).

### 3.2 Naive cohort (for Experiment A only)
The v0.x retrieval, frozen as a function (`mesh_heading__icontains`, all phases, keyword
fallback), with **all** candidates kept (no random truncation) so that it is deterministic.

### 3.3 Disease panel — selection rule (fixed before any modelling)
1. Candidate list (MeSH headings, written down before counting):
   Diabetes Mellitus, Type 2; Hypertension; Asthma; Breast Neoplasms; Epilepsy;
   Depressive Disorder, Major; Schizophrenia; Rheumatoid Arthritis; Hypercholesterolemia;
   Heart Failure; Pulmonary Disease, Chronic Obstructive; HIV Infections; Parkinson Disease;
   Alzheimer Disease; Prostatic Neoplasms; Psoriasis.
2. Run the controlled cohort for each, count curated parents (counts only; no properties, no models).
3. **Main panel:** every candidate with n ≥ 50 curated parents, up to 12, by descending n.
4. Candidates with 20 ≤ n < 50 are reported descriptively as *underpowered*; no ML claims.
5. Record the counts table here before step 4 of the analysis begins.

## 4. Targets

Classified before analysis; only class 2 and 3 can support headline claims.

| Class | Meaning | Targets |
|---|---|---|
| 1 Structural / near-deterministic | Sanity checks | Molecular weight |
| 2 Computed physicochemical | Main analysis | Crippen logP, XLogP3 (PubChem), TPSA, molar refractivity, PubChem complexity |
| 3 External experimental | Main analysis if obtainable at adequate n | `TBD from audit` (candidates: experimental logP, melting point, aqueous solubility) |

Duplicate targets are removed (RDKit and PubChem molecular weight are the same quantity).
Bioactivity values are excluded from v1.0 modelling.

## 5. Feature sets (models compared on identical folds)

| Model | Features |
|---|---|
| M0 null | none (training-fold mean) |
| M_size | heavy-atom count, bond count |
| M_size+ (sensitivity) | + ring count |
| M_TI | all verified TIs (after removing exact duplicates) |
| M_TI-family | one model per index family (Experiment B, family contributions) |
| M_std | RDKit 2D descriptors, **excluding any descriptor that defines or is the target** (e.g. MolLogP when predicting Crippen logP) |

## 6. Models and cross-validation

- Model families: the tool's regressor set. Hyperparameters **and** family are chosen in the
  inner loop only.
- Outer: 5-fold, repeated 5 times with fixed seeds 0–4. Inner: 3-fold.
- All preprocessing (scaling, missing-value handling, constant-feature removal) is inside the pipeline.
- Sensitivity: scaffold-grouped outer folds (Bemis–Murcko), to check analog leakage.

## 7. Endpoints and statistics

- **Primary:** ΔR²(TI − size) per disease × target, averaged over repeats. 95% CI via the
  corrected repeated k-fold t-test (Nadeau & Bengio 2003; Bouckaert & Frank 2004), using
  paired per-fold differences.
- TIs are called *informative beyond size* for a disease × target if the CI lower bound > 0
  after Benjamini–Hochberg control (q = 0.05) across all main-panel disease × target pairs.
- **Secondary:** ΔMAE(size − TI); R² and MAE of every model (mean, SD, median); M_std vs M_TI.
- **Correlations:** Spearman ρ with 95% CI, BH-adjusted q, n; family-level summaries;
  collinearity reported via correlation clustering.

## 8. Experiments

- **A. Cohort sensitivity:** naive vs controlled per disease: Jaccard overlap of parent InChIKeys,
  n, property distributions, top-10 correlated indices overlap, ΔR²(TI − size).
- **B. Size confounding:** R² of M0, M_size, M_TI, M_std; ΔR²(TI − size); family-level ΔR².
- **C. Evaluation optimism:** training R², v0.x-style CV (tuned on all data, best-of-family),
  nested R²; B_train = R²_train − R²_nested and B_CV = R²_CV − R²_nested.
- **D. Replication:** 3–5 published disease-specific TI-QSPR studies, chosen with all co-authors,
  by these criteria: compound list and property values published; index formulas clear.
  Report published R² vs reproduced in-sample R² vs nested R² vs size baseline, with CIs.

## 9. Interpretation (SHAP)

Secondary and descriptive only. SHAP is computed for the model refit in each outer training fold,
evaluated on that fold's test set, summarised by index family, and its stability across folds
is reported. If it is unstable, we say so.

## 10. Provenance

Every result file records: git commit and tag, package version, ChEMBL release, dependency
lock, seeds, run configuration and cohort manifest. Paper numbers come only from a tagged release.

## Amendments

_None yet._
