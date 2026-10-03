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
`max_phase_for_ind` = 4 → small molecules → parent compound → ChEMBL structure-pipeline parent (salts/solvents removed); structures still containing >1 fragment, metals or no carbon are **excluded with a reason**, not reduced to a fragment →
dedup by full InChIKey of the standardised parent → sorted by ChEMBL ID → manifest.
The ChEMBL release is recorded; the study uses one fixed release (ChEMBL 37 unless a newer
release is out at freeze time).

### 3.2 Naive cohort (for Experiment A only)
The v0.x retrieval, frozen as a function (`mesh_heading__icontains`, all phases; its keyword
fallback is omitted because it never ran: the ChEMBL search endpoint returns 404), with **all** candidates kept (no random truncation) so that it is deterministic.

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

### 3.4 Open decision before freezing
The phase-4 filter alone still admits drugs attached by label text mining (e.g. lisinopril for
diabetes). The default evidence filter (none / exclude DailyMed-only / ATC class) is chosen from
the cohort module's measured counts **before** any modelling, and recorded here.

## 4. Targets

Classified before analysis; only class 2 and 3 can support headline claims.

| Class | Meaning | Targets |
|---|---|---|
| 1 Structural / near-deterministic | Sanity checks only | Molecular weight |
| 2 Computed physicochemical | Genre replication and Experiments A/C; **no headline claims** (these are atom-contribution sums, so modelling them is close to circular) | Crippen logP, XLogP3 (PubChem), TPSA, molar refractivity |
| 3 External experimental | **Headline analysis (Experiment B)** | Candidates: aqueous solubility (AqSolDB), melting point (Bradley open set), experimental logP (OPERA/PHYSPROP). Each source is used only after its licence is checked; matched to drugs by InChIKey |

Audit note (2026-10): all v0.x indices are blind to element type and bond order (hexane, 1-pentanol,
pentylamine and 1,4-butanediol get identical values), so polarity-type targets can only be
predicted through their correlation with size.

Duplicate targets are removed (RDKit and PubChem molecular weight are the same quantity).
Bioactivity values are excluded from v1.0 modelling.

## 5. Feature sets (models compared on identical folds)

| Model | Features |
|---|---|
| M0 null | none (training-fold mean) |
| M_size | heavy-atom count, bond count |
| M_size+ (sensitivity) | + ring count (for single-fragment molecules ring count = bonds − atoms + 1, so this only matters for non-linear learners) |
| M_comp | element counts (C, N, O, S, P, halogens, other) + H count |
| M_TI | all verified TIs (after removing exact duplicates) |
| M_TI-family | one model per index family (Experiment B, family contributions) |
| M_std | RDKit 2D descriptors, **excluding every descriptor in the target's derivation family** (e.g. for Crippen logP: MolLogP, MolMR, SlogP_VSA*, SMR_VSA*) |

## 6. Models and cross-validation

- Model families: the tool's regressor set. Hyperparameters **and** family are chosen in the
  inner loop only.
- Outer: 5-fold, repeated 5 times with fixed seeds 0–4. Inner: 3-fold.
- All preprocessing (scaling, missing-value handling, constant-feature removal) is inside the pipeline.
- Every feature set (including size and composition baselines) uses the **same** model families and
  inner tuning, so baseline comparisons are fair.
- Scaffold-grouped outer folds (Bemis–Murcko) as the third step of Experiment C. Acyclic molecules
  (empty scaffold) each form their own group; report the number of groups and the largest group's share.
- Protocol follows Ash et al. 2025 (JCIM 65:9398) for method comparison.

## 7. Endpoints and statistics

- **Primary (confirmatory):** paired ΔR²(TI − size) and ΔMAE(size − TI), for each class-3 target on the
  pooled approved-drug set, and for each main-panel disease × class-3 target with ≥ 50 measured values.
  Averaged over repeats; 95% CI via the corrected repeated k-fold t-test (Nadeau & Bengio 2003;
  Bouckaert & Frank 2004) on paired per-fold differences. TIs are called *informative beyond size*
  if the CI lower bound > 0 after Benjamini–Hochberg control (q = 0.05) over this confirmatory family.
  Sensitivity: paired bootstrap of the whole nested procedure on the pooled set (CV intervals can
  undercover; Bates, Hastie & Tibshirani 2023), and Benjamini–Yekutieli instead of BH (dependence).
  **Practical significance:** ΔR² ≥ 0.05 is the pre-specified minimal meaningful gain
  (to be confirmed by co-authors before freezing).
- **Secondary (exploratory):** ΔR²(TI+comp − comp); ΔR² per TI family; M_std vs M_TI; class-2 targets;
  R² and MAE of every model (mean, SD, median); n for every endpoint.
- **Size dependence of the descriptors themselves:** for every TI and TI family, the cross-validated R²
  of that TI regressed on the size features (how much of each index *is* molecular size).
- **Correlations:** Spearman ρ with 95% CI, BH-adjusted q, n; family-level summaries;
  collinearity reported via correlation clustering.

## 8. Experiments

- **A. Cohort sensitivity:** naive vs controlled per disease: Jaccard overlap of parent InChIKeys,
  n, property distributions, top-10 correlated indices overlap, ΔR²(TI − size).
- **B. Size and composition confounding:** baseline ladder M0 → M_size → M_comp → M_TI → M_TI+comp → M_std;
  ΔR²(TI − size), ΔR²(TI − comp), and ΔR²(TI+comp − comp) (what indices add on top of composition).
  Run primarily on the **pooled** set of all curated approved small-molecule drugs with experimental
  endpoints (per-disease approved cohorts are too small for stable estimates); per-disease results
  are secondary.
- **C. Evaluation optimism:** training R² → ordinary CV (v0.x protocol: tuned on all data, best of
  families) → repeated nested CV (random outer folds) → scaffold-grouped nested CV (Bemis–Murcko).
  Report each step's change: B_train = R²_train − R²_nested, B_CV = R²_CV − R²_nested,
  B_scaffold = R²_nested − R²_scaffold.
- **D. Replication:** 3–5 published disease-specific TI-QSPR studies. Selection is fixed before any
  re-analysis: a pre-written Europe PMC query (recorded here at freeze), screened in date order
  (newest first) against eligibility criteria: (i) compound list resolvable to structures,
  (ii) property values tabulated, (iii) index formulas stated, (iv) n ≥ 10. The first eligible studies
  are taken, subject to co-author agreement, with no regard to their reported R². Chain per study:
  published R² → reproduced in-sample R² → ordinary CV → nested CV → scaffold CV → size baseline,
  with CIs and the optimism gap.

## 9. Interpretation (SHAP)

Secondary and descriptive only. SHAP is computed for the model refit in each outer training fold,
evaluated on that fold's test set, summarised by index family, and its stability across folds
is reported. If it is unstable, we say so.

## 10. Confirmatory vs exploratory

Only the primary endpoint in §7 is confirmatory. Everything else (secondary endpoints, family-level
results, class-2 targets, SHAP, per-disease results below the n threshold) is labelled exploratory
in the paper.

## 11. GO / NO-GO checkpoint

The main experiment starts only after `docs/CHECKPOINT.md` reports GO. It lists: unresolved code or
definition issues; verified vs unverified indices; independent cross-checks (Mordred, ChemGraphX);
cohort validation; leakage audit; frozen endpoints, disease rule, baselines and statistics;
replication candidates; threats to internal and external validity; reproducibility status.

## 12. Provenance

Every result file records: git commit and tag, package version, ChEMBL release, dependency
lock, seeds, run configuration and cohort manifest. Paper numbers come only from a tagged release.

## Amendments

_None yet._
