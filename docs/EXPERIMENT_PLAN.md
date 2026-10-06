# TopoQSPR v1.0 — experiment plan

**Status: FROZEN (2026-10-05).** All checkpoint decisions are recorded in `docs/CHECKPOINT.md`
(decisions 1–2 by Advik, 3–5 by the lead author on delegation from Singati2). The freeze is the git tag
`plan-freeze-v1`. OSF registration is pending (needs the authors' OSF account); until then the tag is
the record. Any later change goes in *Amendments* with date and reason. No main-panel or pooled
measured-target model result had been computed or seen before this freeze.

Study configurations (config_sha256 of the study part):
- `studies/main_confirmatory.json` 63cce0fea1968a8f… (Experiment B, confirmatory family + controls)
- `studies/main_scaffold.json` 0dbd6bc282481ff6… (Experiment C scaffold step)
- `studies/exploratory_families.json` 99687e11396d5fee… (exploratory TI-family ladder)

Frozen at tag: `plan-freeze-v1`  Date: 2026-10-05

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

### 3.3a Panel rule versus measured coverage (counts only, 2026-10-03)
Under §3.3 only Hypertension has ≥ 50 curated approved parents (85); every other candidate has 9–45.
No disease × measured-target cell reaches 50 (largest: Hypertension × melting point, 47).
So the per-disease *modelling* panel is effectively empty for class-3 targets. **Open decision
before freezing** (no model results were seen): either (a) keep per-disease analyses descriptive
(Experiment A, cohort sensitivity) and model only the pooled set, or (b) lower the per-disease
threshold to n ≥ 30 for exploratory analyses only. **Decided 2026-10-03 (Advik): (a).**

### 3.4 Open decision before freezing
The phase-4 filter alone still admits drugs attached by label text mining (e.g. lisinopril for
diabetes). The default evidence filter (none / exclude DailyMed-only / ATC class) is chosen from
the cohort module's measured counts **before** any modelling, and recorded here.
**Decided 2026-10-04:** exclude indications supported only by DailyMed label text; "none" and
"FDA/EMA/ATC only" are Experiment A sensitivity analyses.

## 4. Targets

Classified before analysis; only class 2 and 3 can support headline claims.

| Class | Meaning | Targets |
|---|---|---|
| 1 Structural / near-deterministic | Sanity checks only | Molecular weight |
| 2 Computed physicochemical | Genre replication and Experiments A/C; **no headline claims** (these are atom-contribution sums, so modelling them is close to circular) | Crippen logP, XLogP3 (PubChem), TPSA, molar refractivity |
| 3 External experimental | **Headline analysis (Experiment B)** | Licence-checked sources (docs/EXPERIMENTAL_DATA.md). Pooled approved-drug coverage (1,821 curated parents, full-InChIKey match): melting point Bradley 711 / OPERA 509, logP OPERA 505, logD7.4 AstraZeneca 369, aqueous solubility AqSolDB 351 / OPERA 349. All ≥ 200. **Confirmatory set (decided 2026-10-03, Advik):** mp_bradley, logp_opera, logd74_az, logs_aqsoldb (one source per endpoint; the OPERA MP/WS sets as sensitivity analyses). |

Audit note (2026-10): all v0.x indices are blind to element type and bond order (hexane, 1-pentanol,
pentylamine and 1,4-butanediol get identical values), so polarity-type targets can only be
predicted through their correlation with size.

Target names must be unique across providers: the pipeline rejects a duplicated name rather than
collapsing it. v1.0 ships no PubChem target provider, so the RDKit/PubChem molecular-weight duplicate
of v0.x cannot arise.
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

The baselines (M0, M_size, M_comp) are part of every run: config validation rejects a feature-set
list without `null`, `size` and `composition`, or a reference other than `size`, unless
`evaluation.allow_custom_baselines` is set (recorded in the run record; not used for the frozen
analyses). See PIPELINE.md §2.

## 6. Models and cross-validation

- Model families: the tool's regressor set. Hyperparameters **and** family are chosen in the
  inner loop only.
- Outer: 5-fold, repeated 5 times with fixed seeds 0–4. Inner: 3-fold.
- All preprocessing (scaling, missing-value handling, constant-feature removal) is inside the pipeline.
- Every feature set (including size and composition baselines) uses the **same** model families and
  inner tuning, so baseline comparisons are fair.
- **Stereoisomer grouping (all CV, including random folds):** drugs sharing the InChIKey connectivity
  block (first 14 characters) are always in the same fold. Stereoisomers have identical hydrogen-suppressed
  graphs and therefore identical index, size and composition vectors (e.g. citalopram/escitalopram);
  splitting them would leak (Kapoor & Narayanan 2023). The scaffold-grouped variant groups by scaffold,
  which already contains this grouping.
- Scaffold-grouped outer folds (Bemis–Murcko) as the third step of Experiment C. Acyclic molecules
  (empty scaffold) each form their own group; report the number of groups and the largest group's share.
- Protocol follows Ash et al. 2025 (JCIM 65:9398) for method comparison.

## 7. Endpoints and statistics

- **Primary (confirmatory):** paired ΔR²(TI − size) and ΔMAE(size − TI), for each class-3 target on the
  pooled approved-drug set, and for each main-panel disease × class-3 target with ≥ 50 measured values.
  Averaged over repeats; 95% CI via the corrected repeated k-fold t-test (Nadeau & Bengio 2003;
  Bouckaert & Frank 2004) on paired per-fold differences. TIs are called *informative beyond size*
  for an endpoint if its two-sided corrected-t p-value survives Benjamini–Hochberg control at q = 0.05
  over the confirmatory family **and** the mean ΔR² > 0 (as implemented in `confirmatory.csv`). The
  family size m is fixed in advance (missing p-values count as non-rejections, not as a smaller m), and
  underpowered endpoints (n < 50, or smallest outer test fold < 10 rows; EVALUATION_PROTOCOL §11.2)
  are excluded from the family before testing.
  Sensitivity: paired bootstrap of the whole nested procedure, on the pooled set only (at per-disease
  n the resampled test folds are too small for R² to be defined; pilot: 17–19 of 19 resamples undefined)
  (CV intervals can undercover; Bates, Hastie & Tibshirani 2023), and Benjamini–Yekutieli instead of BH (dependence).
  Frozen per-disease configs set `protocol.bootstrap.enabled = false` unless the code default is already
  pooled-only (PIPELINE.md §9).
  **Practical significance:** ΔR² ≥ 0.05 is the pre-specified minimal meaningful gain
  (decided 2026-10-04).
- **Secondary (exploratory):** ΔR²(TI+comp − comp); ΔR² per TI family; M_std vs M_TI; class-2 targets;
  R² and MAE of every model (mean, SD, median); n for every endpoint.
- **Y-randomised control:** for every feature set in the confirmatory family, the full nested pipeline on
  200 target permutations (Rücker et al. 2007); a feature set is reported as uninformative if its real
  R² does not exceed the 95th percentile of the permutation distribution.
- **Applicability domain:** kNN-distance AD fitted in each outer training fold; the share of outer-test
  molecules outside the domain is reported with every result; primary metrics use all molecules, and
  in-domain-only metrics are secondary.
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
  re-analysis: the Europe PMC query fixed in `tools/replication_screen.py` (run 2026-10-03: 239 open-access candidates, `docs/replication_candidates.csv`), screened in date order
  (newest first) against eligibility criteria: (i) compound list resolvable to structures,
  (ii) property values tabulated, (iii) index formulas stated, (iv) n ≥ 10. Selection (decided 2026-10-04, before any re-analysis): the eligible studies that are exactly
  reproducible from their published data (docs/REPLICATION_DATA.md rows 2, 3, 4, 8), with no regard to
  their reported R². Chain per study:
  published R² → reproduced in-sample R² → ordinary CV → nested CV → scaffold CV → size baseline,
  with CIs and the optimism gap. For comparability with the genre, Q²_F3 and the concordance correlation
  coefficient are also reported (secondary). A permuted-target and a descriptor-pool-randomised control
  (pick the best of the same number of random indices) quantify chance correlation at the study's n.

## 9. Interpretation (attribution deferred)

**Not part of v1.0** (pre-freeze change, 2026-10-03; see *Pre-freeze changes*). The v1.0 package has
no SHAP or other attribution code, so no attribution analysis is run or reported. Fold-internal
attribution (computed on the model refit in each outer training fold, evaluated on that fold's
test set, summarised by index family, with stability across folds) is deferred to a later release
and, if added, will be exploratory only.

## 10. Confirmatory vs exploratory

Only the primary endpoint in §7 is confirmatory. Everything else (secondary endpoints, family-level
results, class-2 targets, per-disease results below the n threshold) is labelled exploratory
in the paper.

## 11. GO / NO-GO checkpoint

After GO, the frozen plan is registered on OSF Registries (secondary-data template) before any
main-panel result is seen; the registration DOI and the frozen commit hash are recorded in the header.


The main experiment starts only after `docs/CHECKPOINT.md` reports GO. It lists: unresolved code or
definition issues; verified vs unverified indices; independent cross-checks (Mordred, ChemGraphX);
cohort validation; leakage audit; frozen endpoints, disease rule, baselines and statistics;
replication candidates; threats to internal and external validity; reproducibility status.

## 12. Provenance

Every result file records: git commit and tag, package version, ChEMBL release, dependency
lock, seeds, run configuration and cohort manifest. Paper numbers come only from a tagged release.

## Pre-freeze changes (logged for transparency; the plan is not yet frozen)

- 2026-10-03, §9: SHAP removed from v1.0; attribution deferred (no attribution code in the v1.0
  package). §10 no longer lists SHAP.
- 2026-10-03, §4: "duplicate targets are removed" replaced by the implemented rule (duplicate target
  names are rejected).
- 2026-10-03, §7: family membership also excludes members whose smallest outer test fold has < 10
  rows (as implemented); bootstrap restricted to the pooled set in the configs.
- 2026-10-03: the paper's "additional annotations" (PubChem properties, UniChem, Lipinski/Veber/PAINS
  filters) are v0.x features that are not in the v1.0 package and are not part of this plan.

## Amendments

- **2026-10-05 (before any Experiment D re-analysis), replication set.** The four studies selected
  on 2026-10-04 were chosen for exact reproducibility of their *data*; one of them, Vimala & Ravi
  Sankar 2026 (shortlist row 3), models only *eccentricity-based* indices (ECI, AECI, EABC, EGA,
  MECI), a family TopoQSPR does not implement, so none of its 48 fits can be recomputed. It is
  replaced by the next shortlist study whose modelled indices TopoQSPR implements, in date order,
  under the same rule; if none qualifies, Experiment D runs on three studies. The converter
  `tools/replication_to_spec.py` refuses any study with an unmapped modelled predictor, so this
  check is mechanical. Two degree indices used by Godlin & Radha 2026 (ReZ1, mM2) were added to the
  registry with closed-form tests; this adds indices to the `TI` feature set used by Experiments B
  and C? No: those runs were launched from tag `plan-freeze-v1` (102 indices) and are unaffected;
  the 104-index registry applies to Experiment D only. Model forms the harness cannot reproduce
  (logarithmic, ridge, elastic net, RF, XGBoost) are listed with their reported statistics and not
  compared.
- **2026-10-05 (before any Experiment D re-analysis), printed index values.** Step (i) of
  Experiment D compares the studies' printed index values with the recomputed ones. The three
  studies print them (Almulla Table 4, Kabeer Table 4, Godlin Tables 4–5); they are now extracted
  from the archived XML like the property values (`index_values.csv`, hashed in `study.json`) and
  passed to the harness. Two converter defects found while doing so were fixed before any
  re-analysis was read: a neighbourhood label `X^nd` with no `S_X` in the registry fell back to the
  degree index `X` (Godlin's `M1^nd` would have been compared with `M1`); and a paper's LOO
  statistics (Kabeer Table 14, Godlin Table 16) were read as in-sample fits. LOO statistics are
  now kept beside the harness's own LOO Q² and not compared as R². Kabeer's printed thousands
  separators ("1,004") are parsed. The three runs were repeated from the corrected specs; the
  earlier outputs were discarded unread except for the index-value check that exposed the defects.
- **2026-10-05 (after the confirmatory run; before any amended control was computed), control learner.**
  The frozen protocol stage (CHECKPOINT decision 4, a budget choice) ran the Y-randomisation and the
  whole-procedure bootstrap with ridge as the control learner, while the confirmatory comparison
  uses the `selected` learner. The run shows that ridge on the TI set is dominated by three cardiac
  glycosides (digoxin, digitoxin, acetyldigitoxin; 54–69 heavy atoms) that the pre-specified kNN
  applicability domain flags as outside the domain in every fold: a linear model extrapolates on
  them to predictions of −1 000 log units, so the pooled R² of ridge-TI is −1.27 (melting point)
  and −651 (solubility) although its in-domain R² is ordinary. Consequences: the permutation
  control declares TI "uninformative" for two targets for a reason unrelated to information, and
  50–79 % of bootstrap resamples explode, so the percentile interval is meaningless. The ridge
  results are reported unchanged as the frozen control, with this explanation. The confirmatory
  statistics themselves (learner `selected`, corrected t, BH/BY) are unaffected and are not re-run.
  Amended sensitivity: the protocol stage is repeated with SVR (RBF; bounded predictions; the
  family `selected` chooses most often for both TI and size) as the control learner, K = B = 200,
  seeds and everything else identical (`studies/main_confirmatory_protocol_svr.json`). Running the
  controls with the `selected` learner itself is out of reach (one nested CV of the six families
  takes 18 min on this machine; 1 600 draws per stage). Both control runs are reported.
- **2026-10-06, amended control run completed.** `main_confirmatory_protocol_svr` (SVR control
  learner; 5 h 54 min). No bootstrap resample explodes (0 of 200 per target); the TI set beats its
  Y-randomised null for all four endpoints; bootstrap intervals of ΔR²(TI − size): melting point
  [+0.02, +0.14], log P [+0.05, +0.19], log D [−0.03, +0.09], log S [−0.03, +0.06]. The run used the
  released registry (104 indices) rather than the 102 of the frozen tag because it was launched from
  the working tree; its CV stage differs from the frozen run by at most 0.005 in ΔR². Both runs are
  archived in `results/expB` and both appear in the paper's controls table.
