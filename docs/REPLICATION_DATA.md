# Experiment D: extracted replication data

Data for the replication arm of the experiment plan (EXPERIMENT_PLAN §8 D) are stored in `data/replication/`. Tooling is in `tools/extract_replication.py` and the tests are in `tests/test_replication_data.py`.
All ten studies are extracted: shortlist rows 1–5 are marked **proposed** and rows 6–10 **reserve**, in the order of
`docs/REPLICATION_SHORTLIST.md`, so the co-authors can choose. Nothing in this package fits a
model. The fit statistics below are copied from the papers, and no new performance number was computed.

## 1. Layout

```
data/replication/
  _pubchem_cache.json            every PubChem PUG-REST response used (keyed by request path)
  <study_key>/
    study.json                   citation, licence, table layout, properties (+ source label),
                                 index formulas, models with verbatim statistics, validation,
                                 consistency notes, resolution issues, reproducibility verdict,
                                 automated_checks (written by the tool)
    drugs.csv                    names as printed + resolution inputs (query_name, rule, note)
                                 + outputs (CID, isomeric SMILES, InChIKey, formula, MW, ...)
    properties.csv               long format: drug_id, row label, property, printed string,
                                 parsed value, missing_reason (NaN + reason, never 0)
    SOURCE.md                    attribution and quoted article licence
```

* Every value in `properties.csv` is re-extracted from the Europe PMC JATS XML (`fullTextXML`)
  by `extract_properties` using the table and column ids recorded in `study.json`. `extract`
  checks that the committed file still matches the XML: 7/7 stored studies match. The SHA-256
  of each XML as retrieved on 2026-10-03 is in `study.json` (`source.xml_sha256`), and so is the
  SHA-256 of the canonical extracted `properties.csv` (`source.properties_sha256`, written by
  `extract --record`, for all ten studies). `validate` checks that hash offline, together with
  row labels and row order (drug order, then property-column order), so an edited, swapped or
  relabelled value is caught without the XML. For the three licence-withheld studies `extract`
  compares its local re-extraction with the same hash.
* `equation_verbatim` is a substring of the rendered XML (`extract` checks this; `validate`
  rejects location remarks such as "(Table 20 caption)" in the field, which belong in `source`
  or `note`).
* Formulas are rendered from MathML (linearised: `(a)/(b)`, `sqrt()`, `_{}`, `^{}`) or from
  TeX (Springer `tex-math`). Formulas that the XML carries only as images are marked `[IMG …]`.
  For Kabeer et al., the index and property symbols that are images in the XML were read from
  the publisher's table pages, which carry the same CC BY 4.0 licence. Those URLs are listed in
  `source.other_sources`. Their equations are stored as in the XML (`equation_verbatim`, with
  `[IMG …]`), and the symbol-filled reading separately (`equation_interpreted`, with
  `interpretation_source`).
* `models[].statistics` hold the printed strings, unchanged (e.g. `"− 1.94"`).

### Commands

```
python tools/extract_replication.py validate              # offline: schema, re-resolution from cache, checks
python tools/extract_replication.py resolve [--offline]   # re-run PubChem resolution, refresh automated_checks
python tools/extract_replication.py fetch-xml [study]     # network: Europe PMC XML into ~/.cache/topoqspr/replication
python tools/extract_replication.py extract [--record] [study]
                                                          # compare properties.csv, its SHA-256 and the
                                                          #  equations with the XML (licence-withheld
                                                          #  studies: write properties.csv to the cache)
```

`resolve --offline` reproduces `drugs.csv` and `automated_checks` byte for byte from the
committed cache. `validate` fails if anything differs. An unknown study key is an error
(exit code 2) for every command.

## 2. Structure resolution rules

Each name is resolved with PubChem PUG-REST: name → CIDs → `SMILES` (isomeric), `ConnectivitySMILES`,
InChIKey, formula, MW, exact mass, heavy atoms, complexity, TPSA and XLogP3. The rule recorded
for each drug in `drugs.csv` is one of the following.

| rule | meaning | used for |
|---|---|---|
| `as_printed` | the printed name is queried unchanged | most drugs |
| `code_mapped` | the table prints a code (BD1, A1). The name comes from the paper's text | Çolakoğlu, Ahmed |
| `typo_corrected` | misspelt name. `query_name` holds the correction and `note` the evidence | 'riseridone' → risperidone (Table 2 MW 410.5 and HAC 30 match) |
| `parent_form` | the printed name denotes an ion or salt. `query_name` is the neutral parent, and `note` says whether the paper's formula or values describe the parent or the salt | 'Valproate' → valproic acid (formula C8H16O2 printed); 'meclofenamate' → meclofenamic acid (values differ from the acid's row, probably the salt) |
| `component_by_printed_index` | combination product whose printed index values identify one component. `query_name` is that component, and `note` gives the evidence | 'IMIPENEM-CILASTATIN' → imipenem (Table 7 EDS 15088 = imipenem's EDS; cilastatin gives 43610) |
| `unresolvable` | combination product, polymer or mixture with no identifying evidence. Status `unresolved`, never guessed | chondroitin sulfate |

`validate` refuses an `as_printed` row whose printed name is an anion ('…ate') but whose PubChem
record is the acid ('… Acid'); such a row must use `parent_form` with a note.

Status values: `resolved`, `ambiguous_first_cid` (several CIDs, first taken, flagged),
`resolved_multicomponent` (salt or mixture record) and `unresolved`. In this data set every
query returned exactly one CID and no record is multi-component. **209 of 210 rows resolve.**
The unresolved row is the polymer (chondroitin sulfate). Salt-form *property values*,
as opposed to names, are flagged in the notes, e.g. fluvoxamine maleate (Kabeer) and diclofenac
sodium (Pavithra).

`automated_checks` (recomputed by `validate`) records the following:
* n listed vs n reported;
* duplicate structures and shared InChIKey connectivity blocks, which must be kept in one CV fold
  (EXPERIMENT_PLAN §6);
* printed formula or CID vs PubChem;
* per property, the count of values and missing values;
* for MW, HAC, complexity, TPSA and logP, agreement with PubChem at the coarser precision: the
  more precise number, rounded to the coarser number of decimals, must give the other (a tie at
  exactly half a unit is accepted). Numbers at the same precision must be equal, so integers
  (HAC, integer-printed TPSA) are compared exactly. PubChem is read at its native precision
  because PUG-REST drops trailing zeros: XLogP3 and TPSA 1 dp (XLogP '3' is 3.0), MW as given.
  Complexity is computed by PubChem to 1 dp (PUG-View: valproic acid 93.4, phenelzine 77.3) but
  PUG-REST returns it truncated to an integer, so a printed decimal must truncate to the PUG-REST
  value and a printed integer must equal it. A flag is set when a printed "MW" is the
  monoisotopic mass;
* identical property rows.

### Source labels

`source_label` per property:
* `computed_pubchem`: confirmed by the automated comparison. `validate` requires at least 80 %
  (`COMPUTED_PUBCHEM_MIN_MATCH`) of the compared values to match.
* `mixed`: part of the column equals PubChem and part does not (salt forms, monoisotopic masses,
  mis-rounding or an unidentified source); `source_evidence` gives the split. Used for Ali MW
  (3/10), Kabeer CO (16/25) and Pavithra MW (7/10).
* `predicted_acd_chemspider`: ACD/Labs Percepta values as displayed by ChemSpider. Godlin & Radha
  state this source. For the other studies the label is inferred from the property set and from
  values identical to Godlin & Radha for shared drugs. ChemSpider itself was not queried because
  its API needs a key.
* `predicted_chemdraw`
* `measured`: melting points spot-checked against PubChem experimental values (PUG-View).
* `measured_or_predicted_unclear`

## 3. Summary

| # | key | status | n (listed / distinct structures) | properties | indices | model forms | licence | reproducible |
|---|---|---|---|---|---|---|---|---|
| 1 | `colakoglu2026_bipolar` | proposed | 18 / 18 | 9 (6 ACD, 3 PubChem) | 13 degree + 13 nbhd-degree-sum | 2-index ridge (LOOCV-selected), RF, XGBoost | CC BY 4.0 | approximately |
| 2 | `kabeer2026_antidepressants` | proposed | 25 / 25 (1 stereo pair) | 8 (4 ACD, 2 PubChem, CO mixed, MP measured) | 9 neighborhood | linear, quadratic, logarithmic (1 index) | CC BY 4.0 | exactly* |
| 3 | `vimala2026_anticancer` | proposed | 23 / 23 | 8 (7 ACD, MP measured) | 13 eccentricity | linear, ridge, elastic net (1 index) | CC BY 4.0 | exactly |
| 4 | `almulla2026_antifungal` | proposed | 20 / 20 | 9 ACD (3 modelled) | 9 neighborhood | linear, quadratic, cubic (1 index) | CC BY 4.0 | exactly* |
| 5 | `ali2026_stiff_person` | proposed | 10 / 10 | 5 (MW only has equations) | 10 degree | linear, quadratic, exponential, power | CC BY 4.0 | approximately |
| 6 | `pavithra2026_skin_cancer` | reserve | 10 / 10 | 8 | 4 chromatic | linear (1 index) | CC BY 4.0 | approximately |
| 7 | `tang2026_arthritis` | reserve | 50 / 48 + 1 polymer | 6 ChemDraw | 9 degree | linear, RF, XGBoost | **CC BY-NC-ND 4.0** | approximately |
| 8 | `godlin2026_mdd` | reserve | 21 / 21 (1 stereo pair) | 9 (6 ACD, 3 PubChem; sources stated) | 12 degree + 12 nbhd | linear, quadratic, logarithmic (1 index) | CC BY 4.0 | exactly |
| 9 | `ahmed2026_antiepileptic` | reserve | 20 / 19 | 10 ACD | 12 (family unclear) | ANN, fuzzy ANN | **CC BY-NC-ND 4.0** | not at all |
| 10 | `supriya2026_tuberculosis` | reserve | 13 / 13 (combination resolved to imipenem by its printed EDS) | 6 ACD (BP/E/FP n = 11) | EDS, DESI | linear, quadratic, log, exp, MLR | **CC BY-NC-ND 4.0** | approximately |

\* exactly in principle: the inputs and method are fully specified. Printed tables contain
internal contradictions, so not every printed number can be matched (details below).

Verdict definitions:
* **exactly**: the printed property and index tables plus the stated model specification determine
  the published statistics, and the indices are uniquely defined for the resolved structures.
* **approximately**: something is underspecified (hyperparameter grid, colouring, split membership,
  fitting method, unidentified model) or one or more drugs have no single structure.
* **not at all**: the reported results depend on unspecified model architecture or protocol.

Whether the printed numbers actually reproduce is the job of the Experiment D analysis and
was not checked here.

**Licence finding.** The shortlist assumed CC BY throughout, but rows 7, 9 and 10 are CC
BY-NC-ND 4.0. For these three studies the folders store names, resolution, table layout,
formulas and reported statistics, but **no property values**. `fetch-xml` + `extract` rebuild
the values locally in the user cache. Rows 1–6 and 8 are CC BY 4.0. The ACS article (row 1)
states only "licensed under CC-BY 4.0", which is quoted.

**Genre observations (descriptive).** All ten studies model mostly ACD/ChemSpider or ChemDraw
predictions and PubChem-computed descriptors; only the melting points in rows 2, 3 and 5 look
measured. MW, HAC or complexity are targets in five studies (rows 1, 2, 5, 6, 8),
and in Godlin & Radha `ReZ1^d` *is* the heavy-atom count (r = 1 is an identity). Stereo pairs
(citalopram/escitalopram, ibuprofen/dexibuprofen) and duplicated rows occur. Pavithra's chromatic
hyper-Zagreb values are 9× and 4× the edge count for 8 of 10 drugs.

## 4. Per study

### 1. Çolakoğlu & Benjamin 2026, *ACS Omega* 11:53404 — bipolar disorder (proposed)
- **n** = 18 (BD1–BD18, names in the text); all resolved. BD16 'riseridone' → risperidone (`typo_corrected`).
- **Properties** (Table 2): BP (°C at 760 mmHg), EoV (kJ mol–1), FP (°C), MR (cm3), POL (10–24 cm3), MV (cm3) → ACD/ChemSpider; MW (g mol–1), COM, HAC → PubChem computed (18/18 equal current PubChem; valproic acid COM 93.4 is PubChem's 1-decimal value, which PUG-REST truncates to 93). The paper's source statement is "obtained from PubChem and ChemSpider".
- **Indices**: 13 degree-based (Table 1: M1, M2, F, SDD, ISI, H, SC, R, σ, SO, ABC, ABS, AZ) and the same 13 with the degree replaced by the neighborhood degree sum (M1* … AZ*). Values are in Tables 5–6 and edge partitions in Tables 3–4.
- **Models**: exhaustive LOOCV search per property over 13 univariate linear models, 91 two-index ridge models, RF and XGBoost. The best (always ridge) is reported in Tables 9–10 (LOOCV) and the in-sample refit with its equation in Table 11.
- **Reported statistics (verbatim, Table 11, in-sample)**: BP R2 0.84898; EoV 0.88086; FP 0.87899; MR 0.97634; POL 0.98205; MV 0.96562; MW 0.98634; COM 0.96336; HAC 0.99984. LOOCV R2 (Tables 9/10) range from 0.80188 (BP, nbhd) to 0.99943 (HAC, degree), e.g. `BP = 5.5172(SDD) + 0.0722(AZ) + 159.021.`
- **Issues**: the ridge α grid and the RF/XGBoost grids are not in the text (Algorithm 1). No external validation.
- **Reproducible: approximately.** Index values and single-index OLS can be reproduced exactly. The ridge fits depend on the unstated α grid.

### 2. Kabeer, Mufti, Almehizia, Assiri & Ganati 2026, *Sci Rep* 16:25505 — SSRI/SNRI/other antidepressants (proposed)
- **n** = 25; all resolved. Citalopram/escitalopram share connectivity.
- **Properties** (Table 3): BP (°C), MV (cm³), MR (cm³), EV (kJ/mol) → ACD (inferred; identical to Godlin & Radha for shared drugs); MW, HAC → PubChem; CO → mixed (16/25 equal PubChem, drugs 1–9 differ); MP (°C) → measured, partly salt forms (fluoxetine 179 and sertraline 246 match PubChem experimental entries; venlafaxine 215 is not the base). No source is stated, and the abstract calls all eight "measured".
- **Indices** (Table 1, images in the XML): NM1, NM2, NABC, NF*, ND1–ND5 (neighborhood degree).
- **Models**: 9 indices × 8 properties × {linear, quadratic, logarithmic} (Tables 5–7). Best equations are in Tables 8–10 and LOOCV in Table 14.
- **Reported statistics (verbatim)**: linear BP `2.579([IMG d33e4555.gif])+105.111`, read as `2.579(ND_1(G))+105.111`, R2 0.773; quadratic BP R2 0.887 (headline); MR quadratic 0.878; LOOCV Q2 BP-quadratic 0.7109, MR-linear 0.867, MW-linear 0.7735.
- **Issues**: fluvoxamine MW 434.4 and HAC 30 are the maleate's; fluvoxamine MV printed 27.3; CO for drugs 1–9 differs from PubChem; citalopram and escitalopram have different index values despite identical graphs; Table 14 Q2 > R2 for MW and two index-symbol mismatches; Table 11 contradicts Tables 5–10.
- **Reproducible: exactly in principle** (all inputs printed, standard formulas, plain OLS). Some printed numbers cannot all hold at once.

### 3. Vimala & Ravi Sankar 2026, *Front Chem* 14:1860531 — anticancer drugs (proposed)
- **n** = 23; all resolved. All 23 printed formulas equal PubChem. 'Azacytidine' (Table 1) and 'Azacitidine' (Table 2) are the same drug; Actinomycin D = dactinomycin.
- **Properties** (Table 2, no units; "taken from Kour and Sankar (2024)"): BP, E, FP, MR, P, ST, MV → ACD (inferred); MP → measured (spot-checked: altretamine, busulfan, carmustine).
- **Indices** (Table 3): 13 eccentricity-based (ECI, AECI, TEI, IEI, AEI, EDS, CEI, EHI, EM1, EM2, EGA, EABC, MECI); values in Table 4.
- **Models**: best single index per property; linear, ridge (α = 1.0) and elastic net (α = 1.0, l1_ratio = 0.5); scikit-learn 1.3.0; LOOCV in Table 8.
- **Reported statistics (verbatim)**: linear R2 BP 0.9207 (EGA), E 0.9269, FP 0.8406 (MECI), MR 0.9841, P 0.9841, MV 0.9561, MP 0.2510 (EABC), ST 0.1404 (AECI). LOOCV Q2 (linear) BP 0.9039, MR 0.9777, ST 0.0299.
- **Issues**: typos in the ridge table (FP intercept, ST R2, MV slope "08.5142"); MP prints p < 0.01 where F = 7.0 on (1, 21) df gives p ≈ 0.015; cyclophosphamide and ifosfamide share all ACD values.
- **Reproducible: exactly.**

### 4. Almulla, Irfan, Mufti, Almulla & Ganati 2026, *Sci Rep* 16:25695 — "antifungal" drugs (proposed)
- **n** = 20; all resolved. Only 10 are antifungals. Levalbuterol = (−)-salbutamol (CID 123600).
- **Properties** (Table 3, source not stated, all ACD set): B.P (°C), D (g/cm3), E.O.V (kJ/MOL), F.P (°C), I.O.R, M.R (cm3), P (cm3), S.T (dyne/cm), M.V (cm3). Only MR, P and MV are modelled.
- **Indices** (Table 1): NM1, NM2, NABC, NF1, ND1–ND5 (neighborhood degree sum); values in Table 4.
- **Models**: linear, quadratic and cubic in one index (best ND4); Tables 6–8, in-sample only.
- **Reported statistics (verbatim)**: MR linear R2 0.914 (R 0.995), quadratic 0.921, cubic 0.922; P 0.913/0.920/0.922; MV 0.840/0.861/0.880. The abstract claims cubic R2 = 0.98.
- **Issues**: R vs R2 mismatch for MR; the quadratic P row repeats d in the e column; the cubic MV g has the wrong sign; Sres looks like an RSS.
- **Reproducible: exactly in principle.**

### 5. Ali, Ali, Malik & Shang 2026, *ChemistryOpen* 15:e70239 — stiff-person syndrome drugs (proposed)
- **n** = 10; all resolved. 'Valproate' → valproic acid (`parent_form`). The printed compound ID 389 is ornithine; the other 9 printed CIDs and all formulas match.
- **Properties** (Table 4, "from PubChem and ChemSpider"): M.W (label `mixed`: 3/10 PubChem average MW, 3 monoisotopic masses, diazepam 278.07 is neither), M.P (unclear), Density (ACD; gabapentin 2.60 implausible), M.V (ACD), logP (4/10 match XLogP3; diazepam 2.82 vs 3.0 does not).
- **Indices** (Table 1): M1, M2, mM2 (image), AZI (printed without the −2), ReZG1–3, HZ, F, SDD. **The printed SDD equals 2 × M2 for every drug.**
- **Models**: MW only, with linear, quadratic (Table 5), exponential and power (Table 6) fits (MATLAB). Table 7 gives R2/Q2 (LOOCV) per property without naming the model.
- **Reported statistics (verbatim)**: MW linear M1 R2 0.9795 (R 0.9595, inconsistent); quadratic M1 0.9754; exponential M1 0.9751. Table 7 Q2: MW 0.7794, Log P 0.3001, Density 0.4989, MP 0.3068, MV 0.1002.
- **Reproducible: approximately** (only MW has equations; exponential and power fitting method not stated; SDD miscomputed).

### 6. Pavithra & Ravi Sankar 2026, *Front Chem* 14:1848388 — skin cancer drugs (reserve)
- **n** = 10; all resolved.
- **Properties** (Table 1): complexity (9/10; diclofenac is the sodium salt's) and PSA (9/10; nicotinamide 57 vs TPSA 56.0) → PubChem; MW → mixed (7/10: diclofenac is the sodium salt's, nicotinamide 122.2 vs 122.12, methotrexate 454.5 vs 454.4); MR, P, MV, ST → ACD; logP matches no XLogP3 value.
- **Indices**: chromatic modified first/second Zagreb and chromatic first/second hyper Zagreb (equations 1–4). They depend on a proper colouring drawn only in Figure 2.
- **Models**: single-index linear fits, in-sample (SPSS/Excel), Tables 5–8.
- **Reported statistics (verbatim)**: HM1Φ: MW R2 0.916, MR 0.9057, P 0.906, C 0.829, MV 0.735; PSA, ST and logP < 0.21.
- **Reproducible: approximately** (the colouring is not unique). With the printed index values the fits are fully specified.

### 7. Tang, Noureen, Ashraf, Azam & Aslam 2026, *ACS Omega* 11:33216 — arthritis drugs (reserve; CC BY-NC-ND)
- **n** = 50 listed. Chondroitin sulfate is a polymer (unresolved); 'meclofenamate' → meclofenamic acid (`parent_form`; its printed values differ from the acid's row, probably the salt), so rows 28 and 39 share a CID; and ibuprofen/dexibuprofen share connectivity.
- **Properties** (Table 3, ChemDraw 19.0 predictions "cross-validated with PubChem", median when sources differ): BP, MP (stated °C, magnitudes fit kelvin), Pc (bar), MR (cm3/mol; missing for 2), tPSA (Å2), CLogP. **Values not stored (licence).** A local re-extraction shows whole-row copies (celecoxib = feprazone, oxaprozin = piketoprofen) and near-copies.
- **Indices** (Table 1): EM1, EM2, M1, M2, AZI, HM, F, SDD, ISZ. The printed ISZ formula equals the vertex count.
- **Models**: single-index linear (Tables 4–9), RF and XGBoost on all 9 indices (Tables 13, 16, full data), XGBoost 70/30 test split (Table 18) and bootstrap/CV (Table 19).
- **Reported statistics (verbatim)**: linear r2 MR (ISZ) 0.8903, BP (EM1) 0.7733; XGBoost full-data r2 0.9866 (BP) to 0.8857 (CLogP); **test split r2 BP 0.44279, MP 0.08008, CLogP 0.11195**; CV R2 BP −0.6311.
- **Reproducible: approximately** (unpublished stratified split; supporting information not in PMC).

### 8. Godlin & Radha 2026, *Front Chem* 14:1784285 — major depressive disorder (reserve)
- **n** = 21; all resolved; citalopram/escitalopram share connectivity and have identical property rows.
- **Properties** (Table 6, sources stated per property): MW, HAC, CO → PubChem computed (21/21 match); BP, EV, FP, MR, PO, MV → ACD/Labs Percepta v14 via ChemSpider.
- **Indices**: 12 functions ξ(l,m) (ABC, Randić, M1, M2, ReZ1, ReZ2, SO, SS, BMG, GBM, GA, mM2), each on degree and on neighborhood-degree edge partitions (24 indices; Tables 4–5).
- **Models**: linear, quadratic and logarithmic in one index; best per property in Table 13; LOOCV in Table 16; 2-drug external set (duloxetine, esketamine) in Table 17.
- **Reported statistics (verbatim)**: Table 13 r2: MW 0.963, HAC 1, CO 0.9362, BP 0.9455, EV 0.9352, FP 0.9062, MR 0.9559, PO 0.9564, MV 0.9101; Q2_LOO 0.8368–0.9442 (HAC 1).
- **Reproducible: exactly** (except that the two external drugs' property values are not printed).

### 9. Ahmed, Zaman, Abd-Elmonem, Siddig & Belay 2026, *Sci Rep* 16:21878 — antiepileptic drugs (reserve; CC BY-NC-ND)
- **n** = 20 codes, 19 drugs: lacosamide is A10 and A18. All resolved.
- **Properties** (Table 5, "Source: ChemSpider"): D, IR, MR, HBA, HBD, ACD/LogD, ACD/BCF, POL, ST, MV, all ACD. **Values not stored (licence).**
- **Indices** (Table 1): M1, M2, H, F, SS, ABC, RI, SC, GA, HZ, ReZ1, ReZ2. They are titled "eccentricity based", but M_u is undefined.
- **Models**: ANN and fuzzy ANN on all indices (Table 9); TOPSIS ranking (not a QSPR result).
- **Reported statistics (verbatim)**: ANN R2 0.7823 (HBD) to 0.9975 (MR); fuzzy ANN 0.7619 to 0.9944. The evaluation protocol is not stated.
- **Reproducible: not at all** (architecture, epochs, seeds and fold use unspecified). The drug list could still serve the baseline/CV chain (n = 19).

### 10. Supriya & Iyer 2026, *Sci Rep* 16:20043 — tuberculosis subset (reserve; CC BY-NC-ND)
- **n** = 13 tuberculosis drugs (the SARS-CoV-2 subset, n = 8, is out of scope). Imipenem-cilastatin is resolved to imipenem (`component_by_printed_index`): the printed Table 7 EDS 15088 equals imipenem's eccentric distance sum (cilastatin: 43610), and the same computation reproduces the printed EDS of ethambutol, isoniazid, pyrazinamide and ethionamide. BP, E and FP are missing for isoniazid and terizidone.
- **Properties** (Table 8, "obtained from Chem Spider", called experimental): BP, E, FP, MR, P, MV, all ACD. **Values not stored (licence).**
- **Indices**: EDS (eccentric distance sum) and DESI (detour eccentric sum, which needs longest-path distances); values in Table 7. Ethambutol (a tree) has DESI = EDS, as it should.
- **Models**: logarithmic, exponential, polynomial (order 1–3) and two-index MLR; best per property in Table 13; LOOCV in Table 16.
- **Reported statistics (verbatim, tuberculosis columns)**: R2 BP 0.85 (P2, EDS), E 0.82 (MLR), FP 0.84, MR 0.94 (P1, DESI), P 0.94, MV 0.92. LOOCV Q2_val: BP **− 1.94**, FP **− 1.803**, MR 0.919, MV 0.72.
- **Reproducible: approximately** (combination product: the indices describe imipenem, but which component the property values describe is not stated; log/exp fitting method not stated; coefficient tables contradict their captions).

## 5. Open points for the co-authors

1. Confirm rows 1–5. Each can be re-analysed, row 5 only partly (MW equations only).
2. Rows 7, 9 and 10 cannot be redistributed with values. If a reserve is needed, row 6 or 8 (both CC BY) avoids this.
3. Agree how the replication handles printed errors (e.g. Ali et al.'s SDD = 2·M2): keep the published index values for "published → reproduced", and our recomputed indices for the CV chain.
4. Stereo pairs and duplicates (Kabeer, Godlin, Tang, Ahmed) must be grouped in all CV (plan §6).
