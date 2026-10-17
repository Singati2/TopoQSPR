# 15. Confirmatory replication: feasibility (doc 12 P1)

**Date:** 2026-10-16.

**Status:** this is a plan plus a licence and power assessment. It records no census results and
reads no target values. Every number marked *estimate* comes from E3 scaling, not from the
candidate data.

## 1. Requirement

The replication needs a measured MP source with three properties:
- **Not used anywhere in this project.**
- **Records not shared with Bradley or OPERA**, because both are already "used". Bradley aggregates
  PHYSPROP (9,687 records), Alfa Aesar (8,785), Karthikeyan (4,083), academic sites (1,480) and
  DrugBank (875); see 01 §1. OPERA MP *is* PHYSPROP.
- **Licence terms that allow publishing derived tables.**

## 2. Candidates

| source | licence (checked) | provenance | overlap with what we used | verdict |
|---|---|---|---|---|
| **Tetko, Lowe & Williams 2016, patents MP** (figshare 2007426 v1, `melting_point.sdf`, 638 MB) | **CC0**: figshare API `license.name = "CC0"`, checked 2026-10-16 | text-mined from US patent experimental sections 1976–2014 (J Cheminform 8:2, doi:10.1186/s13321-016-0113-y) | not in Bradley's source list; not PHYSPROP. Molecule-level overlap unknown until the census | **primary candidate** |
| Karthikeyan 2005 (≈ 4,400) | journal supplement | MDL MDDR-derived | **inside Bradley** (4,083 records) | reject |
| Bergström drug set | journal supplement | literature | inside Bradley's source list | reject |
| OCHEM MP collections | per-upload, often CC-BY | mixed; include the patent data and Bradley | duplicates the above | reject (use the original sources instead) |
| PHYSPROP / EPI Suite MP | terms unknown | PHYSPROP | = OPERA MP, already used | reject |

### Known risks of the patent set
- **Text-mining errors:** for example, a pyrolysis/decomposition value read as a melting point. The
  set has a separate pyrolysis file, but the split is the authors'.
- **Ranges and salts:** patents often report a range, and many values are for salts.
- **Not every value is new:** some patents quote literature values, so reference-level independence
  is **not shown**.
- **Label:** the replication may be called `source_disjoint_data`. It may not be called
  "independent data" (EVIDENCE_SCHEMA in ResearchOS; doc 12 P1).

## 3. Target-blind census (step 1; reads structures only)

1. **Download** figshare file `melting_point.sdf` and record its sha256.
2. **Parse molecule blocks only.** Do not read MP fields. The parser drops every SD property except
   the record id.
3. **Standardise** with the package pipeline (parents, `chembl_structure_pipeline`), as for Bradley.
4. **Remove overlap:** drop any parent whose InChIKey connectivity block (first 14 characters)
   occurs in Bradley (all rows, including `donotuse`) or in OPERA MP/logP. This is constitution
   level, stricter than the E3 external rule.
5. **Count** the remaining parents, formula groups with ≥ 5 constitutions, molecules in those
   groups, and heavy-atom range.
6. **Go/no-go on counts only** (rule frozen here):

| criterion | threshold | why |
|---|---|---|
| molecules in formula groups with ≥ 5 constitutions, after overlap removal | ≥ 3,400 | power, §4 |
| formula groups with ≥ 5 constitutions | ≥ 300 | E3 external train had 305 |
| fraction of parents with ≥ 1 range or salt flag, if available without reading values | reported, no threshold | data quality |

Only after GO are the MP fields read, with the cleaning rules frozen in the protocol (step 2) beforehand.

## 4. Power (*estimate*)

E3a gave SE(ΔR²) ≈ 0.0094 with 6,725 molecules (15 folds, corrected t, t₁₄ = 2.145). Assume SE
scales as 1/√n. That is an approximation: the noise in patent values will be larger than in Bradley.

- **Detecting ΔR² = 0.04** (half of E3a), two-sided α = 0.05, power 0.80, needs SE ≤ 0.04 / (2.145 + 0.868) ≈ 0.0133.
- That gives n ≥ 6,725 × (0.0094 / 0.0133)² ≈ **3,400** molecules in formula groups with ≥ 5 constitutions.
- **Inflation for noisier labels:** if label noise doubled the residual variance, the effect in R²
  units would shrink too. The census therefore reports n, and the protocol fixes a minimum.
  **Power will not be re-estimated from outcome values.**

## 5. Frozen design for the confirmation (to be signed off with the census result)

| item | value |
|---|---|
| primary contrast | `TIa_composition − composition`, E3b-type split (groups = formula), seeds 0–2 × 5 folds |
| minimum effect | ΔR² > 0, and the 95% CI lower bound > 0. The practical threshold is 0.02, the δ of doc 14 |
| references kept at full strength | RDKit 2D (asinh) and ECFP, each with and without TIa (secondary, BH over all secondary contrasts) |
| extra secondary | `edgetypes_composition − composition`, and full vs edge types with the doc 14 fixed learners |
| exclusions | ranges wider than 5 °C; salts not reducible to a parent; values flagged decomposition; replicate SD > 10 °C (E3 rule) |
| reporting | every outcome, including null or negative; absolute R² beside every ΔR² |

## 6. Feasibility judgment

- **The data is feasible:** CC0, large, from a different kind of source.
- **Not yet shown:**
  - that enough non-overlapping isomer groups remain;
  - that the noise level permits the power target.

  The census (§3) answers the first; label noise cannot be assessed until the target-blind gate is passed.
- **Cost:** the census takes about 1–2 h CPU (standardisation of about 250k structures). The
  confirmation itself would cost about the same as one E3 MP cell (≈ 17 h for 10 sets). It needs
  owner approval at each step.
- **If the census fails:** report no-go, do not substitute a weaker source, and state in the paper
  that no confirmatory replication was possible.

## 7. Owner decisions
- [ ] approve the census (download about 640 MB, about 1–2 h CPU, structures only)
- [ ] accept the go/no-go thresholds in §3
- [ ] accept the confirmation design in §5 (frozen only after the census and before reading values)
