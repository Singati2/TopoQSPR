# 01 — Data provenance and chemical identity audit (Gate B)

Detailed licence evidence: `01b_DATA_LICENCES.md` (primary pages read by a sub-agent; access
limits stated there). Census and identity statistics: `05_TARGET_BLIND_DATA_FEASIBILITY.md`.

| Item | Finding | Status | Evidence |
|---|---|---|---|
| ChEMBL release and licence | ChEMBL_37, CC BY-SA 3.0. Redistributed manifests and results must carry the licence plus the attribution "ChEMBL data is from http://www.ebi.ac.uk/chembl - the version of ChEMBL is chembl_37" | CONFIRMED | 01b |
| Bradley MP | CC0. Measured, aggregated: PHYSPROP 9,687 records, Alfa Aesar 8,785, Karthikeyan 4,083, academic sites 1,480, DrugBank 875, others. `donotuse` rows excluded (377) | CONFIRMED | Bradley xlsx `source` column |
| OPERA (PHYSPROP) logP / MP / WS | **Measured** literature values (curated PHYSPROP), not OPERA predictions. Package records "MIT (repository LICENSE)": that is the repo's *software* licence. No data licence in the zip. Upstream SRC/EPA PHYSPROP terms not found | **licence UNRESOLVED**; measurement type CONFIRMED | 01b §1 |
| AqSolDB | CC0 1.0, measured, 9 sources. Two-value conflicts were resolved toward the value closer to an ALOGPS prediction (a model-informed curation step). Shares EPI Suite/PHYSPROP provenance with OPERA WS | CONFIRMED (caveat) | 01b |
| AstraZeneca logD7.4 | CC BY-SA 3.0, measured shake-flask, range −1.5 to 4.5 (assay limits → censored) | CONFIRMED | 01b; `expdata.py` assay_range |
| Biogen | MIT in repo, measured in-house; not used in Experiments B–E | CONFIRMED | 01b |
| Identity rule | Values join drugs by **full** InChIKey of the standardised parent (not the first block). Stereo variants join separately and are grouped by first block in CV (33 groups in the drug cohort; 368 constitutions with stereo variants in Bradley) | VERIFIED | `expdata.py` docstring; 05 |
| Salts and parents | ChEMBL Structure Pipeline parent; salt forms collapse to parent (`keep_salt_forms=False`). **Risk:** MP and solid-state solubility are properties of the solid form (salt or polymorph), not of the parent. A parent-level join can pair a salt's MP with the free base's structure | **METHODOLOGICAL RISK** (HIGH for MP and logS interpretation) | design property; not quantified here |
| Replicates | Median of records. Parents whose replicate SD exceeds 10 °C (MP) or 0.5 log (logP/logS) are dropped (MP: 299 of 19,799) | VERIFIED; wording corrected (06b R1) | 05 |
| Censoring | ChEMBL `standard_relation` ≠ "=" and AZ assay limits → censored, excluded | VERIFIED | audit 2026-10-12 |
| Bradley ↔ OPERA independence | **Not independent at source level** (PHYSPROP inside Bradley). The external test is source-disjoint by construction (06b R3) | CONFIRMED risk, mitigated | 05 §4 |
| Crippen MolLogP in RDKit 2D for measured logP | Fitted to measured logP; possible overlap with PHYSPROP training data | METHODOLOGICAL RISK (F11) | audit F11 |
| Suspect value | Epirubicin MP 344.53 °C (DrugBank-sourced) vs doxorubicin 204.5 °C | UNVERIFIED (kept) | audit F10 |
| Selection bias | Approved drugs with a measured endpoint (959/1821) are not a random sample of chemistry. Bradley is dominated by catalogue chemicals | LIMITATION | — |
| Disease cohort semantics | ChEMBL `drug_indication` with MeSH heading and indication-specific `max_phase_for_ind`. The ChEMBL API exposes these filters directly, so the cohort builder is a convenience, not a new capability | VERIFIED (prior-art impact on novelty D) | 02 §D |

**Not done:** re-download and re-hash of sources (the cached copies' pinned SHA-256 were
trusted); an ID-level check of a random sample of ChEMBL MeSH indications; quantifying the
salt/free-base MP mismatch rate (would need record-level salt metadata, which Bradley lacks).

**Gate B data verdict:** the primary sources are genuine measurements, and identity handling is
explicit and conservative. Proceed, with the salt-form caveat for MP stated as a limitation.
