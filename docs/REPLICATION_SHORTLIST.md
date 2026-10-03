# Experiment D: replication shortlist

Search fixed in `tools/replication_screen.py` (Europe PMC, open access, 2020–2026; run 2026-10-03; 239 candidates in `replication_candidates.csv`).
The 60 newest candidates were screened against the pre-specified criteria (EXPERIMENT_PLAN §8 D) by six screeners **blind to reported performance**; every 'eligible' verdict was then re-checked by an independent auditor who tried to refute it (21 → 19).

Per the plan, the **first eligible studies in date order** are taken, subject to co-author agreement. Proposed: rows 1–5; rows 6–10 are reserves in order.

| # | Date | DOI | Drug set | n | Properties modelled (source) | Arockiaraj co-author |
|---|---|---|---|---|---|---|
| 1 | 2026-08-27 | [10.1021/acsomega.6c06753](https://doi.org/10.1021/acsomega.6c06753) | bipolar disorder drugs | 18 | BP, EoV, FP, MR, POL, MV (ChemSpider ACD predicted); MW, complexity, HAC (PubChem computed). Text says 'obtained from PubChem and ChemSpider'. |  |
| 2 | 2026-08-07 | [10.1038/s41598-026-65586-2](https://doi.org/10.1038/s41598-026-65586-2) | antidepressants (SSRI/SNRI/TCA/MAOI/atypical) | 25 | BP, MV, MR, MW, HAC, complexity, EV, MP (Table 3). The source is not stated (the text cites only DrugBank for drug info). The pattern looks like ChemSpider ACD predicted plus PubChem computed, with MP probably measured. |  |
| 3 | 2026-07-08 | [10.3389/fchem.2026.1860531](https://doi.org/10.3389/fchem.2026.1860531) | anticancer drugs (also melatonin, simvastatin, nitroxoline) | 23 | BP, MP, E, FP, MR, P, ST, MV (Table 2). Taken from Kour and Sankar (2024); the values match the ChemSpider ACD predicted set, with MP likely experimental. |  |
| 4 | 2026-07-07 | [10.1038/s41598-026-61020-9](https://doi.org/10.1038/s41598-026-61020-9) | antifungal drugs (only 10/20 are antifungals; the rest include octacosanol, L-threonic acid, menthyl salicylate, amiloxate, emylcamate, fluindione, levalbuterol, arbutin, marbofloxacin, florfenicol) | 20 | BP, density, EV, FP, refractive index, MR, polarizability, ST, MV (Table 3). Source not stated; this is the ChemSpider ACD predicted property set. |  |
| 5 | 2026-06-01 | [10.1002/open.70239](https://doi.org/10.1002/open.70239) | Stiff-person syndrome drugs | 10 | MW, MP, density, MV, logP (Table 4), from PubChem and ChemSpider, with sources not itemised per property. |  |
| 6 | 2026-05-29 | [10.3389/fchem.2026.1848388](https://doi.org/10.3389/fchem.2026.1848388) | skin cancer drugs | 10 | MW, MR, P, complexity, PSA, MV, ST, logP (Table 1), from PubChem and ChemSpider (mixed computed and ACD predicted). |  |
| 7 | 2026-05-22 | [10.1021/acsomega.6c03590](https://doi.org/10.1021/acsomega.6c03590) | arthritis drugs (NSAIDs, coxibs, DMARDs, immunomodulators; also glucosamine, chondroitin sulfate, minocycline, doxycycline) | 50 | BP, MP (both apparently in K), Pc, MR, tPSA, CLogP (Table 3), from ChemDraw Professional 19.0 prediction cross-checked against PubChem, taking the median when sources differ. |  |
| 8 | 2026-05-19 | [10.3389/fchem.2026.1784285](https://doi.org/10.3389/fchem.2026.1784285) | major depressive disorder drugs | 21 | MW, HAC, CO (PubChem computed); BP, EV, FP, MR, PO, MV (ChemSpider ACD predicted) (Table 6) |  |
| 9 | 2026-05-13 | [10.1038/s41598-026-49153-3](https://doi.org/10.1038/s41598-026-49153-3) | antiepileptic drugs | 20 | D, IR, MR, HBA, HBD, ACD/LogD, ACD/BCF, POL, ST, MV (Table 5), source ChemSpider (ACD predicted). |  |
| 10 | 2026-04-30 | [10.1038/s41598-026-49132-8](https://doi.org/10.1038/s41598-026-49132-8) | tuberculosis drugs (eligible subset); SARS-CoV-2 subset n=8 fails (iv) | 13 | BP, E, FP, MR, P, MV for TB (Table 8); SARS-CoV-2 also has PSA and ST (Table 1). ChemSpider ACD predicted. |  |
| 11 | 2026-04-29 | [10.1038/s41598-026-50541-y](https://doi.org/10.1038/s41598-026-50541-y) | anticancer drugs (alcohol set out of scope) | 22 | MW, HA, AHA, RB, TPSA, MR, HBA, HBD, Fsp3 (SwissADME/PubChem computed, checked with RDKit) and UFF_min (RDKit computed; missing for 2 Pt compounds) (Table 4 and supplement CSV). |  |
| 12 | 2026-04-08 | [10.1038/s41598-026-45777-7](https://doi.org/10.1038/s41598-026-45777-7) | kidney infection drugs (antibiotics) | 12 | BP, E, FP, MR, PSA, P, MV (Table 1), from ChemSpider (ACD predicted). |  |
| 13 | 2026-03-24 | [10.1038/s41598-026-45284-9](https://doi.org/10.1038/s41598-026-45284-9) | influenza/antiviral drugs (also ritonavir, indinavir, atovaquone, azithromycin, nirmatrelvir) | 20 | BP, FP, MR, P, MV (Table 3). The source is stated: 'legacy.chemspider.com', which is ACD predicted. |  |
| 14 | 2026-02-17 | [10.3389/fchem.2026.1763932](https://doi.org/10.3389/fchem.2026.1763932) | WHO AWaRe Access-group antibiotics | 13 | density, BP, EV, FP, MR, P, ST, MV (Table 1), from ChemSpider (ACD predicted). |  |
| 15 | 2026-02-17 | [10.1038/s41598-026-37586-9](https://doi.org/10.1038/s41598-026-37586-9) | renal failure and complication drugs | 18 | MW, TPSA, HAC, CO (PubChem computed); BP, MR, PO, MV (ChemSpider ACD predicted) (Table 6). BP is missing for daprodustat, tenapanor and etelcalcetide. |  |
| 16 | 2026-02-02 | [10.1038/s41598-026-36855-x](https://doi.org/10.1038/s41598-026-36855-x) | anticancer drugs (several are marine bryozoan natural products) | 17 | BP, MP, E, MR (Table 1), described as 'experimental values ... from ChemSpider' (in practice ACD/EPI predicted). MP is missing for 3 drugs and E for 1. |  |
| 17 | 2026-01-14 | [10.1038/s41598-025-34694-w](https://doi.org/10.1038/s41598-025-34694-w) | amphetamine-type CNS stimulants (also furosemide, selegiline, methylphenidate) | 13 | BP, EV, FP, MR, SA (PSA), P, ST, MV (Table 3). The source is not stated; the values follow the ChemSpider ACD predicted format. |  |
| 18 | 2026-01-05 | [10.1038/s41598-025-33532-3](https://doi.org/10.1038/s41598-025-33532-3) | antidepressant drugs (plus adjuncts and non-antidepressants) | 38 | BP, MP, CT, CV, MR (Table 4), 'sourced from ChemDraw and PubChem' (predicted, despite being called experimental). MR is missing for valproate and folate. |  |
| 19 | 2025-12-23 | [10.1038/s41598-025-32906-x](https://doi.org/10.1038/s41598-025-32906-x) | anticancer drugs | 30 | BP, MR, HAC, exact mass, FP, P (Table 5): PubChem computed (HAC, EM) and ChemSpider ACD predicted (BP, MR, FP, P). BP and FP are missing for 11 drugs. | yes |

## What the screen shows about the genre (descriptive, not a modelling result)

- 19 of 60 screened papers are eligible disease/drug-class TI-QSPR studies; median n = 20, range 10–50.
- In **all** eligible studies the modelled 'physicochemical properties' are ChemSpider (ACD/Labs predicted) and/or PubChem/ChemDraw computed values; several describe them as experimental. None models a measured endpoint as its main target.
- Two studies initially judged eligible were rejected by the auditor because the per-drug values could not be traced to the paper or a journal supplement.

## Co-author decisions needed

1. Confirm rows 1–5 (or state a reason to skip one, recorded as an amendment).
2. Kara et al. 2025 (Arockiaraj co-author) is eligible but falls outside the first five; include it only if the group wishes, and label it as an added study.
3. The original R_k paper (Arockiaraj et al. 2023, IJQC) has printed-value errata (INDEX_VERIFICATION.md); if it is re-analysed, use values recomputed from its own definitions.
