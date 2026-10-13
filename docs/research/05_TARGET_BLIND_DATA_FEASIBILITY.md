# 05 — Data feasibility for formula-controlled isomer experiments (Gate C)

**Terminology (corrected, owner's point 4).** Selection here is independent of target *levels*
and of any feature–target relation. It does use the package's pre-specified
replicate-agreement filter (parents whose replicate measurements disagree by > 10 °C for MP,
> 0.5 log units for logP/logS are dropped). That filter reads the *spread* of measurements. So
this census is "target-level-blind", not strictly target-blind. Justification and sensitivity
analysis: 06b R1.

## 1. Approved-drug cohorts (Experiment B data): **NO-GO**

| Endpoint | n | Formulas with ≥ 2 constitutions | ≥ 5 | Molecules in ≥ 2 groups | Largest |
|---|---|---|---|---|---|
| MP | 711 | 39 | 0 | 81 | 4 |
| logP | 505 | 24 | 0 | 49 | 3 |
| logD | 369 | 7 | 0 | 14 | 2 |
| logS | 351 | 14 | 0 | 28 | 2 |

Within-formula statistics on approved drugs are infeasible.

## 2. Full measured datasets (package rules): **GO for MP, GO for logP**

| Dataset (licence) | Parents | Constitutions | ≥ 2 | ≥ 5 | ≥ 10 | Molecules in ≥ 5 groups | Largest | Within-formula pairs | Isomorphic (unlabeled) pairs |
|---|---|---|---|---|---|---|---|---|---|
| Bradley MP (CC0) | 19,500 | 19,066 | 3,163 | **758** | 206 | **6,587** | 55 (C9H10O3) | 40,746 | 2,253 (5.5%) |
| OPERA logP (PHYSPROP; code MIT, data terms unresolved, 01b) | 13,848 | 13,847 | 2,520 | 495 | 112 | 3,994 | 43 | 23,865 | 1,357 (5.7%) |
| AqSolDB logS (CC0) | 8,417 | 8,136 | 1,305 | 222 | 63 | 1,904 | 35 | 12,053 | 600 (5.0%) |

(Pairs counted over all groups ≥ 2. Within groups ≥ 5 the Bradley count is 35,369 pairs,
1,777 isomorphic.)

**Reconciliation of 769 vs 758 (owner's question).** Raw RDKit parsing of the Bradley
spreadsheet (no `donotuse`, single fragment, no standardisation) gives 769 groups / 6,749
molecules. The package's standardisation gives the same 769 / 6,748. The replicate-agreement
filter (299 parents dropped of 19,799) gives **758 / 6,587**. The 11 lost groups exist only
before the filter. No group appears only after it. Duplicates are merged by full standardised
InChIKey (median of records). 368 constitutions have stereo variants, which are kept as separate
parents and grouped in CV.

## 3. Information content (target-level-blind; `scripts/family_information_content.py`)

Within-formula pairs (35,369; 6,587 constitutions; 758 formulas) distinguished:
composition 0%; size 19.9%; degree family = reverse-degree family = (9 edge counts + Δ) 87.4%;
degree entropy 62.9%; Hosoya Z 94.1%; neighbourhood-degree-sum 94.8%; **distance family 95.0%
= all 104 TI = bound for any unlabeled invariant (non-isomorphic pairs)**; ECFP counts 99.95%;
RDKit 2D 100%.

**Consequences.** (i) Unlabeled TI can in principle express 95% of within-formula structural
variation, so a predictive failure would not be a discrimination failure. (ii) Within TI,
discrimination is redundant beyond the 11 distance indices. (iii) **Discrimination is not
prediction (03 M4 C3).** Whether any of this predicts MP is tested only by the held-out paired
contrasts of E3.

## 4. Source independence for external validation (owner's point 5)

The Bradley `source` field names PHYSPROP for 9,687 of 28,645 records. 8,811 of 19,500 kept
parents have ≥ 1 PHYSPROP record. OPERA MP is a PHYSPROP curation. **A Bradley-trained model
tested on OPERA MP is therefore not source-independent.** Design (06b R3): train only on Bradley
parents without any PHYSPROP record (10,689 kept parents; 2,252 in 305 isomer groups), test on
OPERA parents whose constitution is absent from training (1,014 known-formula molecules in 211
formulas; 7,327 unseen-formula molecules). The residual caveat stands: different compilations can
cite the same primary measurement. Source-disjoint is weaker than proven
measurement-independence.

## Verdict

GO (exploratory) for MP isomer study, internal E3a/E3b plus source-disjoint external test.
logP is secondary and budget-dependent. NO-GO on drug cohorts.
