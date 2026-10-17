# 08 — Novelty and submission decision (Gate F), 2026-10-14

Inputs: 02/02b (prior art, read 2026-10-13), 03 (mathematics), 05 (data), 07 (E3 results, all
cells complete 2026-10-16). This document answers the owner's question:

> What genuinely new scientific understanding does TopoQSPR establish that existing mathematical
> graph-theory and QSPR papers have not established?

## 1. Short answer

**No new mathematics.** The empirical contribution is narrow, but it is real and measured. On a
measured endpoint (melting point), with element composition held as the reference and the
molecular formula controlled, unlabeled topology has a reproducible, practically meaningful value
beyond composition, both within and across formulas. That value is carried entirely by 10 numbers:
the degree-pair edge counts and the maximum degree Δ. It is redundant once a standard labelled
representation (RDKit 2D, ECFP) is present. The family that discriminates isomers best (distance
indices) predicts worst. As far as the 2026-10-13 prior-art search reaches, nobody has shown these
three statements together under this design on a measured endpoint.

## 2. Ledger: new vs known

| # | Statement | Evidence | Status |
|---|---|---|---|
| K1 | Every degree-based index is a linear functional of the 9 edge-type counts m_ij (fixed Δ). 81/104 registry indices reduce to (m_ij, Δ) | 03 M1 (P, E) | **ALREADY KNOWN**: Rada, MATCH 96 (2026); the group's own `topological-index-orthogonality` manuscript. Cite; do not claim |
| K2 | The registry matrix [1, n, m, 104 TI] has rank 38 | 03 M2 | **INCREMENTAL** registry fact. Redundancy and low rank of TIs are known (Basak 1988 PCA) |
| K3 | Distance indices discriminate isomers better than degree indices (95.0% vs 87.4% of Bradley within-formula pairs). All 104 reach the isomorphism bound | 03 M3, census | **ALREADY KNOWN** in substance (Dehmer, Grabner & Furtula 2012). The Bradley census is a description of the dataset |
| K4 | Isomorphic unlabeled graphs collide; the Bayes floor is E[Var(Y\|T)] | 03 F1, F3 | **ALREADY KNOWN** (definitional and textbook) |
| K5 | Labelled 2D descriptors and fingerprints beat TIs | 07 §3 | **ALREADY KNOWN** and commonplace; also in the companion repo |
| K6 | MP depends on symmetry and packing that graphs miss; amino acids and zwitterions are hard | 07 §6 | **ALREADY KNOWN** (Dannenfelser & Yalkowsky 1996; Hughes 2008; Lopez & Yalkowsky 2023) |
| **N1** | On measured MP, TI + composition beats an 11-count composition reference by ΔR² +0.085 [0.065, 0.105] within known formulas (E3a). It does so by +0.106 [0.087, 0.125] across held-out formulas (E3b) and by +0.07 / +0.24 on a source-disjoint external set. The result survives the negative control and the S-noflt sensitivity run | 07 §2–4 | **PROBABLY NEW in its controlled form.** The question is known (Altairi 2026 answered "yes" with computed targets, a single split and test-set model selection; Basak's hierarchy). The design (composition reference, formula control, nested repeated CV, source-disjoint external test, negative control) was not found. **Exploratory; one endpoint** |
| **N2** | The whole topology gain is carried by (m_ij, Δ). edgetypes_composition ≥ TIa_composition internally and on external known formulas; on external unseen formulas it is 0.42 vs 0.47 | 07 §3 | **INCREMENTAL / probably new as a predictive statement.** K1 says degree indices *are* the m_ij. N2 says that for MP prediction nothing outside the m_ij span adds much. Overlaps the companion manuscript's theme, so the two papers must be coordinated |
| **N3** | Discrimination ≠ predictive value on a measured endpoint. Distance indices discriminate most but predict worst; on unseen external formulas they add nothing (−0.028 [−0.058, +0.002]) | 07 §3, 03 M4 C3 | **PROBABLY NEW as a measured demonstration.** The general caution is folklore. No paper found tests it with formula-controlled MP and group-held-out CV |
| **N4** | Topology's incremental value beyond RDKit 2D is practically zero (90% CIs inside ±0.05). Beyond ECFP it is zero or *negative* (E3b −0.015 [−0.029, −0.002]; external unseen −0.023) | 07 §3 | **INCREMENTAL.** K5 is known; the equivalence-bounded statement with a significant harm is a clean addition |
| N5 | Within formula, unlabeled topology explains about 18–21% of MP variance; labelled representations explain 52–59% | 07 §5 | **Post hoc.** Descriptive only |
| **N7** | Topology's added value is endpoint-dependent. For logP, composition alone reaches R² 0.77–0.81 and TI adds +0.016 (E3a, q 0.056) or −0.001 (E3b), practically equivalent to 0 (90% CIs inside ±0.05). For MP, which depends on shape and packing, it adds +0.085–0.106 | 07 §3b | **INCREMENTAL.** It is expected chemistry (logP is close to additive in atom and fragment contributions; MP is not). It is still a clean controlled contrast and makes the MP result's scope explicit |
| N6 | Drug-cohort endpoints (E1/E2): no detectable gain beyond composition for logP/logD/logS; MP was underpowered (n = 711) and is consistent with E3 | audit RESULT_RECONCILIATION §7 | Negative and inconclusive results; keep in the paper |

## 3. Decision

**CONDITIONAL GO** for a single evaluation-led research article. Working title:
*"What does unlabeled molecular topology add? A formula-controlled test on measured melting points."*
Target: Journal of Cheminformatics (Research article). Fallback: Digital Discovery.

**Not GO** as a mathematics paper, a new-index paper or a software paper on its own merits.

### Conditions (all required before submission)

1. **Done (2026-10-16):** logP E3 ran, and the BH family is final (41 secondary contrasts; 42 with the primary, which is not in the family). The paper must
   report the endpoint dependence. Topology adds a practically meaningful amount for MP. For logP it
   adds nothing practically meaningful beyond composition (N7).
2. **Companion-manuscript coordination.** `topological-index-orthogonality` (same authors) covers K1
   and the m_ij ceiling. Either merge, or cite it and confine this paper to N1–N4. Undisclosed
   overlap is the largest rejection risk.
3. **Wording.** Use "exploratory, dated amendment committed before fitting", never "preregistered".
   Use "source-disjoint", never "independent", for the external test. Name both the primary result
   and the dominance by RDKit 2D/ECFP in the abstract. Present no discrimination result as a
   predictive claim.
4. **Licensing and repository hygiene** (10, 11): the owner decides on the copyrighted PDF in
   history; DATA_LICENSES/ChEMBL attribution; OPERA terms. Release no per-molecule OPERA data.
5. **Confirmation, preferably.** One confirmatory replication of N1/N2 on a second measured
   endpoint, or an untouched MP source, with the protocol frozen before any fit. Without it, the
   paper is an exploratory evaluation, and the title and abstract must say so.

## 4. What the paper must not claim

"Topological indices improve property prediction" (they do not beyond standard descriptors) ·
"new indices" · "the degree-index linearity / rank result" as new · "independent external
validation" · "preregistered" · "95% isomer discrimination implies predictive power" · any
statement about logP, logD or logS stronger than E2's "no detectable gain; equivalence not shown"
· any disease-specific predictive validity.
