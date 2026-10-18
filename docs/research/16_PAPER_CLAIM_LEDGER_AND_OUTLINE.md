# 16. Paper claim ledger and outline (doc 12 P3)

**Date:** 2026-10-16.

**Sources:** 07, 08, `PUBLICATION_READINESS.md`, doc 13 §2. Every number is from a committed file.
Items marked **PENDING** depend on experiments that have not run, and must not be written as results.

## 1. Allowed wording (exact)

| # | sentence | evidence | label |
|---|---|---|---|
| A1 | "On measured melting points, adding unlabeled-graph topological indices to element counts raised cross-validated R² by 0.085 (95% CI 0.065–0.105) for molecules whose formula was seen in training, and by 0.106 (0.087–0.125) for held-out formulas." | E3a, E3b | exploratory; pre-fit dated protocol, outcome-informed endpoint |
| A2 | "On a source-disjoint external MP set, the gain was 0.072 (0.019–0.138) for known formulas and 0.240 (0.203–0.281) for unseen formulas, where composition alone reached R² 0.23." | E3 external | exploratory, single fit, source-disjoint (not independent) |
| A3 | "Permuting MP within formula removed the gain (−0.028, −0.052 to −0.004)." | negative control | control |
| A4 | "Topology added no practically meaningful information beyond RDKit 2D descriptors (90% CI within ±0.05), and slightly reduced performance when added to ECFP on held-out formulas (−0.015, −0.029 to −0.002; BH q = 0.048)." | E3 secondary | BH over 41 secondary contrasts (42 in total with the primary) |
| A5 | "For logP, the topology gain beyond composition was practically zero (E3a +0.016, q = 0.056; E3b −0.001; 90% CIs within ±0.05)." | E3 logP | exploratory; equivalence |
| A6 | "Within the cross-validated cells, the 9 degree-pair edge counts plus maximum degree matched the full index set (differences −0.006 to −0.001 for MP; all 95% CIs include 0). On unseen external formulas the full set scored higher (0.471 vs 0.423; Δ 0.048, 0.025–0.072), but the two representations used different selected learners, so this gap is not attributed to the indices." | doc 13 §2 | exploratory; learner-confounded |
| A7 | "The index family that discriminates isomers best (distance-based) predicted MP worst under per-representation model selection." | E3 + census | exploratory; fixed-learner check pending |
| A8 | "Degree-based indices are linear functions of the edge-type counts for fixed maximum degree (Rada 2026); on our 104-index registry the feature matrix has rank 38." | 03a, 08 K1–K2 | known / incremental fact. Cite, do not claim |
| A9 | "The pre-declared within-formula R² was negative for every TI representation (E3a TIa+composition −0.67). A post hoc centred variant gives 0.19 for TIa+composition and 0.53–0.59 for labelled representations." | 07 §5 | both reported |

## 2. Forbidden wording

- **Overclaims of the result:**
  - "topological indices improve prediction" (unqualified)
  - "the gain is carried by 10 numbers" (unqualified)
  - "independent external validation" or "independent replication"
  - "preregistered" (the protocol was dated in Git and outcome-informed, not registered by a third party)
- **Novelty claims:** any new index, theorem or "discovery" of linearity or rank.
- **Generalisation claims:**
  - any positive logP, logD or logS claim
  - disease-specific predictive validity
  - "deployed v1.0 web tool"
- **Statistical misreadings:**
  - "equivalent" without a stated margin and interval
  - "non-significant, therefore no effect"

## 3. Novelty and overlap ledger

| item | status | must cite |
|---|---|---|
| Composition-controlled, formula-grouped design for unlabeled topology on measured MP | **plausibly new as a controlled measurement** (08 N1) | Basak hierarchy (design is related); Altairi/Sci Rep 2026 "combined vs physicochemical" |
| Linearity in m_ij, low rank, collisions | known | Rada MATCH 96 (2026); companion `topological-index-orthogonality` manuscript; Basak 1988 PCA |
| Discrimination ≠ prediction | incremental (08 N3) | isomer-discrimination literature (02b) |
| Endpoint dependence (MP vs logP) | expected chemistry, cleanly measured (08 N7) | additive logP literature |
| Companion manuscript overlap | **owner decision before submission**: merge, or cite and narrow | — |

## 4. Abstract (draft; bracketed items are PENDING)

> Topological indices are widely used in QSPR, but how much information unlabeled molecular topology
> adds beyond elemental composition has rarely been isolated. We compared element counts with and
> without 104 verified indices on measured melting points (Bradley, 6,725 molecules in 758 formula
> groups) and lipophilicity (OPERA, 3,994 in 495), using nested, formula-aware cross-validation and a
> source-disjoint external test. For melting point, topology added ΔR² = 0.085 (95% CI 0.065–0.105)
> within known formulas and 0.106 (0.087–0.125) across held-out formulas; a within-formula
> permutation control removed the gain. On the external set the gain was 0.072 for known and 0.240
> for unseen formulas. Topology added nothing practically meaningful beyond RDKit 2D descriptors,
> and slightly hurt ECFP. For logP the gain was practically zero. Nine degree-pair edge counts plus
> maximum degree matched the full index set in cross-validation; [the fixed-learner comparison on
> unseen formulas shows …: PENDING, doc 14]. [A confirmatory test on an untouched patent-derived MP
> source …: PENDING, doc 15.] Unlabeled topology carries real, endpoint-dependent information about
> melting point, but it is subsumed by standard labelled descriptors.

## 5. Outline

1. **Introduction:** what does unlabeled topology add over composition? Prior art: linearity and
   rank are known, and "combined vs physicochemical" studies lack a composition reference.
2. **Methods:**
   - data and licences;
   - parent standardisation;
   - formula-group inclusion (≥ 5 constitutions; the selection rule reads replicate spread only);
   - representations;
   - six learner families and nested CV;
   - E3a/E3b/external estimands;
   - corrected t and cluster bootstrap;
   - BH family (41 secondary contrasts; the primary contrast is outside it);
   - the pre-fit dated protocol, with its outcome-informed endpoint choice disclosed.
3. **Results:**
   - 3.1 MP primary (A1–A3);
   - 3.2 labelled references (A4);
   - 3.3 logP (A5);
   - 3.4 attribution: edge counts vs full set (A6, then doc 14 result, **PENDING**);
   - 3.5 discrimination vs prediction (A7);
   - 3.6 within-formula metrics (A9);
   - 3.7 confirmation (doc 15, **PENDING** or "not possible").
4. **Discussion:**
   - composition as a deliberately weak but interpretable reference;
   - why MP but not logP;
   - collisions in unlabeled graphs;
   - limits: source-disjoint is not independent, and the endpoint was outcome-informed.
5. **Reproducibility:**
   - the public dataset `docs/research/data/E3` (Bradley in full; OPERA values restorable with hash checks);
   - scripts and Colab notebooks;
   - the claims-to-evidence ledger (this file).
6. **Supplement:**
   - all 42 contrasts with absolute R²;
   - learner selections per fold;
   - failure analysis;
   - the 104-index verification (03, 03a, 03b).

## 6. Claim status after 2026-10-17 (Colab-first round; doc 18)

| claim | status | what changed |
|---|---|---|
| A1–A3, A5, A9 | **exploratory, reproduced.** Every number recomputes from the public artifacts at \|diff\| ≤ 1e-9 (notebook 01, 571 checks) | none |
| A4 | exploratory, reproduced; **BH family is 41 secondary contrasts** (42 in total) | wording fixed |
| A6 (edge counts vs full set) | exploratory; external gap learner-confounded; **fixed-learner four-arm test done** (doc 14 r2, frozen 1457bc9): edge counts carry the gain; see the two rows below | the nested contrast C+E+T − C+E replaces C+T − C+E as the test of "beyond edge counts" |
| A7 (discrimination ≠ prediction) | exploratory; **still under per-representation selection**. The four-arm run does not test it (the distance set is not an arm). | unchanged; soften or drop if not tested |
| A8 (linearity, rank 38) | known result: cite, never claim | an outside review (2026-10-17) attacks exactly this in the companion work; keep it as background |
| "TI adds beyond composition, under fixed learners" | **supported (doc 14 r2, frozen 1457bc9, run 2026-10-18).** C+T − C > δ = 0.02 in all 6 of its co-primary tests (ET, SVR × 3 contexts) (Holm p < 0.001; ΔR² 0.087–0.221). Exploratory/post hoc label applies. | `out/doc14r2_full_1457bc9/` |
| "the registry adds beyond edge counts" | **not supported under the co-primary learners.** C+E+T − C+E: 0/6 tests > 0 after Holm (min p 0.30); 3/6 equivalent within ±0.02; external ET worse by > δ (−0.035; 5 seeds −0.036 ± 0.003); external SVR inconclusive (pre-labelled underpowered). **Must also report:** linear learners (ridge, elastic net, PLS, descriptive) gain +0.03 to +0.07 internally and +0.06 to +0.14 externally. | allowed wording: "with flexible learners, the indices add nothing measurable beyond bond-type edge counts and composition" |
| "external result is independent" | **not supported, and must not be claimed.** Source-disjoint only: 3,253 of the 8,341 external test constitutions also occur in the Bradley isomer set used for the internal cells, though not in the PHYSPROP-free training subset. | evidence added |
| "the corrected-t intervals have nominal coverage" | see doc 18 §3.3 (simulation, finite-population approximation) | new |
| confirmation on the patent set | **not attempted**; census awaiting approval | — |

**Wording rule unchanged:** "preregistered", "independent validation" and "deployed v1.0" are forbidden.
