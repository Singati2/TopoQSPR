# 03 — Mathematics: exact statements, proofs, counterexamples (2026-10-13)

Evidence types are kept apart: **(P)** proof, **(I)** implementation agreement between two
codes, **(S)** agreement with a primary source's definition, **(E)** empirical computation on
data. Scripts and outputs are named at each item.

## M1. The nine degree-pair counts (owner's point 1)

**Setting.** G = (V, E) a finite simple connected graph with n ≥ 3 vertices (the
hydrogen-suppressed heavy-atom graph of a standardised, single-component parent; bond
orders and atom labels ignored), maximum degree Δ ≤ 4 (true for every organic structure the
package accepts; a fifth heavy neighbour does not occur). For 1 ≤ i ≤ j ≤ 4, let
m_ij = |{uv ∈ E : {d_u, d_v} = {i, j}}|.

**Restrictions.** (a) m_11 = 0, because an edge between two degree-1 vertices is a component K2,
excluded by connectivity with n ≥ 3. That leaves 9 possible types:
(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,3),(3,4),(4,4). (b) If Δ = 3, all m_i4 = 0 (5 types).
If Δ = 2 (paths, cycles), only (1,2),(2,2) (2 types). (c) Handshake per degree class:
i·n_i = Σ_{j≠i} m_ij + 2·m_ii, so the degree counts n_i, and hence n = Σ n_i and m = Σ m_ij,
are linear functions of the m_ij. Bond orders, heteroatoms, stereo and ring sizes are
invisible to the m_ij.

**Proposition M1 (P).** Every edge-additive degree index I_φ(G) = Σ_{uv∈E} φ(d_u, d_v)
equals Σ_{i≤j} m_ij φ(i, j). It is a linear functional of m = (m_ij). For the reverse-degree
indices, c_v = Δ − d_v + 1, so the edge weight φ_Δ(i, j) depends on Δ, and within a fixed Δ the
index is again linear in m. The modified reverse degree HCM2 is the same. *Proof:* group the
sum by the unordered degree pair of each edge; within a fixed Δ, φ is a fixed function of (i, j). ∎

**Empirical confirmation (E).** On the drug cohort, 81 of 104 registered indices (degree
family 32, reverse-degree family 49, incl. HCM2) are reproduced *exactly* (relative residual
< 1e-9) by a linear fit on the 9 counts plus intercept, separately for Δ = 3 (1024 molecules)
and Δ = 4 (792 molecules). On the Bradley isomer set the degree and reverse-degree families
distinguish exactly the same 30,902 of 35,369 within-formula pairs as (m_ij, Δ). So they are
information-equivalent to 10 numbers. **Prior art (verified 2026-10-13): ALREADY KNOWN and
already written up by the same group.** Rada, MATCH Commun. Math. Comput. Chem. 96 (2026)
841–863 (full text, §3.1) states that every vertex-degree-based index is a linear functional of
the edge-type vector. The owner's group has a public repository and manuscript in preparation,
`github.com/Singati2/topological-index-orthogonality` (Shiwakoti, Natarajan, Arockiaraj, *Ten
Edge-Degree Counts Are All a Degree-Based Topological Index Can See*, v2.1; repo created
2026-05-16, last push 2026-10-08). It turns this into a rank certificate verified on 318,648
MoleculeNet molecules. TopoQSPR must cite both and must **not** present M1 as its own result.
Doing so would be self-duplication.

## M2. Numerical rank 37 vs 38 (owner's point 1)

Matrix A = [1, n, m, 104 indices] (107 columns). **Rank 38** on both the exhaustive atlas of
connected graphs with n ≤ 7 and 3 ≤ Δ ≤ 4 (independent verification agent; pivots listed in
03a) and on the 1,816 drug-cohort molecules with all indices defined (max-abs column scaling,
pivoted QR, relative tolerance 1e-10; σ_38/σ_1 = 1.2e-7, σ_39/σ_1 = 2.0e-16). My earlier
"37" came from **z-scoring**: centring turns the constant column into zeros and removes the
intercept direction. The z-scored spectrum (σ_37/σ_1 = 3.1e-7, σ_38/σ_1 = 2.6e-16) is the
centred rank 37 = 38 − 1. There is no contradiction. The 38 = intercept + n + m + 12 degree
or R_k directions + 10 neighbourhood-sum indices + 11 distance indices + Z + I_deg (03a pivot list).
The gap between σ_38 ≈ 1e-7 and σ_39 ≈ 1e-16 is wide, so the rank is not tolerance-sensitive
between 1e-8 and 1e-14.

## M3. Representation facts (owner's mission §3)

* Discrimination of degree- vs distance-based indices on isomer sets: ALREADY KNOWN in
  substance (Dehmer, Grabner & Furtula, PLoS ONE 7:e38564, 2012, full text: degree-based indices
  have little discrimination power; distance-based do much better). The Bradley census figures
  (87.4% vs 95.0%) are a new *description of this dataset*, not a new result (02b).
* **F1 (P).** If the unlabeled heavy-atom graphs of two molecules are isomorphic, every
  unlabeled graph invariant agrees on them. *Proof:* an invariant is by definition constant on
  isomorphism classes. Atom labels and bond orders are not part of the input. **Prior art:**
  definitional; known. **(E):** all 110 TI collision classes in the drug cohort are isomorphism
  classes. On the Bradley isomer set, all 104 indices distinguish every non-isomorphic
  within-formula pair (33,592 / 33,592) and none of the 1,777 isomorphic pairs.
* **F2 (P).** A deterministic fitted predictor f applied to features T gives f(T(x)) = f(T(x′))
  whenever T(x) = T(x′). That does **not** imply that constitutional isomers collide:
  95.0% of Bradley within-formula pairs have non-isomorphic graphs.
* **F3 (P).** For squared loss, inf_f E[(Y − f(T))²] = E[Var(Y | T)], attained by f = E[Y | T]
  (law of total variance). On a finite sample, the class-mean identity of 06/SCIENTIFIC_RESEARCH_PLAN
  Proposition 1 gives the empirical analogue. A finite collision floor is an in-sample quantity.
  It conflates measurement noise and representation loss, and it says nothing about unseen
  molecules without assumptions on Var(Y | T). **Prior art:** textbook. Zopf (IJCNN 2022)
  computed analogous WL bounds for classification.
* **F4 (P).** For a connected graph, the cyclomatic number is μ = m − n + 1. The package's
  `n_rings` equals m − n + 1 on every tested graph (03a), so `size_rings` adds nothing to
  `size` = (n, m). RDKit's SSSR count differs on 37 molecules but is not used.
* **F5.** `ReZ1` ≡ n for every graph without isolated vertices (Σ_{uv}(1/d_u + 1/d_v) = Σ_v 1).
  (P) and (E): all 1,821 cohort molecules within 1.4e-13. `docs/INDEX_VERIFICATION.md` says this
  index is not registered because it equals n. It *is* registered, so documentation and registry
  disagree (LOW). It adds nothing beyond size.

## M4. Counterexamples

* **C1, representation collision.** Ethanol `CCO` and dimethyl ether `COC`: both P3, all 104
  indices equal (E). NIST normal boiling points 351.5 ± 0.2 K vs 248.2–249.5 K. RDKit 2026.03.6
  Crippen MolLogP differs (it is atom-typed). *Illustrates collision, not prediction.*
* **C2, same formula, different unlabeled topology.** Of the 35,369 within-formula pairs in the
  Bradley isomer set, 33,592 are non-isomorphic, e.g. any chain vs branched isomer pair. So
  "unlabeled TI cannot distinguish isomers" is false as a general statement (owner's plan
  item 1: VALID that the overbroad claim is false).
* **C3, discrimination does not determine prediction.** Exhibited in our own data.
  (a) TI_composition distinguishes strictly more molecules than composition (every pair
  composition separates, plus within-formula pairs), yet in the frozen run it scores lower
  than composition on logP (−0.34 vs 0.60) and logS (−568 vs 0.53). (b) asinh is strictly
  monotone, so asinh(TI) has *exactly the same* discriminating power as raw TI, yet on logP it
  predicts worse on random folds (E1: +0.037 vs +0.153 over size, worst fold R² −1.87).
  Equal discrimination can give different prediction, and more discrimination can give worse prediction.
* **C4, stereoisomers are not constitutional isomers.** Ephedrine and pseudoephedrine share
  the InChIKey first block (same constitution). Their TI collision is a stereo collision, and the
  package groups them in CV. Codeine and hydrocodone are constitutional isomers whose unlabeled
  graphs are isomorphic (an allylic alcohol vs a ketone on the same skeleton).

## M5. Hosoya Z: integer → float (owner's point 3)

`compute_indices` returns Z as an exact Python int. `index_table` stores float64, which rounds
Z > 2^53. **Measured (E):** 9 cohort molecules have Z > 2^53. The relative error is ≤ 1.09e-16
(IEEE bound 1.11e-16), the largest absolute error is 1.2e9 (on Z ≈ 5.2e25), and the
|Δ asinh| ≤ 1.1e-16. No two distinct exact values merge. None of the 9 has a measured value for
any endpoint, so they enter no model. The archived features equal the current floats. In the
Bradley isomer set (smaller molecules) Z > 2^53 does not occur. **Consequence: none for any
numerical or ML conclusion.** Documentation fix (LOW): say "exact in `compute_indices`, float64
in tables".

## M6. Correctness evidence, separated (owner's point 2)

| Evidence | Coverage | Source |
|---|---|---|
| (I) implementation agreement, independent code from documented definitions | 104/104 indices, 111,920 comparisons on 591 graphs; exact for 67 rational indices, ≤ 2.3e-15 rel otherwise | 03a |
| (P) closed forms on paths, cycles, stars (n ≤ 30) | 1,615 values | 03a |
| (S) documented definition = primary-source definition | 44 matched in the original paper; 1 later standard variant (ABC: the agent reports that the 1998 original carries a factor √2; the package follows the universal later convention; **not re-verified by me**; a constant factor is irrelevant to every model); 0 differ; 59 originals inaccessible (54 agree with a restatement; 5 with none: mM2, PI_v, edge Mostar, Hosoya Z, I_deg) | 03b |
| values printed in primary papers | 21 (R_k paper) + 50/51 (hybrid paper) + 199 cases from 7 other papers | `docs/LITERATURE_VALUE_CHECK.md` |

(I) shows the code computes what the documentation says. It is **not** proof that the
documentation matches the literature. (S) plus printed values give that for 44 + (printed-value)
indices. For the 59 inaccessible originals, the claim rests on restatements and closed forms.
Distance indices, Z and the entropy lack a primary-paper printed-value check.

## M7. Proposed new mathematics

None. No new theorem is proposed. The facts above are standard (F1–F5, M1), and the only
possibly reportable item is empirical: on 6,587 measured-MP isomers, 11 distance indices already
exhaust the discriminating power of all unlabeled invariants (95.0%), and 81 indices reduce to 10
numbers. A labelled-index family would be ALREADY KNOWN (Kier–Hall 1976; Barysz 1983; Ivanciuc
2000; 02 matrix), so none is proposed.
