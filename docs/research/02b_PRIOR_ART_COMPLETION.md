# 02b — Prior-art completion for three candidate statements (S1, S2, S3)

Date of all searches and reads: **2026-10-13**. This extends `02_NOVELTY_PRIOR_ART_MATRIX.md` and does not repeat it. Sources already listed there (Bonchev 1981, Balaban 1982, Razinger 1985, Konstantinova 1996, Dehmer 2012 *PLoS ONE* 7:e31214, Mondal 2021, Cruz–Espinal–Rada 2024, Hollas 2004, Basak 1997, Altairi 2026) are cited here only where a verdict depends on them.

Labels are the same as in 02:
- **[FT]** full text read;
- **[PP]** preprint full text;
- **[ABS]** abstract only;
- **[META]** metadata only;
- **[2H]** secondhand (search-engine summary or another source's description), not verified.

No personal identifier was sent to any service. Crossref, Europe PMC, the GitHub public API and publisher PDFs were all queried anonymously.

---

## 0. Verdicts

| ID | Statement (short) | Verdict | Decisive prior art |
|---|---|---|---|
| S1 | Within-formula pairs: the 11 distance indices separate 95.0% (every non-isomorphic pair); the 81 degree and reverse-degree indices separate only 87.4%, because they collapse to (m_ij, Δ) | **ALREADY KNOWN in substance; INCREMENTAL as a measured census on real compounds** | Rada 2026 [FT]; Dehmer, Grabner & Furtula 2012 [FT]; Bonchev 1981 / Dehmer & Mowshowitz 2013 (superindices) [ABS] |
| S2 | 81/104 indices are exactly linear in 9 m_ij (fixed Δ); the 104-index registry has numerical rank 38 | **ALREADY KNOWN** (linearity, redundancy); **INCREMENTAL** (the specific rank census). Also overlaps the group's own public companion repo. | Rada 2026 [FT]; Gutman 2013 / Rada–Cruz–Gutman 2013 [ABS/2H]; Basak et al. 1988 PCA [2H]; Singati2/topological-index-orthogonality [public repo, read via web] |
| S3 | Unlabeled TIs add no detectable held-out gain beyond 11-count composition on measured MP/logP/logD/logS (nested group CV); RDKit 2D and Morgan do better. Within-formula MP: pending. | **INCREMENTAL** (cross-formula part). **UNRESOLVED** (within-formula MP: no prior benchmark found, and the result is not in) | No composition-baseline study found (confirms 02 §2.10). Context: Hughes 2008, Lopez & Yalkowsky 2023, Katritzky 1997, Tetko 2016, Orosz 2022 |

---

## 1. S1: discrimination of within-formula pairs

### 1.1 Evidence

**Rada J. "A Finite-Dimensional Linear Framework for Vertex-Degree-Based Topological Indices." *MATCH Commun. Math. Comput. Chem.* 96 (2026) 841–863, doi:10.46793/match.96-3.32925 — [FT]** (open PDF, CC BY 4.0, match.pmf.kg.ac.rs, Match96/n3)
- §3.1, p. 848–849: "strong discrimination is unattainable for vertex-degree-based (VDB) indices on several natural classes, including chemical trees. The obstruction is simple: a VDB index depends only on the edge-type vector m(G) = (m_ij(G)) … and therefore any two graphs with the same edge-type vector necessarily receive the same value under every VDB index."
- Fig. 1 gives two non-isomorphic chemical trees S, T with identical edge-type vectors. "By attaching the same pendent path to both trees, one obtains such counterexamples with arbitrarily many vertices."
- §3.2 then defines *weak discrimination* (I(G)=I(G′) ⇒ m(G)=m(G′)), following Rada, *MATCH* 82 (2019) 29–41 [META, via Rada 2026 ref. 17].
- So the mechanism behind TopoQSPR's 87.4% ceiling is published, as a theorem-level remark with an explicit counterexample: the degree block cannot exceed the discrimination of the m_ij vector.

**Dehmer M, Grabner M, Furtula B. "Structural Discrimination of Networks by Using Distance, Degree and Eigenvalue-Based Measures." *PLoS ONE* 7:e38564 (2012), doi:10.1371/journal.pone.0038564 — [FT]** (PMC3391207, Europe PMC XML)
- Sensitivity (Konstantinova) of degree-, distance- and eigenvalue-based measures, computed on exhaustively generated graphs, on "all isomers with 11 and 12 vertices" (Table 3) and on chemical trees (Table 4).
- Results: "For the degree-based indices, it is not surprising that these measures have only little discrimination power, as many graphs can be realized by identical degree sequences. This effect is even stronger if the cardinality of the underlying graph set increases."
- The distance-based GA indices showed much higher uniqueness.
- This is the closest published *degree-vs-distance* degeneracy comparison on isomer sets. The qualitative ordering (distance > degree) is the same as S1's.

**Superindices / combinations remove most degeneracy**
- Bonchev, Mekenyan & Trinajstić 1981 (doi:10.1002/jcc.540020202, [ABS], in 02 §2.8): no single index is unique, but a combined "superindex" is far better.
- Dehmer M, Mowshowitz A. "The discrimination power of structural superindices." *PLoS ONE* 8:e70551 (2013), doi:10.1371/journal.pone.0070551 — [ABS] (Europe PMC). It compares superindices "composed of other structural indices" with the individual descriptors.
- S1's "the 11 distance indices jointly separate every non-isomorphic pair in the data" is an instance of this known effect.

**Distance indices alone are not unique**
- Balaban J and the Wiener index are known to be degenerate on larger or general graph sets: Dehmer & Shi 2014, doi:10.1371/journal.pone.0083868 [ABS], and Dehmer 2012 [FT]; for J on N10 graphs, [2H].
- S1's "every non-isomorphic pair" is therefore a property of this dataset, not a general guarantee. It must be worded as an empirical census. Larger isomer classes (for example, alkane trees with ≥ 20 vertices) are where combined distance invariants are expected to start colliding.

**Smoothness, not discrimination**
- Furtula B, Gutman I, Dehmer M. "On structure-sensitivity of degree-based topological indices." *Appl. Math. Comput.* 219 (2013) 8973–8978, doi:10.1016/j.amc.2013.03.072 — [META] + [2H].
- It studies smoothness/structure-sensitivity of 12 degree-based indices on trees, not isomer degeneracy. It is relevant as the standard "testing" reference, not as S1 prior art.

**Recent (2020–2026) work**
- New-index papers keep reporting octane- or decane-level sensitivity: Mondal 2021 (02 §2.4); the *Sci Rep* 2025 second-Davan-index paper, doi:10.1038/s41598-025-20539-z [2H]; arXiv:2605.00044 on hyper-Zagreb [2H].
- Large-scale uniqueness work uses ZINC/PubChem-scale sets (EAID indices; the COMMODE database) [2H], but for purpose-built labeled indices, not the classical unlabeled distance set.
- Not found: a 2020–2026 study that partitions a *measured-property* dataset by molecular formula and reports pairwise discrimination by index family.

### 1.2 Verdict S1: **ALREADY KNOWN in substance; INCREMENTAL in its specific measured form**

- Known, with full-text evidence:
  - the degree/reverse-degree family cannot separate graphs with equal (m_ij, Δ) (Rada 2026);
  - degree-based indices discriminate less than distance-based ones on isomer sets (Dehmer 2012);
  - combining indices removes most degeneracy (Bonchev 1981; Dehmer & Mowshowitz 2013).
- Not found, and therefore incremental:
  - the quantified gap (95.0% vs 87.4%) on 758 real formula groups / 6,587 compounds from a measured-MP dataset;
  - the observation that, *in this data*, the distance set separates every pair that is non-isomorphic as an unlabeled graph.
- State it as a dataset census that motivates the experiment, not as a discovery.

---

## 2. S2: linearity and low dimension

### 2.1 Evidence

**Linearity of VDB indices in m_ij is classical**
- Rada 2026 [FT], §1, p. 843: "If we collect the values of the defining function f into a vector γ = (f(i,j)) … I_f(G) = ⟨γ, m(G)⟩ … Thus every vertex-degree-based index is a linear functional acting on the edge-type vector."
- §2: "the analytic form of the defining function f is irrelevant: on graph classes with bounded maximum degree, a VDB index is completely determined by its values on degree pairs."
- Rada J, Cruz R, Gutman I (2013), on catacondensed hexagonal systems: "A number of vertex-degree-based topological indices can be expressed as a linear combination of the parameters m_ij" ([2H], search-engine quotation of the abstract; DOI not resolved, Crossref rate-limited).
- Gutman I. "Degree-based topological indices." *Croat. Chem. Acta* 86 (2013) 351–361, doi:10.5562/cca2294 — [META] (DOI resolves to hrcak.srce.hr; PDF fetch returned HTML, not read). It is the standard survey that writes VDB indices as Σ m_ij F(i,j) [2H].
- Cruz–Espinal–Rada 2024 (trace form; 02 §2.8, §3).

**Reverse-degree indices**
- With c_v = Δ − d_v + 1, reverse-degree indices are again an edge sum of a function of (d_u, d_v) once Δ is fixed. That makes them linear in m_ij within a fixed-Δ stratum.
- This follows immediately from the same framework (Rada 2026 §2 treats any f on degree pairs). No paper stating it for reverse-degree indices specifically was found. It is a one-line corollary, not a contribution.

**Redundancy and low dimension of TI sets are old**
- Basak SC, Magnuson VR, Niemi GJ, Regal RR. "Determining structural similarity of chemicals using graph-theoretic indices." *Discrete Appl. Math.* 19 (1988) 17–44, doi:10.1016/0166-218X(88)90004-2 — [META] + [2H].
  - PCA of 90 TIs on 3,692 TSCA chemicals: 10 PCs with eigenvalue > 1; the first 4 PCs explain ≈ 78% and 10 PCs ≈ 92% of variance; the PCs are interpreted as size, symmetry, branching and cyclicity.
  - The numbers come from a Basak slide deck (sciforum.net, MOL2NET-1), not the paper.
- Basak SC, Magnuson VR, Niemi GJ, Regal RR, Veith GD. "Topological indices: their nature, mutual relatedness, and applications." *Math. Modelling* 8 (1987) 300–305, doi:10.1016/0270-0255(87)90594-X — [META].
- Basak, Gute & Balaban (2004) clustered 318 TIs over three molecular sets ([2H]; *Croat. Chem. Acta*, not resolved).
- Hollas 2004 (covariance of degree indices with edge count; 02 §2.8).
- The Vukičević "Bond additive modeling" series (*Croat. Chem. Acta*, e.g. part 4, 2011, hrcak.srce.hr/69703 [ABS]) fits variable Adriatic indices to octanes. It was not found to state an m_ij rank result. The "Vukičević MATCH 2018" item in the brief could not be identified.

**Same-group public prior: must be handled**
- GitHub `Singati2/topological-index-orthogonality` (public; created 2026-05-16, last push 2026-10-08; read via web, not peer reviewed).
- The README describes a manuscript by G. Shiwakoti, A. Natarajan, M. Arockiaraj: "Ten Edge-Degree Counts Are All a Degree-Based Topological Index Can See: An Exact Redundancy Certificate…" (v2.1, "in preparation").
- Its claims:
  - every BID index is a linear combination of m_ij, citing Rada *MATCH* 96 (2026) 841–863;
  - a rank certificate: the 18 degree-based indices of a 30-index baseline are in the m_ij span for Δ ≤ 4, verified on 12 of 13 MoleculeNet sets (the HIV set fails because of degree 5–10 organometallics);
  - 14 published degree indices lie in the span (relative residual ≤ 2×10⁻¹²);
  - an information ceiling: 639 FreeSolv molecules map to only 236 hydrogen-suppressed graphs;
  - an RDKit-2D reference outperforms the topological models under a graph-grouped benchmark.
- These are the same authors and lab as TopoQSPR (see `docs/AUDIT_TOPOQSPR.md` on Shiwakoti/Arockiaraj). The overlap with S2 (and partly S3) is substantial. Any TopoQSPR paper must cite this work or merge with it. Otherwise it risks a self-duplication or salami-slicing objection.

### 2.2 Verdict S2: **ALREADY KNOWN (structure); INCREMENTAL (the registry-specific rank)**

- Already known:
  - "degree-based indices are linear in m_ij" (Rada 2026 [FT]; Gutman 2013; Rada–Cruz–Gutman 2013);
  - "TI sets are highly redundant / low-dimensional" (Basak 1987/1988).
- Specific numbers not found elsewhere: "rank 38 for 104 indices with 1, n, m on drug-like graphs" and "81 of 104 in the 9-dimensional span". They are a useful, checkable descriptive fact for the methods section, not a contribution.
- A numerical-rank claim needs its tolerance and conditioning stated. The companion repo itself documents an earlier tolerance artefact.

---

## 3. S3: no gain beyond elemental composition

### 3.1 Evidence

**Formula/element-count baseline vs TIs**
- New searches did not find one, in line with 02 §2.10. They covered: extended web search on "atom counts / molecular formula baseline vs fingerprints / TIs" for logP, solubility and MP; Europe PMC queries combining "topological indices" with "fingerprints" and logP or melting (0 hits); and searches on MP of isomers with QSPR descriptors (0 relevant hits).
- The nearest analogue is in another domain: amino-acid-count baselines vs fingerprints for peptides (arXiv:2501.17901; *Bioinformatics* 42 btag179, 2026) [2H].

**MP is hard, and the reason is outside the graph**
- Hughes LD, Palmer DS, Nigsch F, Mitchell JBO. *J. Chem. Inf. Model.* 48 (2008) 220–232, doi:10.1021/ci700307p — [ABS].
  - MP was "the most difficult (RMSE = 52.8 °C, r² = 0.46)" vs logP r² 0.87.
  - The authors "suggest that deficiencies in descriptors used in melting point prediction contribute significantly to the prediction errors."
- Tetko IV, Lowe DM, Williams AJ. *J. Cheminform.* 8:2 (2016), doi:10.1186/s13321-016-0113-y — [ABS].
  - Consensus MP models in drug-like space reach about the estimated experimental accuracy of 32 °C.
  - This is the noise reference for Bradley-type MP data.
- Lopez DH, Yalkowsky SH. *Pharm. Res.* (2023), doi:10.1007/s11095-023-03576-z — [ABS].
  - Over 541 compounds, external rotational symmetry number is "statistically significant in melting point" but not in boiling point, because "molecules are rotationally restricted in the crystal but not in the liquid."
  - Series of aliphatic hydrocarbons, polysubstituted benzenes and disubstituted naphthalenes illustrate the within-formula symmetry effect.
- Dannenfelser RM, Yalkowsky SH. "Estimation of Entropy of Melting from Molecular Structure: A Non-Group Contribution Method." *Ind. Eng. Chem. Res.* 35 (1996) 1483–1486, doi:10.1021/ie940581z — [META].
  - Melting entropy from symmetry number σ and flexibility. Carnelley's rule (symmetric isomers melt higher) is the standard background [2H].
- Črepnjak et al. "Predicting melting points by the Graovac–Pisanski index." *Fullerenes Nanotubes Carbon Nanostruct.* 26 (2018) 239–245; arXiv:1709.01832 — [ABS/2H].
  - A symmetry-aware (automorphism-based) Wiener variant was proposed precisely for MP.
  - For octane isomers the fit was weaker than for the alkane series.
  - This is direct precedent that ordinary distance TIs miss the symmetry part of within-formula MP variation.

**Within-formula (isomer) MP modelling exists in small form**
- Katritzky AR, Maran U, Karelson M, Lobanov VS. "Prediction of Melting Points for the Substituted Benzenes: A QSPR Approach." *J. Chem. Inf. Comput. Sci.* 37 (1997) 913–919, doi:10.1021/ci970027a — [ABS via search, 2H].
  - 443 mono- and disubstituted benzenes; a 9-descriptor model, R² = 0.837; separate models for the ortho, meta and para subsets, with quantum-chemical and H-bonding descriptors.
  - This partitions by positional isomer class, not by molecular formula, and has no held-out group design.
- Octane MP in the new-index genre (02 §2.11) is n = 18.

**Descriptors vs fingerprints**
- Orosz Á, Héberger K, Rácz A. *Front. Chem.* 10:852893 (2022), doi:10.3389/fchem.2022.852893 — [ABS].
  - 1D/2D/3D descriptors beat Morgan, AtomPairs and MACCS with XGBoost on six ADME-Tox classification targets.
  - This is consistent with S3's ordering (RDKit 2D strong), but it does not test TIs or composition.

**Same-group companion repo (§2.1)**
- It already reports that "a chemistry-aware RDKit-2D reference outperforms the topological models" under graph-grouped evaluation, and that no non-BID candidate keeps residual association once RDKit 2D is added (public README, not peer reviewed).
- This is the closest overlap to S3's second half. Again, it must be cited or merged.

### 3.2 Verdict S3

- **Cross-formula part: INCREMENTAL.**
  - The question has been asked (Altairi 2026, which answers yes; Basak 1997 hierarchy).
  - "TIs < 2D descriptors" is commonplace, and the companion repo reports it.
  - "MP depends on symmetry/packing that graphs miss" is long established (Hughes 2008; Lopez & Yalkowsky 2023; Dannenfelser & Yalkowsky 1996; Graovac–Pisanski 2018).
  - Not found anywhere: an 11-count elemental composition reference model, measured MP/logP/logD/logS, and nested group CV with ΔR² uncertainty. That combination is the niche.
- **Within-formula measured MP: UNRESOLVED.**
  - No prior benchmark with many formula groups and group-held-out CV was found. The nearest are Katritzky 1997 (o/m/p subsets) and octane sets (n = 18).
  - The novelty of the design is plausible. The finding cannot be judged until the result exists.
  - Prior art predicts that unlabeled distance TIs will capture branching but miss symmetry/packing, so a small within-formula gain would be the *expected* result, not a surprise.

---

## 4. What would be genuinely new (stated narrowly)

1. The first measured-property benchmark that holds molecular formula fixed across many groups (758 groups, 6,587 compounds, Bradley MP) and asks whether unlabeled TIs predict within-formula variation under group-held-out CV.
2. An explicit 11-count elemental-composition reference model, rather than a TI or physicochemical baseline, for unlabeled TIs on measured MP/logP/logD/logS under nested, repeated, group-aware CV with ΔR² intervals. It contests Altairi 2026's positive claim.
3. A census showing that, on this real formula-grouped dataset, the unlabeled-graph identity is fully recovered by 11 classical distance indices, while the degree block saturates at the (m_ij, Δ) ceiling (87.4%). This is descriptive support for 1–2, not a stand-alone result.
4. Not new: linearity of degree indices in m_ij, distance > degree discrimination, TI redundancy and low rank, the MP–symmetry link, and "2D descriptors beat TIs". Cite Rada 2026, Dehmer 2012, Basak 1988, Lopez & Yalkowsky 2023 and Hughes 2008 for these.
5. Before any claim: reconcile with the same group's public `Singati2/topological-index-orthogonality` manuscript (m_ij rank certificate; RDKit-2D > topology). Undisclosed overlap is the largest novelty risk.

---

## 5. Query log (2026-10-13)

| # | Tool | Query / URL | Used for |
|---|---|---|---|
| 1 | WebSearch (ext.) | degeneracy degree-based vs distance-based TIs isomer discrimination | S1 leads (Dehmer 2012 PLoS, Davan 2025) |
| 2 | WebSearch | Furtula Gutman Dehmer structure-sensitivity 2013 | S1 (AMC 2013) |
| 3 | WebSearch (ext.) | degree-based indices linear combination edge partition m_ij redundancy rank | S2; found companion repo |
| 4 | WebFetch | github.com/Singati2/topological-index-orthogonality | §2.1 overlap |
| 5 | WebSearch | Rada MATCH 2026 VDB linear | S2 |
| 6 | curl + pdftotext | match.pmf.kg.ac.rs/…/match96n3_841-863.pdf | Rada 2026 [FT] |
| 7 | WebFetch + Europe PMC fullTextXML | PMC3391207 (Dehmer, Grabner & Furtula 2012) | S1 [FT] |
| 8 | WebSearch (ext.) | Konstantinova uniqueness review Wiener/Balaban/Szeged | S1 context |
| 9 | Europe PMC search | DOIs 10.1371/journal.pone.0070551, .0038564, .0083868 | abstracts |
| 10 | WebSearch ×2 | Gutman critique of index proliferation; Vukičević bond additive modeling | S2 (no critique text found) |
| 11 | WebFetch / curl | doi.org/10.5562/cca2294 → hrcak (HTML, not PDF) | Gutman 2013 [META] |
| 12 | WebSearch (ext.) ×3 | MP isomers symmetry Dannenfelser/Yalkowsky; octane MP Randić; atom-count baseline vs fingerprints | S3 |
| 13 | Europe PMC search | Lopez & Yalkowsky 2023; Tetko 2016; MP isomers; Hughes 2008; Orosz 2022; TI+fingerprint logP/MP (0 hits) | S3 abstracts |
| 14 | WebSearch (ext.) ×2 | Mostar/Szeged/PI discrimination; large-database degeneracy | S1 |
| 15 | WebSearch | Basak PCA 90 TIs, 3,692 chemicals | S2 [2H] |
| 16 | WebSearch (ext.) ×2 | ML MP constitutional isomers; Katritzky MP substituted benzenes | S3 |
| 17 | Crossref (anonymous) | Dannenfelser–Yalkowsky; Furtula 2013; Katritzky 1997; Basak 1987/1988; Hughes 2008; ONS MP collection (doi:10.1038/npre.2011.6229) | DOIs (some queries hit rate limits) |
| 18 | GitHub public API | repos Singati2/topological-index-orthogonality, Singati2/LOYOLA-PAPER | creation and push dates |

**Not read in full** (and so not claimed beyond their labels): Gutman 2013 *CCA*, Rada–Cruz–Gutman 2013, Rada 2019 *MATCH* 82, Basak 1987/1988, Basak–Gute–Balaban 2004, Dehmer & Mowshowitz 2013, Katritzky 1997, Hughes 2008, Lopez & Yalkowsky 2023, Dannenfelser & Yalkowsky 1996, Črepnjak 2018, and the companion repo's manuscript PDF (README only).
