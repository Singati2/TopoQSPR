# 02 — Novelty and prior-art matrix (TopoQSPR, candidate framing "What does molecular topology add?")

Date of all searches and reads: **2026-10-13**. Extends `docs/audit-2026-10-12/LITERATURE_2026.md`
(tool matrix, journal fit, benchmark practice, genre papers); nothing there is repeated except
where a verdict depends on it.

Labels:
- **[FT]** full text (or the relevant full sections) read from a primary copy.
- **[PP]** preprint full text read; the peer-reviewed version was not read.
- **[ABS]** abstract only (PubMed, Crossref, Europe PMC or arXiv API).
- **[META]** bibliographic metadata only (Crossref), content not read.
- **[2H]** secondhand: search-engine summary or another paper's description. Not verified.

DOIs were resolved through Crossref on 2026-10-13 unless marked otherwise.

---

## 1. Verdict summary

| ID | Candidate contribution | Verdict | One-line reason |
|---|---|---|---|
| A | Benchmark finding: unlabeled TIs add no detectable gain beyond element composition on **measured** endpoints (nested CV) | **INCREMENTAL, with a real niche** | The question "do TIs add beyond a simple baseline?" has been asked before: Basak's hierarchical QSAR, and Altairi et al. *Sci Rep* 2026, which answers *yes*. No prior study found that uses a composition/formula baseline, measured endpoints and nested CV. |
| B | Representation-limit facts: (i) isomorphic graphs give equal invariants; (ii) Bayes floor E[Var(Y\|T)]; (iii) empirical collision floor | **ALREADY KNOWN** (i, ii); **INCREMENTAL** (iii) | (i) is the definition of a graph invariant. (ii) is a textbook result. (iii) has been done for 1-WL/GNN *classification* on TU molecule sets (Zopf 2022, majority-vote upper bound), not for TI *regression* on measured properties. |
| C | Incremental value of TIs over a baseline: none found | **INCREMENTAL** (finding); the design is **ALREADY KNOWN** | Incremental-value designs exist (Basak hierarchy; Sci Rep 2026 "Combined vs Physicochemical"). The *null* result under nested CV against composition was not found, and it contradicts Sci Rep 2026's positive claim. |
| D | Disease (MeSH) → ChEMBL approved-drug cohort software | **INCREMENTAL** (software convenience) | ChEMBL already exposes `drug_indication` with `mesh_heading` and `max_phase_for_ind` filters (API checked). No TI/QSPR tool packaging this was found (prior audit). |
| E | A new atom/bond-labeled TI family | **ALREADY KNOWN** as a concept; any specific construction would be at most **INCREMENTAL** | Covered by Kier–Hall valence δ (1976), Barysz weighted distance matrices (1983), the Ivanciuc Z/X/Y weighting schemes (2000), alvaDesc valence/Kupchik Zagreb variants, Mordred Barysz matrices, and Nadjafi-Arani 2025 (VEW TIs weighted by six atomic properties). |
| F | Formula-controlled isomer benchmark on Bradley MP | **INCREMENTAL / probably NOVEL in its specific form** | Fixed-formula benchmarks exist: octane isomers since Randić 1975 (measured), and QM9's 6,095 C7H10O2 isomers and ISO17 (computed). No benchmark found that uses measured MP on many formula groups with group-held-out CV. |

**Remaining niche, stated precisely.** A leakage-resistant *re-test* of unlabeled TI value:
- on measured endpoints;
- against **size and element-composition baselines**;
- with repeated nested, group-aware CV and uncertainty on ΔR²;
- with an explicit decomposition of what an unlabeled graph cannot see (isomorphism collisions, within-formula variance).

It directly contests the only recent paper with the same question, Altairi et al. 2026 (§2.1). The tool (D) and any new index (E) are not publishable novelties on their own.

---

## 2. Paper-by-paper evidence

### 2.1 Altairi A, Alhaj Z, Alsharafi M, Zeren Y. "From graph theory to chemoinformatics: modified bond-based indices and a hypothesis-driven multi-task QSAR/QSPR benchmark." *Sci Rep* 16 (2026), doi:10.1038/s41598-026-40969-7 — **[FT]** (Europe PMC PMC13021970, full-text XML; received 2025-11-26, accepted 2026-02-17)

This is the strongest overlap. Formulas in the XML are images, so the index definitions below come from the prose.

**Indices.** There are 16 "modified bond-based indices". Each is an edge sum with a vertex kernel f(d_u, d_v) in the numerator and a "bond factor" in the denominator. The factor is zero only for K₂ ("both endpoints have degree 1"), and the paper calls it "bond degree", so it is d_u + d_v − 2 (inferred from that prose). Section "Computational complexity" says each summand "depends only on endpoint degrees and the bond factor", so the indices are **unlabeled degree-based invariants**. The paper gives a unified edge-partition representation, a sum over degree classes (i, j) (Abstract; "Chemoinformatics perspective" section).

**Data.**
- 3,219 ChEMBL "antibacterial" molecules, de-duplicated by (ChEMBL ID, InChIKey) ("Dataset and target properties").
- **Endpoints:** 10 continuous targets. Nine are "taken directly from ChEMBL-provided physicochemical fields": MW, AlogP, PSA, HBA, HBD, rotatable bonds, QED, heavy atoms, NP-likeness (Table 4 and its abbreviations; "AlogP = calculated logP"). The tenth is pMIC from heterogeneous MIC records. **Nine of ten targets are computed descriptors, not measurements.** Only pMIC is measured, and it is pooled across organisms and assays.

**Baselines.**
- "Physicochemical descriptors" scenario: the same 10 ChEMBL fields, with the column equal to the target removed ("leakage prevention", eq. 7).
- No element-composition baseline. No fingerprint baseline: fingerprints appear only as future work in "Conclusions and future directions". No size-only baseline, although HA and MW sit inside the physicochemical set.

**Splits and protocol.**
- "a single train–test split was defined at the molecule level and reused for every target and scenario" ("Multi-task formulation and data splitting"). It is a random molecule-level split; the ratio sits in image-only equations, so I could not read it. There is no scaffold or group split.
- "three-fold cross-validation" is used for "internal stability". The best model per target and scenario is then **"selected as the regressor minimising test-set RMSE"** ("Training protocol…").
- So there is **no nested CV**, and model selection uses the test set.
- No confidence intervals, repeats or significance tests are reported for scenario differences.

**Result and claim.**
- Combined vs Physicochemical: macro test R² **0.861 vs 0.852**, global zRMSE 0.373 vs 0.385. TI alone: 0.673 (Table 5).
- The paper concludes the TIs "provide complementary signal that yields consistent, multi-target improvements" (Discussion).
- Per-target, TI-only reaches R² 0.99 for heavy-atom count and 0.94 for MW (Table 7). Those targets are size, and an unlabeled-graph index recovers them trivially.
- Data and code: "processed dataset is included in the Supplementary Information". Scripts are "available from the corresponding author upon reasonable request" (Data availability). The SI was **not** read.

**Overlap verdict.** Same question as TopoQSPR's (A)/(C): do unlabeled TIs add beyond a standard descriptor baseline? The opposite answer comes from a weaker design:
- computed targets;
- no composition baseline;
- one split;
- test-set model selection;
- a +0.009 macro-R² gain with no uncertainty.

TopoQSPR should cite it as the direct counterpoint. Report that under nested, repeated, group-aware CV on measured endpoints with a composition baseline, the gain is not detectable. Do not claim the question is new.

### 2.2 Ivanciuc O. "QSAR Comparative Study of Wiener Descriptors for Weighted Molecular Graphs." *J Chem Inf Comput Sci* 40(6):1412–1422 (2000), doi:10.1021/ci000068y — **[ABS]** (PubMed 11128100; ACS full text 403)

From the abstract:
- Wiener-type descriptors, including the "Wiener graph operator", computed for "vertex- and edge-weighted molecular graphs, corresponding to organic compounds with heteroatoms and multiple bonds";
- applied to acute toxicity of 47 nitrobenzenes to *Tetrahymena pyriformis*, with hydrophobicity as a co-descriptor.

The methods and supplement were not accessible.

### 2.3 Ivanciuc O, Ivanciuc T, Cabrol-Bass D, Balaban AT. "Comparison of Weighting Schemes for Molecular Graph Descriptors: Application in QSRR Models for Alkylphenols in Gas–Liquid Chromatography." *J Chem Inf Comput Sci* 40(3):732–743 (2000), doi:10.1021/ci990129w — **[ABS]** (PubMed 10850777)

From the abstract:
- three atom/bond parameter sets: atomic number **Z**, relative electronegativity **X**, relative covalent radius **Y**;
- descriptors from Wiener, hyper-Wiener, min/max eigenvalue, Ivanciuc–Balaban and information-on-distance operators, applied to the distance and reciprocal-distance matrices;
- endpoint: GC retention of alkylphenols.

**Relevance to E:** atom- and bond-labeled ("weighted") molecular graphs with interchangeable weighting schemes were systematised by 2000.

### 2.4 Mondal S, Dey A, De N, Pal A. "QSPR analysis of some novel neighbourhood degree-based topological descriptors." *Complex Intell Syst* 7(2):977–996 (2021), doi:10.1007/s40747-020-00262-0 — **[ABS]** for the journal version (Crossref abstract; CC BY 4.0; Springer PDF blocked by a login redirect) + **[PP]** for arXiv:1906.06660 (6-page preprint, read in full)

- Indices ND1–ND6 are built on neighbourhood degree sums (preprint §2–3).
- **Data:** 18 octane isomers (entropy, acentric factor, HVAP, DHVAP, Cp; "taken from www.moleculardescriptors.eu", preprint Table 1) and 67 alkanes (bp, mv, mr, hv, st…; "compiled from [16]", Table 9). These are measured literature values, all **hydrocarbons**. Labeled vs unlabeled is moot for them.
- **Degeneracy (preprint §5 and Table 19):** uses Konstantinova's sensitivity S_I = (N − N_I)/N, where N_I is the number of isomers the index cannot distinguish. It is evaluated on octane and decane isomers.
- Correlation with known indices (Table 17): ND1 vs F, r = 0.995.

**Relevance:**
- isomer degeneracy of exactly TopoQSPR's NDS family is already measured on fixed-formula alkane sets (B, F);
- the new indices are near-collinear with classical ones (r ≥ 0.95 for most).

The journal abstract also mentions "isomer discrimination". Whether the journal version changes tables or sections was not verified.

### 2.5 Landrum GA, Beckers M, Lanini J, Schneider N, Stiefl N, Riniker S. "SIMPD: an algorithm for generating simulated time splits…" *J Cheminform* 15:119 (2023), doi:10.1186/s13321-023-00787-9 — **[FT]** (PMC10712068)

- A multi-objective GA that splits public data so that train/test differ the way early and late compounds differ in >130 NIBR lead-optimisation projects. There are eight objectives, e.g. Δ median SA_Score = 0.28 (section "Objectives").
- Built and validated on **bioactivity classification**: pAC50 threshold 6.3, RF classifier on Morgan r=2 fingerprints, balanced accuracy/F1/κ (Methods "Data sets", "Model building").
- 99 ChEMBL-derived split sets are released.

**Relevance:** the simulated-time concept targets medicinal-chemistry project series. A cross-sectional cohort of approved drugs from many diseases has no project or time axis. SIMPD is cited as the leakage-realism reference, but it is not a natural primary split for TopoQSPR's question.

### 2.6 Wu Z et al. "MoleculeNet: a benchmark for molecular machine learning." *Chem Sci* 9:513–530 (2018), doi:10.1039/C7SC02664A — **[FT]** (PMC5868307)

- Table 1 recommends **random** splits for all physical-chemistry regression sets: ESOL 1,128, FreeSolv 643, **Lipophilicity 4,200**. Scaffold splits are recommended only for HIV, BACE and BBBP (classification).
- MoleculeNet Lipophilicity comes from the same AstraZeneca ChEMBL deposition that TopoQSPR uses for logD7.4 (the acknowledgements credit its contribution; the count of 4,200 matches the repo's `logd74_az` raw count, `docs/EXPERIMENTAL_DATA.md`).

**Relevance:** the default benchmark for logD is random-split. TopoQSPR's group-aware nested CV on the same source is a stricter replication axis.

### 2.7 Guo Q, Hernandez-Hernandez S, Ballester PJ — scaffold splits overestimate

- Preprint arXiv:2406.00873 (2024) and **ICANN 2024**, LNCS pp. 58–72, doi:10.1007/978-3-031-72359-9_5 — [META] / [2H]. Already in the prior audit.
- **Peer-reviewed successor:** "UMAP-based clustering split for rigorous evaluation of AI models for virtual screening on cancer cell lines", *J Cheminform* 17:94 (2025), doi:10.1186/s13321-025-01039-8 — **[FT]** (PMC12153141; received 2024-12-28, accepted 2025-05-27).
  - NCI-60, 60 cell lines with ~33k–54k molecules each; 7-fold merged CV × 5 seeds; LR, RF, GEM (GNN) and Transformer-CNN.
  - **Continuous vs classification:** models regress continuous pGI50. The **primary metric is hit rate** among the top-100 predictions after dichotomising at pGI50 > 6, with MCC, ROC-AUC and RMSE as secondary ("Test set performance metrics").
  - The difficulty ordering (UMAP > Butina > scaffold ≈ random) is therefore established mainly in a thresholded virtual-screening sense. That is not the same as R² on physicochemical regression.
  - TopoQSPR should cite this as evidence that scaffold grouping is still optimistic. It should not claim the ordering transfers unchanged to logP/MP regression.

### 2.8 Classic degeneracy and discrimination literature

| Reference | Status | What it establishes |
|---|---|---|
| Randić M. "Characterization of molecular branching." *JACS* 97:6609–6615 (1975), doi:10.1021/ja00856a001 | [META] | Connectivity index; alkane/octane branching as the canonical fixed-formula test. |
| Bonchev D, Mekenyan O, Trinajstić N. "Isomer discrimination by topological information approach." *J Comput Chem* 2:127–148 (1981), doi:10.1002/jcc.540020202 | [ABS] | "No index is found to discriminate isomers uniquely". Proposes a "superindex" tested on 427 acyclic, mono- and bicyclic graphs with 4–8 vertices. |
| Balaban AT. "Highly discriminating distance-based topological index." *Chem Phys Lett* 89:399–404 (1982), doi:10.1016/0009-2614(82)80009-2 | [META] | J index designed for low degeneracy. Still degenerate on larger graph sets (see Dehmer 2012). |
| Razinger M, Chrétien JR, Dubois JE. "Structural selectivity of topological indexes in alkane series." *JCICS* 25:23–27 (1985), doi:10.1021/ci00045a006 | [META] | Selectivity/degeneracy of TIs over alkane series. |
| Konstantinova EV. "The Discrimination Ability of Some Topological and Information Distance Indices for Graphs of Unbranched Hexagonal Systems." *JCICS* 36:54–57 (1996), doi:10.1021/ci9502461 | [META]; the sensitivity definition is [2H] via Mondal preprint §5 | Defines the sensitivity S = (N − N_I)/N. |
| Dehmer M, Grabner M, Varmuza K. "Information Indices with High Discriminative Power for Graphs." *PLoS ONE* 7:e31214 (2012), doi:10.1371/journal.pone.0031214 | [ABS] (Europe PMC) | Uniqueness tested on ~12 million exhaustively generated non-isomorphic unweighted graphs. A degree-degree information index "outperforms the Balaban J index tremendously". |
| Kier LB, Hall LH. "Molecular Connectivity VII: Specific Treatment of Heteroatoms." *J Pharm Sci* 65:1806–1809 (1976), doi:10.1002/jps.2600651228 | [META] | Valence-delta (δᵛ) connectivity, i.e. **atom-labeled** connectivity indices. |
| Barysz M, Jashari G, Lall RS, Srivastava VK, Trinajstić N. "On the distance matrix of molecules containing heteroatoms." In King RB (ed.) *Chemical Applications of Topology and Graph Theory*, Elsevier 1983, pp. 222–230 | [2H] (VCCLab descriptor help page; no DOI found) | Z-weighted distance (Barysz) matrix: atom- and bond-labeled Wiener-type indices. Mordred implements Barysz matrices weighted by Z, mass, electronegativity, etc. (prior audit S5). |
| Basak SC, Gute BD, Grunwald GD. "Use of Topostructural, Topochemical, and Geometric Parameters in the Prediction of Vapor Pressure: A Hierarchical QSAR Approach." *JCICS* 37:651–655 (1997), doi:10.1021/ci960176d | [META] + [2H] | **Closest conceptual precedent for A/C.** Unlabeled ("topostructural", TS) indices are entered first, then labeled ("topochemical", TC), then 3-D. A group slide deck (not the paper) gives TS R² ≈ 0.48/0.58 vs TC 0.79/0.86 for 476 chemicals [2H, unverified]. The baseline is a hierarchy of index types, not element composition. |
| Hollas B. "The Covariance of Topological Indices that Depend on the Degree of a Vertex." arXiv:math/0404512 (2004) | [ABS] | Correlation between edge-additive degree indices and the edge set (size). Already in the prior audit. |
| Cruz R, Espinal C, Rada J. "A Matrix Approach to Vertex-Degree-Based Topological Indices." *Mathematics* 12:2043 (2024), doi:10.3390/math12132043 | [ABS] | Every VDB index equals trace(φᵀ α(G)), where α(G) is the matrix of edge counts by degree pair. See §3. |

### 2.9 Representation collisions, WL limits, information and Bayes floors

| Reference | Status | Relevance |
|---|---|---|
| Xu K, Hu W, Leskovec J, Jegelka S. "How Powerful are Graph Neural Networks?" ICLR 2019, arXiv:1810.00826 (no DOI) | [ABS] | MPNNs are at most as discriminative as 1-WL. GIN reaches the 1-WL bound. |
| Morris C et al. "Weisfeiler and Leman Go Neural: Higher-Order GNNs." *AAAI* 33:4602–4609 (2019), doi:10.1609/aaai.v33i01.33014602 | [META] | Same bound; k-GNNs. |
| Morris C, Fey M, Kriege NM. "The Power of the Weisfeiler-Leman Algorithm for Machine Learning with Graphs." IJCAI 2021, pp. 4543–4550, doi:10.24963/ijcai.2021/618 | [META] | Survey. |
| **Zopf M. "1-WL Expressiveness Is (Almost) All You Need." IJCNN 2022, pp. 1–8, doi:10.1109/IJCNN55064.2022.9892655** | **[PP]** (arXiv:2202.10156v1, read §VI and Tables I–VII) | **Closest prior art for B(iii).** (1) "Fraction of identifiable graphs" per dataset. (2) An upper bound on classification accuracy computed by majority vote within each cluster of identical WL representations (§VI-B). (3) The same analysis *without node labels* (§VI-C, Table IV) and for k = 0, a pure node-label histogram (Table V). The data are TU benchmarks including the molecule sets **MUTAG, NCI1, NCI109**. Example: NCI1 without labels at k=0 has 0.28% identifiable graphs and a 63.70% accuracy bound; with labels, 56.34% and 91.34% (Table V). It is classification only; no regression or E[Var(Y\|T)] analogue was found. |
| Adamczyk J, Czech W. "Molecular Topological Profile (MOLTOP) — Simple and Strong Baseline for Molecular Graph Classification." ECAI 2024, doi:10.3233/FAIA240663 | [ABS] (arXiv:2407.12136) + [2H] | Topological descriptors plus one-hot atom/bond types with RF are a strong GNN baseline. Uses the **decalin vs bicyclopentyl** 1-WL-indistinguishable molecule pair, attributed to earlier work (primary source not located; [2H]). |
| Virany W, Tripp A. "Hash Collisions in Molecular Fingerprints: Effects on Property Prediction and Bayesian Optimization." arXiv:2511.17078 (2025) | [ABS] | Exact (collision-free) fingerprints give "small yet consistent" gains on five DOCKSTRING regression tasks. Representation non-injectivity from hashing, measured on regression. |
| Gütlein M, Kramer S. "Filtered circular fingerprints…" *J Cheminform* 8:60 (2016), doi:10.1186/s13321-016-0173-z | [META] | Folding/bit-collision effects in ECFP. |
| Teufel J et al. "Hyper-Dimensional Fingerprints as Molecular Representations." arXiv:2604.27810 (2026) | [ABS] + [2H] | Hash-based fingerprints lose size- and composition-dependent information. |
| Brown SP, Muchmore SW, Hajduk PJ. "Healthy skepticism: assessing realistic model performance." *Drug Discov Today* 14:420–427 (2009), doi:10.1016/j.drudis.2009.01.012 | [META] + [2H] | Assay error caps attainable R². |
| Kramer C, Kalliokoski T, Gedeck P, Vulpetti A. "The Experimental Uncertainty of Heterogeneous Public Ki Data." *J Med Chem* 55:5165–5173 (2012), doi:10.1021/jm300131x | [META] + [2H] | Experimental uncertainty sets the upper limit for in-silico models. |
| Crusius D, Cipcigan F, Biggin PC. "Are we fitting data or noise? …" *Faraday Discuss* 256:304–321 (2025), doi:10.1039/D4FD00091A | [META] | Noise ceilings of common datasets. |

**Search result on "Bayes error of molecular representations / non-injective representation floor for regression":**
- no paper found that computes E[Var(Y|T)] over representation-collision classes for molecular *regression* (WebSearch extended, three phrasings; Consensus search).
- The idea exists for WL classification (Zopf) and for assay noise (Brown, Kramer, Crusius). The regression analogue for TIs is a small, legitimate addition.
- `AUDIT_FINDINGS.md` F9 states the empirical ceilings are **not binding** on TopoQSPR's data (TI-alone ceiling R² 0.94–0.997). So B(iii) is a diagnostic, not a headline.

### 2.10 Prior tests of TIs against an element-composition or formula baseline

Searches:
- WebSearch, extended: "QSPR baseline atom counts OR elemental composition outperforms topological indices";
- WebSearch, extended: "topological indices compared with molecular formula baseline";
- Consensus: "Do topological indices add predictive information beyond molecular formula or atom counts in QSPR?" (20 results).

Result: **no study found that uses element counts or formula as the reference model for TIs on a property set.**

Nearest neighbours:
- Basak's TS→TC hierarchy (labeled vs unlabeled, §2.8).
- The in-genre admission of size confounding (Çolakoğlu 2026, prior audit §3c).
- Sci Rep 2026's physicochemical baseline (§2.1).
- Fixed-formula isomer benchmarks (§2.11), which control composition by design.

Materials informatics routinely uses composition-only models as the baseline that structure must beat; the TopoQSPR discussion can cite this analogy:
- Goodall & Lee, "Predicting materials properties without crystal structure: deep representation learning from stoichiometry", *Nat Commun* 11:6280 (2020), doi:10.1038/s41467-020-19964-7, [META].

### 2.11 Formula-controlled isomer benchmarks (for F)

- **Octane/alkane isomers.** Measured thermodynamic properties of the 18 octanes, used since Randić 1975 and in every new-index paper (e.g. Mondal preprint Table 1, from moleculardescriptors.eu). Fixed formula, hydrocarbons only, n = 18.
- **QM9 C7H10O2 subset.** 6,095 constitutional isomers, figshare article 1057645 (CC0), Ramakrishnan et al. *Sci Data* 1:140022 (2014), doi:10.1038/sdata.2014.22. Computed (DFT; G4MP2 for this subset per [2H]).
- **ISO17.** MD snapshots of C7O2H10 isomers with "known vs unknown connectivity" test sets. Computed energies and forces ([2H]).
- **J Cheminform 2025, doi:10.1186/s13321-025-01062-9.** Fingerprints plus a TI improve isomer discrimination for boiling point on 155 molecules ([2H], not opened).
- **Not found:** a benchmark of measured MP (or another measured property) over many formula groups with group-held-out CV (unseen constitution within a seen formula; unseen formula).
- The repo's own census (`docs/research/out/isomer_census.csv`, not re-verified here) reports for `mp_bradley`:
  - 758 formula groups with ≥5 constitutions (6,587 molecules);
  - 40,746 within-formula pairs, of which 2,253 (5.5%) are unlabeled-isomorphic.
  
  So most within-formula variation *is* visible to unlabeled TIs in principle. That makes F a fair test, not a rigged one.

---

## 3. Exact mathematical equivalences under different notation

1. **All edge-additive vertex-degree-based (VDB) indices are linear functionals of one small vector.**
   - Any index of the form Σ_{uv∈E} f(d_u, d_v) equals Σ_{i≤j} m_ij f(i,j), where m_ij is the number of edges joining degree-i and degree-j vertices.
   - Rada's matrix form (Cruz–Espinal–Rada 2024) writes this as trace(φᵀα(G)). Altairi 2026 states the same "edge-partition representation", and the whole Sci Rep/Heliyon genre computes indices from m_ij tables.
   - For hydrogen-suppressed organic graphs (max degree 4) there are at most **10** m_ij classes. TopoQSPR's 30 degree-based indices, and Altairi's 16 "modified bond-based" indices (bond factor d_u+d_v−2 is also a function of (i,j)), are therefore all **deterministic functions of the same ≤10 integers**.
   - TopoQSPR's 49 reverse-degree indices use c_v = Δ − d_v + 1 and are functions of (m_ij, Δ).
   - Consequence: the degree-based and reverse-degree blocks cannot exceed the information in (m_ij, Δ). Any "new" VDB index (a new f) adds no information beyond the m_ij vector, only a different nonlinear view of it.
   - This is **ALREADY KNOWN** (edge-partition method). It is worth stating explicitly in the paper as the reason the TI block is low-dimensional. Optionally, the vector itself (m_ij plus Δ) can be added as a transparent baseline. Not verified in the TopoQSPR code here; check whether all 30+49 indices indeed reduce to f(i,j) sums.
2. **Neighbourhood-degree-sum (NDS) indices** (Mondal ND1–ND6, TopoQSPR's 10) have the same edge-partition form with S_u (the sum of neighbour degrees) in place of d_u. They depend on the NDS edge-partition counts.
3. **Altairi's "modified" Zagreb/Forgotten forms** are f(i,j)/(i+j−2) re-weightings. The edge degree d(e) = d_u + d_v − 2 is the vertex degree in the line graph, as in the reformulated Zagreb index (Miličević, Nikolić & Trinajstić 2004, not re-checked here). Their "bond factor" is the line-graph degree under a new name. Note: equivalence by definition; not checked against each of Altairi's 16 kernels because the formulas are images.
4. **Valence-weighted (labeled) Zagreb indices** already exist under several names: alvaDesc ZM1V, ZM1Kup, ZM1Mad, ZM1Per… (prior audit S8); Mordred "Zagreb"; Barysz/Ivanciuc weighted matrices. A TopoQSPR "labeled family" (E) would most likely coincide with one of these for some weight choice. Check each proposed weight against alvaDesc and Mordred before claiming novelty.

---

## 4. Claim-by-claim detail

### (A) Unlabeled TIs add no detectable gain beyond composition on measured endpoints

- **Prior art:**
  - Altairi 2026 (opposite conclusion, weak design, computed targets);
  - Basak 1997 hierarchy (TS weaker than TC; no composition baseline);
  - Çolakoğlu 2026 and Kara 2025 (size confounding admitted; prior audit);
  - Niane & Djagba 2026 re-test (prior audit);
  - Dearden 2017 chapter, which says TIs are criticised as "inferior to physicochemical descriptors in modelling" ([ABS] via Consensus, doi:10.1007/978-3-319-56850-8_2).
- **Not found:** composition (element counts) as the reference; measured MP/logP/logD/logS; nested repeated group CV with ΔR² uncertainty.
- **Verdict: INCREMENTAL.** The niche is the protocol and the baseline, not the question. Frame it as a replication/re-test that contests Altairi 2026 and the genre.

### (B) Representation-limit facts

- (i) Isomorphic graphs give equal invariants: definitional. **ALREADY KNOWN.** Stating it as a "result" would be criticised.
- (ii) The Bayes risk under squared loss is E[Var(Y|T)] (law of total variance; the regression function is the L²-optimal predictor). Textbook. **ALREADY KNOWN.**
- (iii) An empirical collision floor on real data: Zopf 2022 did the classification analogue on molecular TU sets, with and without node labels. TI regression on measured endpoints was not found. **INCREMENTAL.** It is non-binding on TopoQSPR's data (F9), so present it as a diagnostic.
- The degeneracy and sensitivity literature (Bonchev 1981, Razinger 1985, Konstantinova 1996, Dehmer 2012, Mondal 2021) covers *non-isomorphic* degeneracy. TopoQSPR's finding that all in-data collisions are isomorphism collisions (no non-isomorphic degeneracy among 102 indices) is consistent with the "superindex" result (combining indices removes degeneracy). It is a minor empirical note.

### (C) Incremental value: none found

- The design is known (Basak hierarchy; Altairi "Combined vs Physicochemical").
- The null finding against composition with nested CV was not found. It directly contradicts Altairi 2026.
- **INCREMENTAL.** It is the publishable core together with (A).

### (D) Disease → ChEMBL cohort software

- ChEMBL API `drug_indication` supports `mesh_heading__icontains` and `max_phase_for_ind`. Checked on 2026-10-13: "parkinson", phase 4 returned total_count 47.
- The prior audit found no TI/QSPR tool packaging this.
- **INCREMENTAL** convenience; not a scientific novelty.

### (E) New atom/bond-labeled index family

- Prior constructions:
  - Kier–Hall 1976 (δᵛ);
  - Barysz 1983;
  - Ivanciuc 2000 (Z, X, Y schemes);
  - Basak TC indices;
  - alvaDesc and Mordred labeled variants;
  - Nadjafi-Arani et al. 2025: "degree-distance-based topological indices … derived from vertex edge weighted (VEW) molecular graphs, weighted by atomic number, mass, radius, density, electronegativity, and ionization energy", 166 molecules ([ABS] via Consensus, doi:10.46793/match.95-2.18225).
- **ALREADY KNOWN** as a concept.
- If pursued, it is only defensible as an **ablation**: does labeling close the gap to composition+? Not as "a new family".

### (F) Formula-controlled isomer benchmark on Bradley MP

- Octane and QM9/ISO17 precedents fix the formula, but they are tiny (octanes) or computed (QM9, ISO17).
- **INCREMENTAL, likely NOVEL in specific form:** measured MP, 758 formula groups ≥5, E3a/E3b estimands.
- Caveats:
  - the amendment itself flags MP as outcome-informed;
  - Bradley contains the E2 drugs;
  - MP has a known symmetry dependence (Carnelley's rule, [2H] via a search summary) that unlabeled TIs partly encode.

---

## 5. Query log (2026-10-13)

| # | Tool | Query / URL | Used for |
|---|---|---|---|
| 1 | Europe PMC REST | `DOI:10.1038/s41598-026-40969-7` → PMC13021970 `fullTextXML` | §2.1 [FT] |
| 2 | Crossref works | 10.1021/ci000068y, 10.1021/ci990129w, 10.1007/s40747-020-00262-0, 10.1186/s13321-023-00787-9, 10.1039/C7SC02664A | metadata, licence fields |
| 3 | WebFetch | pubs.acs.org/doi/10.1021/ci000068y | 403 |
| 4 | WebFetch / curl | link.springer.com article + PDF for 10.1007/s40747-020-00262-0 | login redirect / HTML stub |
| 5 | Semantic Scholar API, PubMed efetch 11128100, 10850777 | Ivanciuc abstracts | §2.2–2.3 |
| 6 | Unpaywall API | 10.1007/s40747-020-00262-0 | OA location (publisher only). Note: this request sent the account email as Unpaywall's required `email` parameter. |
| 7 | WebSearch | Mondal Dey De Pal neighbourhood… octane degeneracy | found arXiv:1906.06660 |
| 8 | arxiv.org/pdf/1906.06660 | full preprint (pdftotext) | §2.4 [PP] |
| 9 | WebSearch | Guo Hernandez-Hernandez Ballester scaffold UMAP 2025 | §2.7 |
| 10 | Europe PMC | PMC12153141 fullTextXML | §2.7 [FT] |
| 11 | WebSearch (extended) | topological indices vs molecular formula baseline | §2.10 (none found) |
| 12 | Crossref | Konstantinova, Balaban, Randić, Kier–Hall, Dehmer, Morris, Gütlein, Crusius DOIs; bibliographic queries for Bonchev 1981, Razinger 1985, Bemis–Murcko, Zopf, Ramakrishnan, MOLTOP, Morris 2021, Goodall–Lee, Basak 1997, Brown 2009, Kramer 2012 | §2.8–2.11 |
| 13 | WebSearch | Razinger…; Bonchev…; Barysz… 1983 | bibliographic confirmation ([2H] for Barysz) |
| 14 | WebSearch | WL-indistinguishable decalin bicyclopentyl | §2.9 |
| 15 | WebSearch (extended) | molecular representation collisions … irreducible error | §2.9 |
| 16 | WebSearch (extended) | QSPR baseline atom counts / elemental composition vs TIs | §2.10 (none found) |
| 17 | arXiv API | 2202.10156, 2407.12136, 2511.17078, 2604.27810, 1810.00826, 2509.01254, math/0404512 | abstracts |
| 18 | arxiv.org/pdf/2202.10156v1 | Zopf full preprint | §2.9 [PP] |
| 19 | Consensus | "Do topological indices add predictive information beyond molecular formula or atom counts in QSPR?" | §2.10, §4(E) |
| 20 | WebSearch | constitutional isomers melting point ML benchmark | §2.11 |
| 21 | WebSearch | heavy-atom/MW trivial baseline (Walters) | nothing usable |
| 22 | WebSearch (extended) ×2 | Bayes error upper bound molecular representations; Brown/Muchmore/Hajduk | §2.9 |
| 23 | Europe PMC | PMC5868307 (MoleculeNet), PMC10712068 (SIMPD) fullTextXML | §2.5–2.6 [FT] |
| 24 | Europe PMC | Dehmer 2012 abstract | §2.8 |
| 25 | WebSearch | Basak Gute Grunwald 1997 vapor pressure hierarchical | §2.8 ([2H] numbers) |
| 26 | ChEMBL REST | `drug_indication.json?mesh_heading__icontains=parkinson&max_phase_for_ind=4` | §4(D) |
| 27 | WebSearch (extended) | ML isomer discrimination experimental MP/BP benchmark | §2.11 |

**Not accessed** (and therefore not claimed): ACS full texts (Ivanciuc ×2, Razinger, Konstantinova, Basak, Kier–Hall), the Springer journal version of Mondal 2021, Altairi 2026 SI, the Barysz 1983 book chapter, and the primary source for the decalin/bicyclopentyl example.
