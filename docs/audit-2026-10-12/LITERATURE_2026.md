# TopoQSPR: independent literature and landscape audit (2026-10-12)

Scope: web search (WebSearch, standard and extended), Crossref metadata API, Europe PMC full-text XML, GitHub API (`gh`), PyPI JSON API, arXiv API and official project pages, all queried on 2026-10-12. Repo files were read but not edited (`docs/LITERATURE_GAP_2026.md`, `docs/CHEMGRAPHX_COMPARISON.md`, `docs/PAPER_VS_CHEMGRAPHX.md`, `paper/main.tex`, `paper/refs.bib`).

Labels used in this report:
- **confirmed**: read in a primary source (paper full text, official documentation, source code or package metadata).
- **confirmed absent**: the primary source lists the full feature or input set, and the item is not in it.
- **not found**: not seen in the sources I searched. This does not mean the item is absent.
- **secondhand**: taken from a search-engine snippet or an aggregator, not from the primary page.

---

## 1. Feature / research comparison matrix

### 1a. Matrix

Cell sources are coded by bracketed numbers (S#) and listed in 1b.

| Tool (version, year) | Topological indices: number and kind | Graph entropies | Labeled (atom/bond-weighted) vs unlabeled graphs | Disease → ChEMBL cohort | Nested CV / leakage-safe evaluation | Scaffold split | Web UI | License | Maintenance |
|---|---|---|---|---|---|---|---|---|---|
| **TopoQSPR** (repo `fecb800`, 2026-10-06; authors' claims) | 104 indices: 30 degree-based, 49 reverse-degree, 10 neighbourhood-degree-sum, 11 distance-based, exact Hosoya Z, degree entropy (S1) | 91 bond-wise entropies; distance entropies not implemented (S1) | Unlabeled hydrogen-suppressed graph; the authors state the indices are "blind to element type and bond order" (S1) | Yes: MeSH heading → phase rule → parent and salt handling → manifest (S1) | Yes: repeated nested CV with grouped folds, corrected t, BH/BY (S1) | Yes: Bemis–Murcko-grouped outer folds (S1) | Yes, an exploratory mode; the confirmatory protocol runs via CLI (S1) | MIT (S1, LICENSE) | Pre-release; not public yet (S1) |
| **ChemGraphX** v1.0 (J Cheminform 18:39, 2026; doi:10.1186/s13321-025-01141-x) | 52 degree-based (incl. degree-sum and hybrid) + 24 distance-based; the text also says "150 different TIs" (S2) | 52 degree + 21 distance entropies (S2) | Unlabeled: indices "defined on the basis of degree and distance measures"; .mol/.pdb parsed into graphs; "weighted" means graph-weighted Mostar/Szeged variants, not atom-type weights (S2). Atom-type weighting not found | Confirmed absent: only 4 input methods (.mol, .pdb, adjacency list, adjacency matrix) (S2) | Confirmed absent: index calculator only; no modelling (S2) | Confirmed absent (no modelling) (S2) | Yes: chemgraphx.vercel.app (React + AWS Lambda) (S2) | Code MIT; article CC BY-NC-ND 4.0 (S2, S3) | Repo: 77 commits, last push 2025-08-26 (before acceptance on 2025-12-11); no test suite seen (S3) |
| **Mordred** (original) / **mordredcommunity** 2.0.7 | 1,826 descriptors (1,613 2D). Includes ABC, ABCGG, Zagreb1/2 and modified, Wiener, Wiener polarity, Balaban J, eccentric connectivity, Chi (path/cluster/chain), information content IC/TIC/SIC, Barysz, detour, topological charge (S4, S5). In the descriptor list, confirmed absent: Hosoya Z, Szeged, Sombor, Mostar, PI, Harary, Schultz (S5) | Information-content (neighbourhood-symmetry) indices only; bond-wise entropies of degree-based indices not found (S5) | Both: unlabeled (Zagreb, Wiener) and labeled (Chi "v/d" variants, Barysz matrices weighted by Z, mass, electronegativity) (S5) | Not applicable (calculator); not found | Not applicable (calculator) | Not applicable | Not found (CLI/Python) | BSD-3-Clause (S4, S6) | Original: README says "no longer maintained", last push 2024-02-07. Community fork: 2.0.7 released 2026-01-22, active pushes in Oct 2026 (S4, S6) |
| **PaDEL-Descriptor** 2.21 | 1,875 descriptors (1,444 1D/2D) via CDK, incl. 42 information-content, 20 topological-charge, 11 detour-matrix, path counts (S7). Zagreb/Wiener/eccentric classes are CDK descriptors (secondhand; PaDEL descriptor xls not opened) | Information-content descriptors (S7) | Both (valence chi paths, E-state, ETA) (S7) | Not applicable | Not applicable | Not applicable | Desktop Java GUI + CLI (S7) | "Free for all" freeware; CDK parts under their own licenses (S7) | Last version 2.21, 21 Jul 2014 (S7) |
| **alvaDesc** 3.0 | Vendor says "almost 6,000" descriptors. The topological-indices block has 90 (ZM1/ZM2 incl. valence/Kupchik variants, Wap, Gutman MTI, Schultz MTI, Narumi, Balaban centric, Xt, eccentric connectivity); matrix-based Wiener-/Harary-/Balaban-/Hosoya-like indices; information indices block 51 (S8). Confirmed absent in descriptor list: Szeged, ABC, GA, Sombor, Mostar (S8) | Information indices (51) (S8) | Both (valence-weighted Zagreb etc.) (S8) | Not applicable | Not applicable (alvaModel, a separate product, not examined) | Not applicable | Desktop GUI + CLI, KNIME/Python wrappers (S8) | Proprietary, commercial (S8; Wikipedia secondhand) | Active (3.0 release noted on vendor page) (S8) |
| **RDKit** 2026_09_1 descriptors | GraphDescriptors: BalabanJ, BertzCT, Ipc/AvgIpc, Chi0–4 (n/v), Kappa1–3, HallKierAlpha (S9). Zagreb, Randić-as-such, ABC, Wiener not found in GraphDescriptors | Ipc (information content of characteristic polynomial) (S9) | Both (Chi v) (S9) | Not applicable | Not applicable | Not applicable (no splitter in RDKit core) | None | BSD-3-Clause (S10) | Release 2026_09_1, 2026-10-03 (S10) |
| **QSPRpred** v3.2.1 (J Cheminform 16:128, 2024) | Via DescriptorSet plug-ins: RDKit, Mordred, Mold2, PaDEL, fingerprints (S11) | Only via plug-ins | Inherits from plug-ins | Not found; ChEMBL-based workflows are possible via Papyrus, but no disease-cohort builder was found (S11) | Train/test plus CrossValAssessor; grid/Optuna HPO uses a ModelAssessor on training data. Built-in nested-CV routine not found (S11) | Yes: Random/Cluster/Scaffold/Temporal splits via BalanceSplit (S11) | Not found (CLI + Python API) (S12) | MIT (S12) | Last tagged release 2024-10-02; commits to 2026-09-11 (S12) |
| **DOPtools** 1.3.9 (Digital Discovery 4:1188, 2025) | CircuS, ChyLine fragment descriptors, fingerprints; Mordred calculator **removed in 1.3.7 (2025-06-19)** (S13) | Not found | Labeled (atom-typed fragments) (S13) | Not found | Optuna trials scored by repeated KFold `cross_val_predict` on the same data; no outer loop in `optimizer.py` (S14). Non-nested per code reading | Not found in optimizer code (S14) | Not found (CLI) | LGPL-3.0 (S13) | Last PyPI release 2025-06-23; last push 2026-05-25 (S13) |
| **QSAR in the Browser (QITB)** (JCIM 66:7805, 2026; doi:10.1021/acs.jcim.6c01010) | None described; Morgan, MACCS, RDKit fingerprints; "physicochemical descriptors may be added" for display (S15) | Not found | n/a (fingerprints) | ChEMBL **target**-based bioactivity fetch, not disease-based (S15) | Single 10-fold CV of RF/XGBoost with user-set hyperparameters; no nesting described (S15) | Scaffold *analysis* (scaffold tree, MMS) but no scaffold *split* described (S15) | Yes: client-side static app (RDKit-MinimalLib WASM + Pyodide) (S15) | Code AGPL-3.0; article CC BY 4.0 (S15, S16) | Repo created 2023-12-30, last push 2026-09-09 (S16) |
| **RANQSAR** (J Comput Aided Mol Des 40:169, 2026; doi:10.1007/s10822-026-00881-8) | Abstract: Morgan fingerprints + RDKit descriptors (secondhand, abstract only; S17) | Not found | Not found | Not found | Abstract: CV, external validation, Y-randomization, AD (secondhand) | Not found | Desktop app (abstract) | "Open-source" per title; license not found | Not open access; no repository found by me or by the authors' audit |
| **QSARtuna** (JCIM 64:5365, 2024; doi:10.1021/acs.jcim.4c00457) | ECFP, MACCS, etc.; README shows fingerprints (S18) | Not found | n/a | Not found | Optuna HPO over algorithm × descriptor with k-fold CV on training data, then build/test. Nesting not stated in README (S18) | Not found in README | Not found | Repo license field empty on GitHub API (S18) | Last push 2026-06-30 (S18) |
| **scikit-fingerprints** (SoftwareX 28:101944, 2024) | `skfp.descriptors.topological`: average Wiener, Balaban J, Burden, diameter, graph distance index, Hall–Kier alpha, Petitjean, polarity number, Kappa1–3, radius, Wiener, Zagreb M1/M2 (S19) | Not found | Both (Burden/Hall–Kier labeled; Wiener/Zagreb unlabeled) | Not found | sklearn-compatible grid/random search (nesting is up to the user) (S19) | Yes: scaffold, randomized scaffold, Butina, MaxMin, PubChem splits (S19) | None | MIT (S19) | Active (push 2026-10-11) (S19) |
| **PyTIE** (arXiv 2510.26813, 2025) | Degree-based and neighbourhood-degree-sum indices from edge partitions; symbolic expressions (S20) | Shannon entropy (S20) | Unlabeled (edge partitions) (S20) | Not found | Not found (calculator) | Not found | Not found | "Open-source" (abstract); license not found | Preprint only (abstract read; full text not read) |
| GitHub one-study pipelines (2025–2026): HakanKucuk38/Topological-Indices-PhysChem; JKMMaths/TIsCalculator; avimallick/topoindex; aliberkanbektas/drug-indices-project_1 | 24 vertex/edge-weighted indices + RF/XGB/NN + SHAP (Kucuk); degree-based from SMILES (TIsCalculator, Streamlit); Wiener/Zagreb/hyper-Wiener (topoindex); PubChem name → indices (Bektas) (S21) | Not found | Kucuk: "vertex- and edge-weighted" (labeled); the others unlabeled | Bektas: drug name → PubChem, no disease/ChEMBL cohort | Kucuk: "cross-validation scores" mentioned; nesting not found | Not found | TIsCalculator: Streamlit | topoindex and Bektas MIT; others no license | Single-author scripts, 0–5 stars (S21) |

### 1b. Sources for the matrix

- S1: TopoQSPR repo: `paper/main.tex` §2.3–2.5 and Table 2; `docs/LITERATURE_GAP_2026.md` §4; `LICENSE`. These are the authors' own claims, which I did not re-execute.
- S2: ChemGraphX full text, Europe PMC PMC13020044 (fetched `fullTextXML`). Code-metadata table, Table 1 ("52 degree-based indices and 52 degree-based entropies; 24 distance-based indices and 21 distance-based entropies"), Introduction, Conclusion (RAM limit). Received 2025-08-26, accepted 2025-12-11, published 2026-02-11 (Crossref).
- S3: `gh api repos/aditya-jindal/chemgraphx` (MIT, pushed_at 2025-08-26) and https://github.com/aditya-jindal/chemgraphx (77 commits; top level has client, server, LICENSE, README).
- S4: https://github.com/JacksonBurns/mordred-community (BSD-3, 1,826 descriptors) and https://github.com/mordred-descriptor/mordred (last push 2024-02-07).
- S5: https://mordred-descriptor.github.io/documentation/master/descriptors.html (full page grepped).
- S6: https://pypi.org/project/mordredcommunity/ (2.0.7, 2026-01-22).
- S7: http://www.yapcwsoft.com/dd/padeldescriptor/ (fetched; "1875 descriptors", "Version 2.21 [21 Jul 2014]", licence text). Yap 2011, doi:10.1002/jcc.21707.
- S8: https://www.alvascience.com/alvadesc-descriptors/ (full list grepped) and https://www.alvascience.com/alvadesc/ ("alvaDesc 3.0"). Proprietary per https://en.wikipedia.org/wiki/AlvaDesc (secondhand).
- S9: https://www.rdkit.org/docs/source/rdkit.Chem.GraphDescriptors.html
- S10: `gh api repos/rdkit/rdkit/releases/latest` (Release_2026_09_1, 2026-10-03; BSD-3).
- S11: QSPRpred full text, Europe PMC PMC11566221 (sections "Data splitting", "Hyperparameter optimization", "Descriptors"). doi:10.1186/s13321-024-00908-y.
- S12: https://github.com/CDDLeiden/QSPRpred and the GitHub API (MIT, v3.2.1 2024-10-02, pushed 2026-09-11).
- S13: https://pypi.org/pypi/DOPtools/json and README at https://github.com/POSidorov/DOPtools (changelog: "Mordred ... has been removed" in 1.3.7, 2025-06-19; LGPL-3.0).
- S14: `doptools/optimizer/optimizer.py` at HEAD (lines 135–252: `cross_val_predict(model, X, Y, cv=KFold(cv_splits))` over `cv_repeats`, inside the Optuna objective).
- S15: QITB full text, Europe PMC PMC13417879 ("QITB evaluates the model via 10-fold cross validation"; fingerprint options; client-side architecture; GitHub URL; CC BY 4.0).
- S16: `gh api repos/syedzayyan/qsar-in-browser` (AGPL-3.0; created 2023-12-30; pushed 2026-09-09).
- S17: Crossref record for doi:10.1007/s10822-026-00881-8 and a WebSearch snippet of https://link.springer.com/article/10.1007/s10822-026-00881-8 (abstract only; Europe PMC isOpenAccess = N).
- S18: https://github.com/MolecularAI/QSARtuna README (via gh). Paper doi:10.1021/acs.jcim.4c00457 (not OA, not read).
- S19: https://github.com/MLCIL/scikit-fingerprints (`skfp/descriptors/topological.py`, `skfp/model_selection/splitters/`); SoftwareX doi:10.1016/j.softx.2024.101944.
- S20: arXiv API abstract, arXiv:2510.26813 (2025-10-27).
- S21: GitHub API READMEs, queried 2026-10-12.

### 1c. Observations relevant to the authors' claims

1. **"No tool builds disease-anchored ChEMBL cohorts."** I found no such feature in the tools above, in a GitHub search for "topological indices" (30 most recently updated repos), or in WebSearch for 2025–2026 TI-QSPR software. Scope limit: I did not search PyPI exhaustively, ChemRxiv, or non-English venues. Write it as "not found in a search of …", not "no tool does X". QITB does ChEMBL retrieval, but by target, not by disease (S15).
2. **The nested CV claim needs care.** QSPRpred provides leakage-aware splitters (scaffold/cluster/temporal) and HPO on training folds with a separate test assessor (S11). That is leakage-safe if used as documented, even though it is not packaged as "repeated nested CV". The TopoQSPR comparison table (Table 3, `main.tex` l.177–193) does not list a "leakage-safe evaluation" column, which is fair. But the text "none evaluates index-based models against size and composition baselines by default" is narrow enough to be true. The paper should not imply that QSPRpred lacks leakage-safe evaluation.
3. **DOPtools footnote is correct.** Mordred removal in 1.3.7 (2025-06-19) is confirmed (S13). DOPtools' optimizer reports the same repeated-CV score it optimises, which is non-nested (S14). If the paper says this, it should be worded as a code reading at HEAD, not as a claim about the published article.
4. **ChemGraphX counts are internally inconsistent in the source.** Table 1 lists 52+24 indices and 52+21 entropies (149 values), while the text says "ChemGraphX computes 150 different TIs" (S2). The authors' docs already note the 149/150 discrepancy. Good.
5. **The ChemGraphX article contains an error about PaDEL.** It says PaDEL "does not support graph-theoretical indices" (S2). PaDEL's own page lists 42 information-content, 20 topological-charge and 11 detour-matrix descriptors, and chi path indices (S7). Mention this only if needed. It concerns the senior co-author's paper.
6. **Labeled-graph blind spot.** Mordred, PaDEL, alvaDesc, RDKit and scikit-fingerprints all include valence- or atom-weighted variants (S5, S7, S8, S9, S19). TopoQSPR and ChemGraphX do not. Reviewers will see this as the obvious reason TIs lose to "composition" in TopoQSPR's own ladder (`main.tex` l.202: composition > TI on every endpoint). Pre-empt it by citing Kier–Hall valence connectivity (already cited) and stating the scope.
7. **Hosoya Z, Szeged, Mostar, PI and Sombor are absent from Mordred and from alvaDesc's lists** (S5, S8; alvaDesc has only "Hosoya-like" matrix indices). The claim that TopoQSPR implements indices that common calculators lack holds against these two. ChemGraphX has Szeged/Mostar/PI (S2).

---

## 2. Modern benchmark practice (verified citations)

All DOIs below were resolved through Crossref on 2026-10-12 unless marked otherwise.

### 2a. Split strategy and its limits

- **Scaffold splits are still optimistic.** Guo, Hernandez-Hernandez & Ballester, "Scaffold Splits Overestimate Virtual Screening Performance", ICANN 2024, LNCS pp. 58–72, doi:10.1007/978-3-031-72359-9_5. Extended in "UMAP-based clustering split for rigorous evaluation of AI models for virtual screening on cancer cell lines", *J Cheminform* 17:94 (2025), doi:10.1186/s13321-025-01039-8 (OA, PMC12153141). Difficulty ranking: UMAP > Butina > scaffold ≈ random (search snippet; abstract-level). A 2026 survey (arXiv:2604.16586, secondhand) relays Walters' caution that UMAP splits give imbalanced clusters.
- **Time and simulated-time splits.** Sheridan, *JCIM* 53:783 (2013), doi:10.1021/ci400084k. Landrum et al., "SIMPD", *J Cheminform* 15:119 (2023), doi:10.1186/s13321-023-00787-9.
- **Leakage-aware splitting.** Joeres et al., "Data splitting to avoid information leakage with DataSAIL", *Nat Commun* 16:3337 (2025), doi:10.1038/s41467-025-58606-8. Tossou et al., "Real-World Molecular Out-Of-Distribution", *JCIM* 64:697 (2024), doi:10.1021/acs.jcim.3c01774. General ML: Kapoor & Narayanan, *Patterns* 4:100804 (2023), doi:10.1016/j.patter.2023.100804.
- **Implication for TopoQSPR.** Bemis–Murcko-grouped outer folds are a reasonable stricter check, but cite them as "still optimistic". The paper already does this (`main.tex` l.144). Approved drugs from many disease areas are structurally diverse, so random and scaffold CV may differ less than on congeneric sets. Report the observed random-vs-scaffold gap and avoid assuming it.

### 2b. MoleculeNet / TDC critiques

- MoleculeNet: Wu et al., *Chem Sci* 9:513 (2018), doi:10.1039/C7SC02664A.
- TDC: Huang et al., NeurIPS 2021 Datasets & Benchmarks, arXiv:2102.09548 (not resolved via Crossref; arXiv ID only).
- Walters, "We Need Better Benchmarks for Machine Learning in Drug Discovery", Practical Cheminformatics blog, 2023-08-03 (read in full). Points: invalid SMILES and duplicates with conflicting labels in BBBP; about 71% undefined stereocentres in BACE; mixed-assay IC50s; **ESOL spans >13 log units while pharma solubility assays span ~2.5–3, so wide dynamic range flatters simple models**; no agreed splits. This bears directly on TopoQSPR's AqSolDB endpoint (AqSolDB: Sorkun et al., *Sci Data* 6:143 (2019), doi:10.1038/s41597-019-0151-1), which is also wide-range. Report MAE/RMSE alongside R² and discuss dynamic range.
- Landrum & Riniker, "Combining IC50 or Ki Values from Different Sources Is a Source of Significant Noise", *JCIM* 64:1560 (2024), doi:10.1021/acs.jcim.4c00049.
- Wognum, Ash et al., "A call for an industry-led initiative to critically assess ML for real-world drug discovery", *Nat Mach Intell* 6:1120 (2024), doi:10.1038/s42256-024-00911-w.
- Steshin, "Lo-Hi: Practical ML Drug Discovery Benchmark", NeurIPS 2023 D&B (search result only; not opened).

### 2c. Method-comparison protocol

- **Ash, Wognum, Rodríguez-Pérez, Aldeghi, Cheng, Clevert, Engkvist, Fang, Price, Hughes-Oliver, Walters, "Practically Significant Method Comparison Protocols for Machine Learning in Small Molecule Drug Discovery", *JCIM* 65(18):9398–9411 (2025), doi:10.1021/acs.jcim.5c01609** (preprint doi:10.26434/chemrxiv-2024-6dbwv-v2). The paper is not OA and I did not read it. From the authors' official repo (github.com/polaris-hub/polaris-method-comparison, MIT; notebook `ADME_example/ML_Regression_Comparison.ipynb`):
  - 5×5 repeated CV as the sampling mechanism (for about 2k compounds);
  - repeated-measures ANOVA, then **Tukey HSD**;
  - effect-size-aware "MCSim" plots and CIs on differences;
  - an explicit "Don't do this" against bold-the-best mean tables;
  - a dynamic-range example;
  - scaffold splits used in the example.
  
  **Implication.** TopoQSPR uses a corrected repeated k-fold t-test (Nadeau–Bengio/Bouckaert) with BH. That is defensible but differs from Ash et al.'s recommended default. Say why: pairwise ΔR² against one reference is not all-pairs. Also report a practical-significance threshold, which TopoQSPR appears to have ("practical threshold", `main.tex` l.202). The citation year in the user's prompt ("2024/2025?") resolves to 2025 (journal); the preprint is 2024.
- Bender et al., "Evaluation guidelines for machine learning tools in the chemical sciences", *Nat Rev Chem* 6:428 (2022), doi:10.1038/s41570-022-00391-9.
- Artrith et al., "Best practices in machine learning for chemistry", *Nat Chem* 13:505 (2021), doi:10.1038/s41557-021-00716-z. The repo's `refs.bib` key `artrith2021best` should be checked for the 505–508 page range.
- Classic QSAR validation: Tropsha, *Mol Inform* 29:476 (2010), doi:10.1002/minf.201000061; Golbraikh & Tropsha, "Beware of q2!", *J Mol Graph Model* 20:269 (2002), doi:10.1016/S1093-3263(01)00123-1; Alexander, Tropsha & Winkler, "Beware of R2", *JCIM* 55:1316 (2015), doi:10.1021/acs.jcim.5b00206; Hawkins, "The Problem of Overfitting", *JCICS* 44:1 (2004), doi:10.1021/ci0342472.

### 2d. Recommended baselines

- **ECFP + RF/GBM and RDKit 2D descriptors are strong baselines.**
  - Deng et al., "A systematic study of key elements underlying molecular property prediction", *Nat Commun* 14:6395 (2023), doi:10.1038/s41467-023-41948-6.
  - Praski, Adamczyk & Czech, arXiv:2508.06199 (v4, Feb 2026; preprint): 25 pretrained embedding models on 25 datasets, "nearly all neural models show negligible or no improvement over the baseline ECFP".
  - van Tilborg et al., "Exposing the Limitations of Molecular Machine Learning with Activity Cliffs", *JCIM* 62:5938 (2022), doi:10.1021/acs.jcim.2c01073 (descriptor/fingerprint ML beats deep learning on cliffs).
  - Fang et al., Biogen prospective ADME validation, *JCIM* 63:3263 (2023), doi:10.1021/acs.jcim.3c00160.
- **Chemprop / D-MPNN.**
  - Yang et al., *JCIM* 59:3370 (2019), doi:10.1021/acs.jcim.9b00237.
  - Heid et al., "Chemprop: A Machine Learning Package for Chemical Property Prediction", *JCIM* 64:9 (2024), doi:10.1021/acs.jcim.3c01250. Current PyPI chemprop 2.3.1, MIT.
  - Low-data foundation model: Burns et al., "Deep Learning Foundation Models for Low-Data Regimes from Classical Molecular Descriptors" (CheMeleon), *JCIM* 66:10692 (2026), doi:10.1021/acs.jcim.6c01546 (not OA; not read).
- **Implication.** TopoQSPR's ladder (null → size → composition → TI → RDKit 2D) covers the descriptor baseline. It omits **ECFP + RF/GBM**, which is the single most-expected baseline in 2023–2026 practice and is cheap to add. A learned baseline (Chemprop/CheMeleon) is deferred in the repo, which is defensible for a software paper. ECFP is not.
- Other 2026 practice paper cited in the repo, verified: Alshehade et al., "Systematic Multi-Target QSAR Benchmarking", *JCIM* 66:8923 (2026), doi:10.1021/acs.jcim.6c01237 (not read).

---

## 3. Criticisms and replications of small-sample TI-QSPR

### 3a. Index-quality and proliferation critiques (chemical graph theory community)

- Gutman & Tošović, "Testing the quality of molecular structure descriptors. Vertex-degree-based topological indices", *J Serb Chem Soc* 78:805–810 (2013), doi:10.2298/JSC121002134G. Twenty degree-based indices against octane-isomer heat of formation and boiling point; "correlation ability of many of these indices is either rather weak or nil" (search snippet). Note: octane isomers have fixed size, so this tests discrimination, not size confounding (the repo already notes this).
- Furtula, Gutman & Dehmer, "On structure-sensitivity of degree-based topological indices", *Appl Math Comput* 219:8973 (2013), doi:10.1016/j.amc.2013.03.072.
- Gutman, "Degree-Based Topological Indices", *Croat Chem Acta* 86:351 (2013), doi:10.5562/cca2294 (review).
- **Vukičević, Li, Sedlar & Došlić, "Lanzhou Index", *MATCH Commun Math Comput Chem* 80:863–876 (2018)** (read p. 1–2 of the PDF). The introduction opens with the proliferation critique ("hundreds and thousands of them are clogging the literature"). It reports Gutman's criterion that new indices deserve study only if they outperform the best known indices on a chemically relevant property, and refers to Randić's twelve requirements for new descriptors. No DOI found (MATCH has no DOIs for that volume; URL: https://match.pmf.kg.ac.rs/electronic_versions/Match80/n3/match80n3_863-876.pdf).
- Redžepović, "Chemical applicability of Sombor indices", *J Serb Chem Soc* 86:445 (2021), doi:10.2298/JSC201215006R (benchmark-style test on octanes/benzenoids).
- Hollas, "The Covariance of Topological Indices that Depend on the Degree of a Vertex", arXiv:math/0404512 (2004). Correlation between edge-additive indices is driven by the edge set (that is, size). This is the closest formal statement of size domination for degree-based indices I found.
- Size dependence: Gao & Hosoya, "Topological Index and Thermodynamic Properties. IV. Size Dependency of the Structure-Activity Correlation of Alkanes", *Bull Chem Soc Jpn* 61:3093 (1988), doi:10.1246/bcsj.61.3093 (metadata only). Basak, Magnuson, Niemi et al., *Math Model* 8:300 (1987), doi:10.1016/0270-0255(87)90594-X. Rücker & Rücker, "On Topological Indices, Boiling Points, and Cycloalkanes", *JCICS* 39:788 (1999), doi:10.1021/ci9900175.
- Debate context: Randić, "On History of the Randić Index and Emerging Hostility toward Chemical Graph Theory", *MATCH* 59:5–124 (2008) (search snippet; a defence that documents long-running criticism of connectivity indices).

### 3b. QSAR-general pitfalls (multiple testing, chance correlation)

- Topliss & Edwards, "Chance factors in studies of QSAR", *J Med Chem* 22:1238 (1979), doi:10.1021/jm00196a017.
- Rücker, Rücker & Meringer, "y-Randomization and Its Variants in QSPR/QSAR", *JCIM* 47:2345 (2007), doi:10.1021/ci700157b.
- Dearden, Cronin & Kaiser, "How not to develop a QSAR/QSPR relationship", *SAR QSAR Environ Res* 20:241 (2009), doi:10.1080/10629360902949567.
- Cronin & Schultz, "Pitfalls in QSAR", *J Mol Struct THEOCHEM* 622:39 (2003), doi:10.1016/S0166-1280(02)00616-4.

### 3c. Recent genre papers that state the problems themselves (primary full text read)

- **Çolakoğlu & Benjamin, *ACS Omega* 11:53404 (2026), doi:10.1021/acsomega.6c06753** (PMC13563987). n = 18 bipolar drugs. Properties are "obtained from PubChem and ChemSpider" (computed or predicted values, not measured). RF/XGBoost deteriorated under LOOCV. The authors state that MW and HAC "are strongly related to molecular size … their excellent predictive performance should be interpreted primarily as reflecting descriptor consistency with molecular size". HAC LOOCV R² = 0.9994. This is the clearest in-genre admission of size confounding and of computed targets.
- **Kara et al., *Sci Rep* 16:3005 (2025/2026), doi:10.1038/s41598-025-32906-x** (PMC12830717). n = 30 anticancer drugs; properties from PubChem and ChemSpider; BSO₁ vs HAC R = 0.999; limitations section concedes the small dataset.
- Almulla et al., *Sci Rep* (2026), doi:10.1038/s41598-026-61020-9 (OA, PMC13483912; Crossref rate-limited, so DOI confirmed via Europe PMC only).
- Other 2026 genre hits (WebSearch only, not read): SARS-CoV-2/TB distance-based indices with 8 and 13 drugs, negative LOOCV Q² (*Sci Rep*, doi:10.1038/s41598-026-49132-8, secondhand); SSRI/SNRI antidepressants (*Sci Rep*, doi:10.1038/s41598-026-65586-2, secondhand).

### 3d. Re-tests and replications

- Niane & Djagba, "Graph-Theoretic Models for the Prediction of Molecular Measurements", arXiv:2604.19840 (2026-04-21; preprint). Re-test of the Mukwembi–Nyabadza D(G)–ζ(G) polynomial model, originally fit on a small flavonoid set, on five MoleculeNet sets: average R² = 0.24. Enhanced models with Morgan fingerprints and physicochemical features reach about 0.79 (abstract read). This is the only external re-test of a small-sample graph-index model on larger data that I found.
- **No published replication, Matters Arising, or meta-analysis specific to drug TI-QSPR (n ≈ 20–30) was found** (WebSearch for Matters Arising/comment on TI-QSPR in Sci Rep; Europe PMC not exhaustively queried). TopoQSPR's Experiment D would be, as far as I found, the first systematic re-analysis. Frame it as "not found in a search of …".

### 3e. Retractions in the genre (integrity context, not methodological critique)

- "RETRACTED: QSPR Predictions of SARS-CoV-2 Drug Materials via Degree-Based Topological Indices", *EPJ Web Conf* 359:01021 (2026), doi:10.1051/epjconf/202635901021. The "RETRACTED" title prefix is confirmed in Crossref. The reason was not readable (publisher page returned 403).
- "Retracted: On Vertex Degree-Based Topological Indices for Fixed Branching Vertices of Trees", *J Chem* (Hindawi) 2023, doi:10.1155/2023/9836740. This sits within Hindawi's 2023 mass retraction of >8,000 paper-mill articles (Retraction Watch, 2023-12-19, secondhand). The specific reason was not read.
- Use these sparingly, if at all. They do not bear on TopoQSPR's argument and could read as hostile to the genre the co-authors work in.

---

## 4. Journal fit

Official pages for all three journals (ScienceDirect/Elsevier, link.springer.com, rsc.org) returned HTTP 403 or a JavaScript bot challenge to both WebFetch and curl. Chrome automation needed a browser choice that only the user can make, so it was not used. Findings below therefore come from (i) other official artefacts I could read, (ii) official pages seen through search-engine snippets (secondhand), or (iii) aggregators (flagged).

### 4a. SoftwareX (Elsevier; ISSN 2352-7110)

**Article types.** Original Software Publication (OSP) and Software Update.
- Source: the official Elsevier "SoftwareX software-update article template Version 6 (March 2026)", https://legacyfileshare.elsevier.com/promis_misc/softwarex-software-update-template.docx (read in full).
- The template says OSPs must use the OSP template "found on the Guide for Authors" and that "Any submission that deviates widely will likely face desk rejection."

**Length.** Unresolved conflict:
- the SoftwareX aims page (via an aggregator mirror of the aims text, myhuiban.com/journal/711.md) says "a short descriptive paper of 3000-word limit";
- a WebSearch summary of the Guide for Authors says 4,000 words (not verifiable in the primary source);
- the repo's own `docs/LITERATURE_GAP_2026.md` says 4,000, also from secondary sources;
- the 2015 author pack said 6 pages.

My crude count of `paper/main.tex` body text, excluding tables, figures and references, is about 4,700 words. **That is over either limit.** Verify manually in a browser.

**Code requirements** (official template, confirmed):
- "We only accept GitHub repositories within Metadata Tables … required for forking to our SoftwareX repository". The ElsevierSoftwareX GitHub org has 2,539 forked repos, most recently created 2026-09-28, which confirms the practice is current.
- "GitHub repositories must have a well-documented README.md file, and a Licence.txt file."
- Code metadata table C1–C8.
- License from a pre-approved list (Apache, GPL, MIT …). MIT is fine.

**APC.** USD 760 per DOAJ/LetPub/Researcher.Life (secondhand). One French university guide says USD 970 (secondhand). Confirm on ScienceDirect.

**Topical fit.** Chemistry and cheminformatics descriptor software is published there:
- scikit-fingerprints, *SoftwareX* 28:101944 (2024), doi:10.1016/j.softx.2024.101944;
- CombinatorixPy (mixture descriptors), 29:102060 (2025), doi:10.1016/j.softx.2025.102060;
- ODAMNet (chemicals → rare diseases), doi:10.1016/j.softx.2024.101701.

Impact factor about 1.9 (2025 JCR; secondhand, journalmetrics.org).

**For.**
- Format matches a workflow tool with a code metadata table.
- MIT plus GitHub is already in place.
- There is precedent for descriptor software in the journal.
- The low APC.
- Does not require a "significant advance over existing software" in J Cheminform's sense.
- Avoids publishing in the same journal and year as ChemGraphX.

**Against.**
- The strict word limit forces Experiments A–D, which are the paper's scientific substance, into the supplement.
- Reviewers judge software impact, not the methodological finding. The "TIs don't beat size on 3/4 endpoints" result would get less weight.
- GitHub-only, plus the forking requirement (fine).
- The repo must be public with README.md and **Licence.txt**. Note that the repo has `LICENSE`, not `Licence.txt`. Check whether the name matters; the template wording is literal.
- The repo's "purge of `related_papers/`" must happen before the code is public.

### 4b. Journal of Cheminformatics (Springer Nature/BMC)

**Article types.** Research, Methodology, Software, Database, Data Note, Brief Report, Review, Educational, Comment and others (WebSearch snippet of the official guidelines; secondhand).

**Software-article requirements** (journal-specific guidelines maintained in the journal's GitHub org, https://jcheminform.github.io/jcheminform-author-guidelines/guidelines/software.html; last commit to that repo 2023-01-11, so possibly stale):
- the tool must be "likely to be of broad utility" and represent "a significant advance over previously published software (usually demonstrated by direct comparison with available related software)";
- OSI-approved license with no registration barrier;
- source code archive included at submission;
- reviewer-anonymous access to the software;
- Zenodo archiving encouraged;
- web links archived in the Internet Archive;
- "will only publish research or software that is entirely reproducible by third parties".

**Scientific Contribution statement.** Mandatory abstract section, at most 3 sentences: Zdrazil et al., "Are new ideas harder to find? …", *J Cheminform* 16:6 (2024), doi:10.1186/s13321-023-00798-6. The TopoQSPR abstract already has one, apparently drafted with J Cheminform in mind.

**Length.** No word limit found for Software articles (BMC style).

**APC.** £1,690 / $2,390 / €1,990 (official fees page via search snippet; secondhand). University of Melbourne agreement listing, updated 2026-07-04, gives EUR 1,990 (secondhand).

**For.**
- Best topical fit: QSPR methodology plus open software plus reproducibility.
- No tight word limit, so Experiments A–D can stay in the main text.
- The "Scientific Contribution" statement is already written.
- The readership (Walters/Landrum-style evaluation culture) values negative and size-baseline results.
- Direct comparison with related software is expected, and the repo has done it thoroughly.

**Against.**
- ChemGraphX, from the senior co-author's group, appeared there in Feb 2026, so the editors may see overlap. TopoQSPR must be clearly a workflow and evaluation contribution, not "more indices".
- The "significant advance" bar for a Software article is high. The more natural article type may be **Research or Methodology**, with the software as the vehicle.
- Highest APC of the three in USD.
- Reviewers will ask for ECFP and Chemprop baselines and for Ash et al.-style comparisons.
- The unlabeled-graph limitation will be criticised.

### 4c. Digital Discovery (RSC)

**Article types.** Communications (normally ≤4 printed pages), Full papers (no page limit; "appropriateness of length to content will be taken into consideration"), Commit (technical notes), Reviews, Perspectives, Opinions, Comments, Tutorial Reviews. Source: WebSearch snippets of https://www.rsc.org/publishing/publish-with-us/publish-a-journal-article/digital-discovery (secondhand; direct fetch 403).

**Data and code.**
- Code and data deposited in repositories and accessible to referees during review; custom code in a persistent repository (Zenodo, Code Ocean …) on acceptance (official guidelines via snippet; secondhand).
- An additional **data reviewer** assesses data and code for original research (RSC blog, 2021-12-16, https://blogs.rsc.org/dd/2021/12/16/data-reviewer/; secondhand via snippet).
- An explicit open-source license mandate was not found.

**APC.** £2,200, gold OA (Digital Discovery newsletter, 2026-04-15, https://blogs.rsc.org/dd/2026/04/15/q2-2026-news/, read; RSC payments page via snippet). The same newsletter gives IF 5.6 (2024) and a 67% acceptance rate.

**Topical fit.** ML/AI/automation for chemistry. DOPtools (descriptor + optimisation software) was published there as a Full paper (doi:10.1039/d4dd00399c).

**For.**
- Software plus benchmarking papers fit, with DOPtools as a direct precedent.
- No page limit for Full papers.
- The data reviewer will value the run records, frozen plan and regenerable tables.
- High acceptance rate and IF.
- Avoids the ChemGraphX journal overlap.

**Against.**
- The audience expects ML baselines (ECFP+GBM, Chemprop). A study whose learners are classical regressors on 104 unlabeled-graph indices may look dated unless framed as an evaluation and replication study.
- The APC is near J Cheminform's.
- The disease-cohort feature is less central to that readership.

### 4d. Recommendation

**First choice: Journal of Cheminformatics, as a Research or Methodology article ("Software" only if the editors prefer).** Lead with the scientific result: on curated approved-drug cohorts with nested, scaffold-grouped CV, topological indices rarely beat molecular size and never beat element composition. Published TI-QSPR R² values collapse under the genre's own procedure. TopoQSPR is the reproducible instrument that shows this. This uses the paper's strongest asset, Experiments B–D, which SoftwareX's word limit would push into the supplement.

Conditions before submitting:
1. Add an ECFP + RF/GBM rung to the ladder.
2. Cite Ash et al. 2025 and justify the corrected-t + BH choice against RM-ANOVA + Tukey.
3. Address the unlabeled-graph limitation explicitly.
4. Make the ChemGraphX relationship and author overlap transparent in the cover letter.

**Fallback: Digital Discovery** (Full paper), if J Cheminform declines on overlap with ChemGraphX. **SoftwareX** only if the authors want a short tool paper and are willing to cut to about 3,000–4,000 words. In that case Experiments A–D must move to a separate methods paper, which would be a reasonable two-paper plan.

---

## 5. Access limitations

**Could not read in full** (403, bot challenge, or paywall):
- SoftwareX Guide for Authors and OSP template;
- J Cheminform official Springer guideline and fees pages;
- RSC Digital Discovery author guidelines;
- the DOPtools article page (RSC 403; facts taken from PyPI/GitHub instead);
- the JCIM landing pages for Ash et al. 2025, QSARtuna 2024, CheMeleon 2026 and Alshehade 2026 (not OA in Europe PMC);
- RANQSAR (JCAMD, not OA; abstract via snippet only);
- the EPJ Web Conf retraction notice;
- Lo-Hi and TDC papers (citations only);
- PaDEL's descriptor spreadsheet (not downloaded).

**Secondhand** (search snippets or aggregators):
- SoftwareX word limit (3,000 vs 4,000 conflict unresolved) and APC (USD 760 vs 970);
- J Cheminform APC and article-type list;
- Digital Discovery article types, length rules and code/data policy wording;
- RANQSAR features;
- Hindawi mass-retraction context;
- the UMAP-split critique attributed to Walters (via a 2026 survey);
- alvaDesc being proprietary (Wikipedia, consistent with the vendor's licence-request flow);
- PaDEL's Zagreb/Wiener classes (inferred from CDK).

**Read in full or in the relevant primary sections:**
- ChemGraphX, QSPRpred, QITB, Çolakoğlu 2026 and Kara 2025 full texts (Europe PMC XML);
- the Walters 2023 blog post;
- the Vukičević et al. 2018 introduction;
- the SoftwareX update template (March 2026);
- J Cheminform journal-specific software guideline (GitHub-hosted, last updated 2023);
- Mordred, alvaDesc and RDKit descriptor lists;
- DOPtools/QSPRpred/QITB/ChemGraphX/QSARtuna/scikit-fingerprints repository metadata and code;
- the Ash et al. official repository notebook.

**Not verified:**
- that TopoQSPR's own numbers (104 indices, test counts, experiment results) are correct. This audit takes the repo's text as the authors' claims.

**Tooling note.** Claude in Chrome was available but needed the user to pick between two connected browsers. As a subagent I could not ask, so the journal pages behind bot challenges stay unverified. A human should open the three guideline pages and confirm the word limit, APC and current templates before choosing a venue.
