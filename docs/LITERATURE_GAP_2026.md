# Literature and software gap (October 2026)

Sources: the landscape audit (`docs/audit/expert-reports-2026-10.json`), the skeptic's re-tests,
and the lead author's own checks of ChemGraphX (Crossref + Europe PMC full text, 2026-10-02).
Items marked † come from the landscape audit only and were **not** independently re-verified.

## 1. What TopoQSPR is not new for

| Tool | What it does | Inputs | Topological indices | QSPR / ML | Disease → drug set | Nested CV + size baseline | Provenance |
|---|---|---|---|---|---|---|---|
| **ChemGraphX** (Jacob, Clement, Arockiaraj, Jindal, *J Cheminform* 18:39, 2026; MIT; github.com/aditya-jindal/chemgraphx) | Web tool for topological indices and graph entropies | .pdb, .mol, adjacency list/matrix (no SMILES) | ~150 per the text (its own table sums to 149): degree, degree-sum, hybrid, distance, entropies | no | no | no | no |
| Mordred / mordredcommunity 2.0.7 | General descriptor calculator | SMILES/SDF | Zagreb, Randić, ABC, Wiener, Balaban J and more (matches TopoQSPR exactly on these, verified on 5 drugs); no Sombor, GA, ISI, reverse-degree, NDS, Szeged, Mostar, PI, Hosoya | no | no | no | no |
| PaDEL, alvaDesc | General descriptor calculators | SMILES/SDF | partial | no | no | no | no |
| DOPtools (Byadi et al., *Digital Discovery* 2025) | Descriptor calculation + hyperparameter optimisation | SMILES | ColorAtom/fragment descriptors | yes | no | no | partial |
| QSPRpred (*J Cheminform* 2024, 16:128) | Reproducible, serialised QSPR modelling | SMILES | via descriptor plug-ins | yes | no | CV, no size baseline | yes (model serialisation) |
| "QSAR in the Browser" (*JCIM* 2026, doi:10.1021/acs.jcim.6c01010) | No-code web QSAR with ChEMBL access | ChEMBL / upload | general descriptors | yes | no (target-based) | not stated | not stated |
| PyTIE (Jeyaraj et al., arXiv 2510.26813, 2025), MathChem (MATCH 71:657-680, 2014), topoindex (PyPI) † | Topological-index calculators | various | degree / NDS / entropy | no | no | no | no |
| Genre pipelines on GitHub †: aliberkanbektas/drug-indices-project_1 (PubChem → RDKit → indices); HakanKucuk38/Topological-Indices-PhysChem (indices + RF/XGB + SHAP, 4,711 antibiotics) | Scripts for one study | names / SMILES | yes | yes (latter) | no | no | no |

**Conclusion.** Computing many topological indices, and combining indices with ML, are both
already available. TopoQSPR must not claim novelty for either. "135 indices" is removed from the
selling points.

## 2. What no tool found does

1. Turn a disease concept into a **controlled, versioned drug manifest** (MeSH ID, phase rule,
   parent/salt normalisation, exclusion reasons, ChEMBL release).
2. Evaluate topological-index models with **leakage-free repeated nested CV against size and
   composition baselines**, as a default, not an option.
3. Quantify **evaluation optimism** (training → ordinary CV → nested → scaffold) and
   **re-analyse published disease-specific TI-QSPR studies** under that protocol.
4. Record **run provenance** suitable for archiving (commit, release, manifest, seeds, environment).

ChemGraphX is complementary: an index engine for arbitrary graphs and materials. TopoQSPR is the
study and evaluation layer for drug QSPR, and it uses ChemGraphX as an independent cross-check.

## 3. The genre this paper addresses

- Europe PMC hits for ("topological indices" OR "topological index") AND QSPR AND (drug OR drugs):
  16 (2022), 40 (2023), 62 (2024), 113 (2025), 78 (2026 to 2026-10-03); papers also mentioning machine
  learning, random forest or XGBoost: 5, 13, 18, 59, 47. (Recomputed on 2026-10-03 with `tools/genre_counts.py`; query:
  `("topological indices" OR "topological index") AND QSPR AND (drug OR drugs) AND PUB_YEAR:<year>`;
  the landscape audit's counts agree.)
- Typical studies use 9–50 manually listed drugs and in-sample or LOOCV R², often on properties that
  are computed (PubChem) or predicted (ChemSpider/ACD Labs) rather than measured. Verified
  examples: Almulla et al., *Sci Rep* 16:25695 (2026), 20 antifungal drugs, polynomial regression;
  Kara et al., *Sci Rep* 2025 (doi:10.1038/s41598-025-32906-x), n = 30, no CV, linear R² 0.999 (full text, PMC12830717, read by the audit skeptic) for
  heavy-atom count; Ramachandran et al., *Chem Pap* 2026, anti-glaucoma drugs with ML.
- The Arockiaraj group has several drug-QSPR papers (lung cancer *Sci Rep* 2025; TB *Heliyon* 2024;
  cardiovascular and anticancer *EPJE* 2024/2025; COVID *IJQC* 2021) †. Any re-analysis of these
  is chosen with the group and framed constructively.
- Methodological critiques to cite: Dearden, Cronin & Kaiser 2009 (21 QSAR error types);
  Gutman & Tošović 2013 (*J Serb Chem Soc* 78:805, doi:10.2298/JSC121002134G) and Furtula et al. 2013 on index quality (tested on octane isomers, where size is fixed); Ash et al. 2025 (*JCIM* 65:9398)
  for method comparison. No meta-analysis specific to drug-TI-QSPR was found.

## 4. What the evidence so far says (preliminary, not paper results)

Quick runs by two auditors, **before** any frozen protocol. They only motivate the design:

- All 135 v0.x indices are blind to element type and bond order.
- Whether indices beat size **depends on the learner**: on measured solubility (ESOL, n = 1,128),
  ridge gives size 0.34 vs TI 0.20, but random forest gives size 0.40 vs TI 0.53; element
  composition alone gives 0.72–0.85 and RDKit 2D descriptors 0.89–0.90. On ~800 approved drugs,
  TI beat size by +0.14 (logP) and +0.17 (TPSA).
- Consequence: results are reported as a **baseline ladder per learner**, never as one number.

## 5. 2026 opportunities (ranked; only 1–2 are in scope for v1.0)

1. A benchmark over every ChEMBL MeSH indication with enough approved small-molecule parents,
   reporting ΔR² over size per family (headline figure candidate).
2. Export the run record as a Workflow Run RO-Crate (PLOS ONE 2024, doi:10.1371/journal.pone.0309210)
   with a Zenodo DOI.
3. One learned baseline (Chemprop/CheMeleon) under the same protocol. **Deferred** per
   RECOMMENDATIONS.md unless the core result calls for it.
4. A thin MCP server (resolve_disease, build_manifest, run_qspr) for LLM agents. **Deferred**; at most
   one paragraph in the paper.

## 6. Venue

- **SoftwareX**: best fit; survived the red team's refutation attempt. Word limit reported as 4,000
  (excluding title, authors, references and metadata tables) by secondary sources; the primary guide
  returned HTTP 403, so confirm manually.
- **JOSS**: not eligible before about April 2027 (needs six months of public history) and web apps are
  generally out of scope.
- **J Cheminform**: has just published ChemGraphX from the same group; only with a strong
  methodological result.
- **Digital Discovery**: possible if the benchmark/replication result is strong.
