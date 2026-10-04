# TopoQSPR vs ChemGraphX: scorecard

Goal (Singati2, 2026-10-04): the TopoQSPR paper and website must be demonstrably better than
ChemGraphX (Jacob, Clement, Arockiaraj, Jindal, *J. Cheminform.* 18:39, 2026; site
https://chemgraphx.vercel.app/, code github.com/aditya-jindal/chemgraphx @ d8fea04). ChemGraphX is
the senior co-author's group's tool; the comparison is factual and neutral, and the two tools are
complementary (ChemGraphX: index calculator for arbitrary graphs; TopoQSPR: QSPR study workflow).
Every "ChemGraphX" cell below was checked in its source code or against its live API (2026-10-04).

Legend: ✅ has it · ❌ does not · ⏳ TopoQSPR work in progress (target).

## Website

| Capability | ChemGraphX | TopoQSPR (target) |
|---|---|---|
| Input: .mol file | ✅ | ⏳ |
| Input: .pdb file | ✅ | ⏳ |
| Input: adjacency list / matrix (general graphs) | ✅ | ⏳ |
| Input: SMILES, drug name, SDF, batch CSV | ❌ | ⏳ |
| Input: disease → curated drug set (ChEMBL, manifest preview) | ❌ | ⏳ (Phase 3) |
| Degree-, degree-sum-, distance-based indices | ✅ ~97 | ✅ 102 verified (reverse-degree R_k, hybrids, exact Hosoya Z not in ChemGraphX) |
| Entropy measures per index | ✅ 52 degree + 21 distance | ⏳ |
| Formula and primary reference shown for every index | ❌ | ⏳ |
| 2D structure depiction / edge-partition view | partial (graph plot) | ⏳ |
| Download results (CSV/JSON) | ❌ | ⏳ |
| QSPR modelling with baselines and confidence intervals | ❌ | ⏳ (Phase 3, fast web mode) |
| Run record / provenance download | ❌ | ⏳ |
| Undefined values reported with a reason | partial ("indeterminate", inf) | ✅ (NaN + reason) |
| Output rendering without innerHTML injection risk | ❌ (chart via innerHTML) | ⏳ |
| Persistent jobs, rate limits, health endpoint | n/a | ⏳ |
| Automated test suite / CI | ❌ (none found) | ✅ 786 offline tests, CI |
| Accessibility (labels, keyboard, contrast) | to assess | ⏳ |
| Cold / warm response for one molecule | 3.4 s / 0.6 s (live, 2026-10-04) | ⏳ measure |

## Paper

| Element | ChemGraphX paper | TopoQSPR paper (target) |
|---|---|---|
| Scientific contribution | new web index calculator | reproducible disease-anchored QSPR workflow + methodological study |
| Index validation | comparison with existing tools on example structures | five independent routes incl. closed forms, brute force, Mordred, ChemGraphX 41/41, printed values; errata found in primary sources; independent-literature check |
| Statistical evaluation | none (descriptor tool) | repeated nested CV, baseline ladder, corrected-t CIs, BH/BY, Y-randomisation, bootstrap, scaffold CV |
| Experiments on real data | illustrative structures | cohort sensitivity (A), size/composition confounding (B), evaluation optimism (C), replication of four published studies (D) |
| Reproducibility | code on GitHub | tagged release, Zenodo DOI, run records, lock file, Docker, tests, preregistered plan |
| Related-work comparison table | ✅ | ✅ (more tools, each cell verified) |
| Figures / tables | 10 / 3 | target ≥ 6 informative figures from result files |

The goal is met when every ⏳ is ✅ and verified, and the paper's results come from the frozen,
tagged analysis.
