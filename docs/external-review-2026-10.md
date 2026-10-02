# External review (October 2026): claims checked

An external assessment of TopoQSPR rated it about 4.5/10 submission-ready, and about 8/10 for
a SoftwareX paper after the planned v1.0 work. Every factual claim in it was checked against
the code or the primary source on 2026-10-02. Verdicts: **holds**, **holds with caveat**,
**not checked**.

## Code claims

| Claim | Verdict | Evidence |
|---|---|---|
| CV is not truly nested (GridSearchCV fitted on all data, then outer CV) | holds | `ml.py` `run_ml_qspr`: `gs.fit(X_prop, y)` before `cross_val_score` |
| Selecting the best of 10 models by the same CV scores adds selection bias | holds | Best model chosen by `mean_r2 > best_r2` over the same folds |
| Hosoya Z wrong for >18 bonds | holds | `_hosoya_z`: enumeration only runs if `m <= 18`; 19–25 bonds return 1; >25 an ad-hoc formula |
| Reverse degree differs between code and paper | holds | Code `ml.py:225` uses n+1−d(v); draft used Δ−d(v)+1 (the literature form). The paper is now marked to reconcile this. |
| `DS1` and `DS2` are identical | holds | `ml.py:260-261` both add `s*s`; **missed by our own audit** |
| Failed indices are silently set to 0 | holds | `compute_topological_indices`: `df[key].fillna(0.0)` |
| Cohort is not reproducible (substring match, set truncation, no phase filter, keyword fallback) | holds | `fetcher.py` `fetch_chembl_drugs` |
| Caco-2 value is an invented formula | holds | `fetcher.py` `_compute_admet_rdkit`: `0.5 - 0.01*tpsa + 0.02*logp` |
| 135 × ~20 correlation tests use raw p-values | holds | `run_correlation` assigns stars from raw Pearson p |

## Literature and landscape claims

| Claim | Verdict | Evidence |
|---|---|---|
| ChemGraphX (J. Cheminform. 2026) computes many topological indices; M. Arockiaraj is a co-author; not cited in our draft | holds | Jacob, Clement, Arockiaraj, Jindal, *J Cheminform* 18:39 (2026-02-11), doi:10.1186/s13321-025-01141-x. Full text: input .pdb/.mol/adjacency; no SMILES, ChEMBL, disease cohorts, regression or ML. It is complementary, and is now cited. |
| "about 150 descriptors" in ChemGraphX | not checked | The count was not stated in the full-text passages searched |
| 2026 Sci. Rep. study: 9 neighbourhood-degree indices on 20 antifungal drugs | holds | Almulla et al., *Sci Rep* 16:25695 (2026), doi:10.1038/s41598-026-61020-9; uses linear/quadratic/cubic regression |
| 2026 anti-glaucoma QSPR with ML | holds | Ramachandran, Kavitha, Naresh Kumar, *Chemical Papers* (accepted 2026-02-22) |
| 2026 cancer paper with 23 drugs | not checked | No link given |
| ChEMBL stresses MeSH/EFO mapping, drug forms and phase for indication data | holds | Hunter et al., *J Med Chem* 68(19):19800–19827 (2025), doi:10.1021/acs.jmedchem.5c00920; now cited |
| DOPtools automates descriptor calculation and QSPR model optimisation | holds | Byadi et al., *Digital Discovery* 4(5):1188–1198 (2025); now cited |

## Recommendations

| Recommendation | Verdict | Note / action |
|---|---|---|
| Frame as a reproducible workflow, not "N indices + ML"; compete with ChemGraphX on rigour, not index count | **adopt** | Positioning paragraph added to the paper |
| True nested CV with model-family selection inside the outer fold | **adopt** | Roadmap M4; paper protocol updated |
| Benjamini–Hochberg FDR; effect sizes and CIs over stars | **adopt** | Roadmap M4; paper protocol updated |
| Failures recorded as missing, never zero | **adopt** | Roadmap M3 |
| M0 null / M_size / M_TI / M_standard (RDKit 2D), headline ΔR²(TI − size) with CI | **adopt, with caveat** | CV folds are not independent, so naive CIs are too narrow; use repeated nested CV with a corrected resampled t-test or a paired bootstrap of the whole procedure |
| 8–12 diseases with large enough approved cohorts | **adopt, with caveat** | Set the minimum cohort size before running (Alzheimer has only 13 approved drugs in ChEMBL 37) |
| Replicate 5–10 published disease-specific QSPR studies | **adopt, reduced to 3–5** | Needs the original drug lists **and** property values (often ACD/Labs-predicted, from ChemSpider); with n ≈ 20 every estimate is noisy, so report CIs; choose studies together with all co-authors and frame constructively |
| Title "…reveals the contribution of graph-topological descriptors beyond molecular size" | **defer** | Fix the title only after results; the finding may be positive or negative |
| SoftwareX first; J. Cheminform. only with a strong methodological result; Digital Discovery possible; Software Impacts as fallback | **adopt** | Consistent with our plan |

## Additions not in the external review

- Cross-check TopoQSPR index values against **ChemGraphX** where the two overlap. This is an independent implementation from the senior co-author's group, and the strongest available validation for M3.
- Cite the source of the hybrid geometric-harmonic-Zagreb indices (Arockiaraj et al., *SAR QSAR Environ Res* 34:569–589, 2023), which TopoQSPR implements.
- Cite the reverse-degree and R_k definitions from their primary sources, or drop variants that lack one.
