# TopoQSPR vs ChemGraphX: scorecard

Goal (Singati2, 2026-10-04): the TopoQSPR paper and website must be demonstrably better than
ChemGraphX (Jacob, Clement, Arockiaraj, Jindal, *J. Cheminform.* 18:39, 2026; site
https://chemgraphx.vercel.app/, code github.com/aditya-jindal/chemgraphx @ d8fea04). ChemGraphX is
the senior co-author's group's tool; the comparison is factual and neutral, and the two tools are
complementary (ChemGraphX: index calculator for arbitrary graphs; TopoQSPR: QSPR study workflow).
Every "ChemGraphX" cell below was checked in its source code or against its live API (2026-10-04).

Legend: ✅ has it · ❌ does not · ⏳ TopoQSPR work in progress (target).

**Similarity requirement (Advik, 2026-10-05):** the website should look and work like ChemGraphX and
the paper should follow the ChemGraphX article's structure. Website: done on `main` (three-tab shell,
hero, input-format card, cite card, serif/light/blue; side-by-side in `docs/WEB.md`). Paper: section
structure mapped in `docs/PAPER_FORMAT.md` and applied to `paper/main.tex`.

**Status 2026-10-06:** every experiment is run from the frozen plan and written into the paper from
generated files: A (`results/expA`), B (confirmatory, scaffold ladder, index-family ladder, frozen
ridge controls and the amended SVR controls, `results/expB`), C (`results/expC`), D (`results/expD`).
Supplement (S1 index table from the registry, S2 OECD mapping, S3 checklist) compiles. Independent
literature value check finished (199 cases, 7 papers; `docs/LITERATURE_VALUE_CHECK.md`). Full offline
suite green (1 101 tests). Open: author items in the paper (URLs, DOIs, CRediT, funding, AI wording),
website adversarial review (agent limit until 2026-10-07) and Render deployment (needs an account),
purge of `related_papers/` before the repository goes public.

## Website

| Capability | ChemGraphX | TopoQSPR (target) |
|---|---|---|
| Input: .mol file | ✅ | ✅ (`inputs.py`; also SDF multi-record) |
| Input: .pdb file | ✅ | ✅ |
| Input: adjacency list / matrix (general graphs) | ✅ | ✅ (validated: symmetric, 0/1, no loops, connected) |
| Input: SMILES, SDF, batch lists | ❌ | ✅ (batch of 50 molecules in 7 s) |
| Input: disease → curated drug set (ChEMBL, manifest preview) | ❌ | ✅ (`/api/headings`, `/api/cohort`, manifest review and CSV before any computation) |
| Degree-, degree-sum-, distance-based indices | ✅ ~97 | ✅ 104 verified (reverse-degree R_k, hybrids, exact Hosoya Z not in ChemGraphX) |
| Entropy measures per index | ✅ 52 degree + 21 distance | ✅ 91 bond-wise entropies; 1,414 values agree with ChemGraphX's code (max rel. diff 8e-16); distance entropies not implemented |
| Formula and primary reference shown for every index | ❌ | ✅ (`/api/catalog`, 191 entries, searchable) |
| 2D structure depiction / edge-partition view | partial (graph plot) | ✅ (RDKit SVG depiction; degree edge partition per structure) |
| Download results (CSV/JSON) | ❌ | ✅ (zip with indices, reasons, catalogue; on-page table) |
| QSPR modelling with baselines and confidence intervals | ❌ | ✅ (exploratory web mode: ladder null/size/composition/TI/RDKit 2D, ΔR² and ΔMAE with 95 % CI; confirmatory protocol via CLI) |
| Run record / provenance download | ❌ | ✅ (`run_record.json` and results zip per job) |
| Undefined values reported with a reason | partial ("indeterminate", inf) | ✅ (NaN + reason) |
| Output rendering without innerHTML injection risk | ❌ (chart via innerHTML) | ✅ (no innerHTML; textContent only) |
| Persistent jobs, rate limits, health endpoint | n/a | ✅ (SQLite jobs, per-IP limits, caps, /health) |
| Automated test suite / CI | ❌ (none found) | ✅ ~880 offline tests incl. 35 web tests, CI |
| Accessibility (labels, keyboard, contrast) | not assessed | partial ✅ (every input labelled, live-region status, lang/viewport; tested); contrast audit pending |
| Response for one molecule | 3.4 s cold / 0.6 s warm, one of three views (live AWS Lambda, 2026-10-04) | 0.73–0.82 s for all 191 values (local, machine under full load, 2026-10-05); **not comparable until TopoQSPR is deployed**: Render free tier sleeps after 15 min idle, so its cold start will be slower than Lambda's |

## Paper

| Element | ChemGraphX paper | TopoQSPR paper (target) |
|---|---|---|
| Scientific contribution | new web index calculator | reproducible disease-anchored QSPR workflow + methodological study |
| Index validation | comparison with existing tools on example structures | five independent routes incl. closed forms, brute force, Mordred, ChemGraphX 41/41, printed values; errata found in primary sources; independent-literature check |
| Statistical evaluation | none (descriptor tool) | repeated nested CV, baseline ladder, corrected-t CIs, BH/BY, Y-randomisation, bootstrap, scaffold CV |
| Experiments on real data | illustrative structures | **A, B, D done** (A: 4 diseases, naive query returns 4–16× more structures, Jaccard 0.06–0.28; B: pooled 1,821 approved drugs, 4 measured endpoints, frozen plan `plan-freeze-v1`; D: 3 published studies replicated on their own data); C running |
| Reproducibility | code on GitHub | tagged release, Zenodo DOI, run records, lock file, Docker, tests, preregistered plan |
| Related-work comparison table | ✅ | ✅ (more tools, each cell verified) |
| Figures / tables | 10 / 3 | generated only from result files (`tools/make_results.py`, `make_expA.py`): so far 4 figures + 7 tables (workflow, homepage, inputs, A, B forest, D) ; C adds one more of each |

The goal is met when every ⏳ is ✅ and verified, and the paper's results come from the frozen,
tagged analysis.
