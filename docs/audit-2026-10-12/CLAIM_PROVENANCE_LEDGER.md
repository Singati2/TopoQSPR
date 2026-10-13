# Claim provenance ledger (audit 2026-10-12)

Manuscript claim → data → run (old → new commit) → test → effect → source → **permitted wording**.
Commit mapping: `sha_map_2026-10-12.tsv`. "Frozen" = `main_confirmatory` (102 TI, random,
grouped, commit 659fac3 → 61f63d8, clean tree).

| # | Claim (location) | Data / run | Test, effect | Status | Permitted wording |
|---|---|---|---|---|---|
| 1 | "On 1821 approved drugs with measured end points" (abstract, `main.tex:36`) | frozen `targets.csv` | counts: 1821 parents, 959 with ≥ 1 endpoint, 711/505/369/351 | **FAIL** (F3) | "1821 approved drugs, of which 351–711 per endpoint have a measured value" |
| 2 | TI beat size for 1 of 4 endpoints after BH; log P ΔR² 0.15 [0.08, 0.23], q = 0.002 (abstract, `:202`) | frozen `confirmatory.csv` | corrected t, BH m = 4; PR #2 test recomputes BH | **PASS** | as written, plus "random grouped folds, 102 indices" |
| 3 | TI "never reach the element-composition baseline" (abstract) | frozen `summary.csv` | fold-mean R² composition > TI for 4/4; E1 both splits | **PASS** (descriptive; E1, E2) | as written |
| 4 | Paired TI+composition vs composition "not run" (`:202`) | — | E1 (exploratory, asinh): no 95% CI excludes 0 in 4 endpoints × 2 splits; E2 pending | **LIMITED**: true when written; superseded by E1/E2 | "an exploratory paired test (dated amendment) found no detectable gain …; not an equivalence claim" → final: "in a dated exploratory amendment (similarity-cluster split, BH over four endpoints) no endpoint showed a detectable gain beyond composition (q ≥ 0.45); equivalence within ±0.05 R² was not shown" |
| 5 | Pooled R² < −1 cells explained by Hosoya Z (`:202`) | frozen `predictions.csv`, `features.csv` | trace: TI cells → Z; **rdkit2d cell → RDKit Ipc** | **FAIL** in part (F1) | "…Hosoya's Z for the index sets and RDKit's Ipc for the descriptor set…" |
| 6 | Scaffold: log P gain "the same size but wider" (`:235`) | `main_scaffold` (104) vs frozen (102); E1 scaffold with 102 | +0.153 → +0.100 [−0.050, +0.249] with the registry held fixed | **FAIL** (F5) | "about a third smaller on scaffold folds, with an interval that includes zero" |
| 7 | Amended SVR control: bootstrap CIs agree with corrected t (`:208`) | `main_confirmatory_protocol_svr` (104, fe92837 → 1a504a9, **dirty**) | archived values match | **LIMITED** (F12, F16) | add "family-specific control (SVR), not the selected procedure; run from an uncommitted tree whose package source was unchanged" |
| 8 | "leakage-free, preregistered statistics" (`:38`, `:265`; supplement `:63`) | `EXPERIMENT_PLAN.md`, tag `plan-freeze-v1` | no OSF record | **FAIL** (F4) | "statistics pre-specified in a plan frozen as a git tag before the main runs (not externally registered)" |
| 9 | Free-text retrieval gives 4–16× more structures (`:223`, `:265`) | `results/expA` | macros `expAminRatio` 4, `expAmaxRatio` 16 | **LIMITED**: ratios 3.59–15.97 re-derived; "4" rounds 3.59 up (F18) | "3.6–16 times" |
| 10 | v0.x protocol overstates training R² by up to 0.60 (`:265`) | `results/expC` | macro `expCmaxBtrain` 0.60 | **PASS** (re-derived: 0.909 − 0.313 = 0.595) | as written |
| 11 | Three published studies "statistically indistinguishable from the size baseline" (abstract) | `results/expD/*/nested.csv` | CI widths: Godlin 0.24–0.27, Kabeer 0.7–0.9, Almulla 1.9–2.1; 0/78 significant | **FAIL** as worded (F6) | "none shows a significant gain over size; only the Godlin set is precise enough to bound it (|ΔR²| ≤ 0.2 in 13 of 27 fits)" |
| 12 | Index values agree with ChemGraphX on all 41 shared indices (`:173`) | `docs/CHEMGRAPHX_CROSSCHECK.md`, tests | test suite passes (`test_chemgraphx_crosscheck.py`) | **PASS** | as written, **plus** disclosure that M. Arockiaraj co-authored ChemGraphX |
| 13 | Web interface workflow and homepage figure (`:156`, `:168`) | live site | serves v0.x | **UNRESOLVED** (F7) | keep the `\todo`; do not submit until the v1.0 site is deployed |
| 14 | Run reproduced from its run record (`:86`, `:200` `\runCommit`) | run records | cited commit 659fac3 absent from GitHub | **FAIL** until the map or archive is cited (F0) | cite `61f63d8` and the SHA map, or a Zenodo archive |
