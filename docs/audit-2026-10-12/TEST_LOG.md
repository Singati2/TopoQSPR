# Test and execution log, audit 2026-10-12

Environment: Linux 6.8, 16 cores, venv `~/.venvs/topoqspr` (Python 3.10.12, RDKit 2026.03.6,
scikit-learn 1.7.2, NumPy 2.2.6, pandas 2.3.3). Threads pinned (`OMP/OPENBLAS/MKL_NUM_THREADS=1`,
`PYTHONHASHSEED=0`) for the experiment runs.

| # | Command | Commit | Expected | Actual |
|---|---|---|---|---|
| 1 | `python -m pytest -q` | `fecb800` (HEAD, main) | all pass | 1103 passed, 3 skipped, 8 deselected, 92 s |
| 2 | `python -m pytest -q tests/test_stats.py` | audit branch, after F2 fix | new test passes; old p=0 assertion replaced | 23 passed |
| 3 | `python -m pytest -q` | audit branch, after F2 fix | all pass | 1104 passed, 3 skipped, 8 deselected, 120 s |
| 4 | `python -I docs/audit-2026-10-12/scripts/catastrophic_predictions.py ~/topoqspr_runs/main_confirmatory docs/audit-2026-10-12/out` | audit branch | trace R² < −1 folds | 8 folds, 24 rows → `out/catastrophic_main_confirmatory.csv` (F1) |
| 5 | `python -I docs/audit-2026-10-12/scripts/representation_collisions.py ~/topoqspr_runs/main_confirmatory docs/audit-2026-10-12/out` | audit branch | floors and collision classes | `out/collision_floor.csv`, `out/collision_classes.csv`, `out/ethanol_dimethyl_ether.csv` (104/104 indices equal) |
| 6 | `python -I docs/audit-2026-10-12/scripts/sensitivity_composition_asinh.py ~/topoqspr_runs/main_confirmatory ~/topoqspr_runs/audit_2026-10-12_sensitivity` | audit branch | archived raw sets reproduce exactly; E1 results | first attempt stopped at its reproduction assert (default CSV float parsing, F15); rerun with round-trip parsing 17:07–18:27: 8 target×split cells, all reproduction asserts passed; results in RESULT_RECONCILIATION.md §3 |
| 7 | `pdflatex`×3 + `bibtex` of `paper/main.tex`, `paper/supplement.tex` (out of tree) | `fecb800` | compiles | 34 pp / 14 pp; 0 errors; 0 undefined refs/cites; 1 BibTeX warning (`lundberg2017` empty pages) |
| 8 | `curl https://moleculens.onrender.com/{,health,openapi.json}` | — | v1.0 app | v0.x app served; `/health` reports hard-coded "1.0.0" (F7) |
| 9 | PR #2 tests (`tests/test_adversarial_audit_20261012.py` from `origin/audit-tests-20261012` @ `5fc799e`), run as a temporary copy | audit branch / `fecb800` | — | 5 passed, 1 failed (`test_a4`, asserts old p = 0; F2) / 6 passed |
| 10 | same, with `pr2_test_a4.patch` applied | audit branch | all pass | 6 passed |
| 11 | `python tools/audit_permutation_exchangeability.py` (PR #2) | `5fc799e` tool, installed package | anti-conservative | 96/200 rejections, rate 0.480 (nominal 0.05) |
| 12 | `python -I docs/audit-2026-10-12/scripts/permutation_block_check.py ~/topoqspr_runs/main_confirmatory ~/topoqspr_runs/main_confirmatory_protocol_svr` | audit branch | 0 row-wise rows | 0 row-wise rows in all 8 target×run cells; exit 0 (F17) |
| 13 | `python -I docs/audit-2026-10-12/scripts/e2_cluster_ecfp.py --smoke …` (permuted y, logD, cluster, seed 0) | `533dcbb` | completes | completed in 251 s; numbers discarded unread |
| 14 | `python -I docs/audit-2026-10-12/scripts/e2_cluster_ecfp.py ~/topoqspr_runs/main_confirmatory ~/topoqspr_runs/audit_2026-10-12_E2` | `533dcbb`, clean | 12 cells; E1 reproduction asserts | 12 cells, 0 failed folds, all asserts passed, 14 371 s |
| 15 | independent recomputation of E2 primary CIs, BH q and group integrity (inline script, see RESULT_RECONCILIATION §7) | — | match | max abs diff 1.1e-16; q identical; 0 of 23 945 group×repeat combos straddle folds |
| 16 | `python -I docs/audit-2026-10-12/scripts/permutation_replay_and_calibration.py ~/topoqspr_runs/main_confirmatory ~/topoqspr_runs/main_confirmatory_protocol_svr` | `8b9d95b` + script | 0 split groups | 0 / 6400 split; calibration table in `out/permutation_calibration.csv` (36 min) |
| 17 | `python -m pytest -q tests/test_protocol.py tests/test_pipeline.py tests/test_stats.py` | `90ca304` (pre-commit tree) | pass | 129 + new tests pass |
| 18 | `python -m pytest -q` | `90ca304` (pre-commit tree) | all pass | **1111 passed, 3 skipped, 8 deselected**, 221 s |
| 19 | PR #2 test file + `pr2_test_a4.patch`, against `90ca304` | — | 6 pass | 6 passed |
| 20 | Crossref lookups of 5 DOIs; PyPI DOPtools metadata | — | match report | all match (F18 section) |

Not run (blocked or out of scope): Lean formalisation (no toolchain); isomer-controlled study (SCIENTIFIC_RESEARCH_PLAN.md E3); clean rerun of the amended SVR control (~6 h); end-to-end raw-input → table reconstruction (`topoqspr rerun`, ~5 h);
re-download of source datasets (cached copies with pinned SHA-256 were used); manual check of
epirubicin's melting point against a primary source (F10).
