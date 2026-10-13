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

Not run (blocked or out of scope): Lean formalisation (no toolchain); similarity-cluster split,
ECFP baseline and isomer-controlled study (planned, SCIENTIFIC_RESEARCH_PLAN.md E2–E4);
re-download of source datasets (cached copies with pinned SHA-256 were used); manual check of
epirubicin's melting point against a primary source (F10).
