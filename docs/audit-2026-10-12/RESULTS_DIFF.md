# Results diff (audit 2026-10-12)

## Archived and published values: **none changed**

`git diff --stat fecb800 HEAD -- results/ paper/ studies/` is empty on the audit branch. Every
archived CSV, run record, study config, generated macro and the manuscript are byte-identical
to `fecb800`. The run directories under `~/topoqspr_runs/main_*` and `exploratory_families`
were only read (inputs hashed in EXPERIMENT_AMENDMENT.md §1).

The only code change is `src/topoqspr/stats.py` (F2). It affects `corrected_ttest` only when the
paired fold differences have zero variance and a non-zero mean. No archived paired group meets
that condition (all zero-variance groups are reference-vs-itself, mean 0, p = 1, unchanged).
Re-running the generators would therefore reproduce every published number. This is argued
from the condition, not by re-running `tools/make_*.py`, which was not done.

## New values (exploratory; separate output directories)

| analysis | directory | status |
|---|---|---|
| E1 sensitivity (asinh; composition-conditional; identical-registry scaffold) | `~/topoqspr_runs/audit_2026-10-12_sensitivity`, copies `out/sensitivity_*.csv` | done; archived raw sets reproduced to 1e-9 |
| E2 amendment run (cluster split, ECFP) | `~/topoqspr_runs/audit_2026-10-12_E2`, copies `out/E2_*.csv` | done; E1 reproduced to 1e-9; results RESULT_RECONCILIATION §7 |
| catastrophe trace, collision floor, permutation block check | `out/` | done |

Hashes: E1 `summary.csv` 4a21b4ac…, `paired_tests.csv` f6b04197…, `folds.csv` 4b35aa73….
E2 `summary.csv` 5f42656f…, `paired_tests.csv` 2a1dc1b4…, `folds.csv` 9d1ada8f….

## Code changes after the audit started and their effect on archived outputs

| commit | change | effect on archived/published values |
|---|---|---|
| `33e8198` | corrected t: zero-variance, non-zero-mean → NaN p/CI | none (branch never reached; F2) |
| `90ca304` | permutation p withheld when the row-wise fallback fires; `rowwise_fallback_rows` recorded; dirty runs archive the diff | none: fallback rows = 0 in every archived target (replayed); diff only for dirty runs; config hashes unchanged; `rerun` byte-identity test passes |
