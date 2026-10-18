# Doc 14 r2 full run (four-arm fixed-learner protocol)

- Protocol frozen at `1457bc9358fad88e449a13a218d60d1f0f985077` after the owner's approval on 2026-10-18. The run made no deviations from it.
- Run on the workstation, 2026-10-18, from 13:40 to 15:51 EDT.
  - Clean worktree at that commit.
  - Python 3.12.15 with the pinned set; `n_jobs=8`, nice 10.
  - 760/760 units verified; 46,666 CPU-seconds of training (98.6 min wall).
- Checks:
  - audit 32/0;
  - build 25/0;
  - robustness 3/0.
- `run_full.sh` is the exact driver. Per-unit outputs, feature matrices and the restored OPERA
  labels stay local and are not committed (licence).
- Label: exploratory, post hoc (doc 14 header). Not prospective, not independent validation.

## Co-primary results (Holm over 12, δ = 0.02 R²; 95% CI)

| context | learner | C+T − C | verdict | C+E+T − C+E | verdict |
|---|---|---|---|---|---|
| mp_E3a | extra trees | +0.087 [0.068, 0.106] | > δ | +0.006 [−0.001, 0.013] | equivalent |
| mp_E3a | SVR | +0.087 [0.069, 0.104] | > δ | −0.007 [−0.014, −0.000] | inferior, not shown beyond δ |
| mp_E3b | extra trees | +0.098 [0.074, 0.122] | > δ | −0.005 [−0.020, 0.011] | equivalent |
| mp_E3b | SVR | +0.100 [0.078, 0.122] | > δ | +0.000 [−0.010, 0.010] | equivalent |
| external, unseen formula | extra trees | +0.122 [0.099, 0.145] | > δ | −0.035 [−0.051, −0.019] | inferior by more than δ |
| external, unseen formula | SVR | +0.221 [0.183, 0.260] | > δ | −0.003 [−0.026, 0.020] | inconclusive (pre-labelled underpowered) |

- **Secondary:** C+E − C exceeds δ in every co-primary cell (q < 0.001).
- **Linear learners** (descriptive, not co-primary): C+E+T − C+E is +0.03 to +0.07 internally and up to +0.14 externally.
  - Reading: the indices supply nonlinear functions of the edge counts that linear models cannot build themselves.
- **External extra trees, training seeds 0–4:** C+E+T − C+E = −0.036 (SD 0.003); C+T − C = +0.123 (SD 0.009).
- **Robustness:**
  - Corrected-t coverage at nominal 95%: ridge 0.998, SVR 0.990. The naive t-test gives 0.858 and 0.850.
  - The seed-0 external refit reproduces the archived predictions.
