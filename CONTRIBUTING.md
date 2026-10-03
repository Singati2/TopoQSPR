# Contributing to TopoQSPR

Thank you for helping. TopoQSPR is research software behind a paper, so correctness and
provenance come first.

## Ground rules

- **Every index needs a primary source.** A new or changed topological index must cite the paper
  that defines it, state the formula in its docstring, and come with tests: a closed form on a
  small graph and, where possible, a value printed in the source. Record it in
  `docs/INDEX_VERIFICATION.md`.
- **No silent failures.** Undefined or failed values are `NaN` with a reason, never `0`.
- **Determinism.** Outputs must not depend on hash seeds, dict/set order or unseeded randomness.
- **No peeking.** Do not look at model performance on the main-panel diseases before the
  experiment plan (`docs/EXPERIMENT_PLAN.md`) is frozen; use the quarantined pilot disease or
  synthetic data for debugging.
- **Data licences.** ChEMBL-derived files are CC BY-SA 3.0. Never commit publisher PDFs or data
  whose licence forbids redistribution.

## Workflow

1. Branch from `v1.0-dev`.
2. `pip install -e ".[test]"`, then `pytest` (network tests are opt-in: `pytest -m network`).
3. Open a pull request describing what changed, why, and which tests cover it.
