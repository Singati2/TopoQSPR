# ResearchOS quickstart

Requirements: Python ≥ 3.10 and Git. ResearchOS uses only the standard library. pytest is needed only for the tests.

```bash
git clone https://github.com/Singati2/TopoQSPR.git && cd TopoQSPR
git checkout feature/researchos-orchestrator

# isolated environment (nothing to install for the CLI itself)
python3 -m venv /tmp/ros-venv && /tmp/ros-venv/bin/pip install -q pytest
export PY=/tmp/ros-venv/bin/python

# optional: keep state outside the repo while trying it out
export RESEARCHOS_HOME=/tmp/ros-state

$PY -m researchos doctor            # read-only environment check (exit 0 = required checks pass)
$PY -m researchos status            # branch, HEAD, claims, tasks, protected processes
$PY -m researchos agents list       # which roles exist (3) and which are interfaces only

# the worked example on the real E3 artifacts (no training)
$PY -m researchos claims import --path researchos/examples/e3_claims.json
$PY -m researchos task propose --from-file researchos/examples/melting_point_gap.md   # prints T-…
$PY -m researchos task dry-run T-…  # what would run, what is blocked and why
$PY -m researchos task run T-…      # runs read-only steps; stops at approvals; state -> blocked
$PY -m researchos claims show C-EDGE10
$PY -m researchos report T-… --output /tmp/decision.md
$PY -m researchos evidence verify docs/research/out/E3/external/SHA256SUMS
$PY -m researchos audit --limit 20

# tests
$PY -m pytest tests/researchos -q
```

## Approving a step (owner only)

```bash
# interactive: you are asked to type the task id
$PY -m researchos task approve T-… --scope analysis-only
# non-interactive: repeat the task id
$PY -m researchos task approve T-… --scope analysis-only --confirm T-…
```

The example's `fold_contrast` step also reads a run archive outside the repository. It runs only
after the owner both approves `analysis-only` **and** adds that directory to `allow_roots` in a
`researchos.json` at the repo root:

```json
{"allow_roots": [".", "~/topoqspr_runs/research_2026-10-13_E3/cells"]}
```

`compute` approvals are recorded, but v0.1 never launches compute.

## Exit codes

`0` ok · `1` a check failed, a step failed or an action was refused · `2` usage, configuration or
unknown repository state.
