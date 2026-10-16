# ResearchOS v0.1: implementation report (2026-10-16)

Branch `feature/researchos-orchestrator`, created as a separate worktree (`~/moleculens-researchos`)
from `origin/research-2026-10-13` at `8e5bc2e`. The research checkout (`~/moleculens`, branch
`research-2026-10-13`) was not touched: still clean, still at `8e5bc2e`, checked before and after.
Remote `main` was unchanged at `3ef0b50` throughout. Commit IDs for this work are listed in
`git log 8e5bc2e..feature/researchos-orchestrator`.

## 1. What was built (real) vs not

**Built and tested:**
- Standard-library Python package `researchos/` (about 1,840 lines; tests about 510) with:
  - typed claim/evidence/task schema;
  - task state machine;
  - SQLite store with migrations and append-only evidence and audit;
  - realpath allowlist with a private-path deny list;
  - redaction and untrusted-text screening;
  - read-only Git, process and health inspectors;
  - sha256 and Git-provenance manifest verification.
- Seven registered adapters:
  - five read-only;
  - one analysis-only (approval required);
  - one compute adapter, which is never launched;
  - plus an LLM stub.
- Approval gate bound to the spec hash, dry-run planner, locked idempotent executor with crash resume, Markdown decision report.
- CLI: `doctor`, `status`, `agents list`, `claims list|show|import`, `tasks list`, `task propose|dry-run|approve|run`, `report`, `evidence verify`, `audit`.
- 37 tests. The ADR and the six docs requested.

**Not built:** see `LIMITATIONS.md`. That covers the dashboard, SLURM/HPC, the math/stats/publication agents, the LLM adapter, compute execution and signed approvals.

## 2. Orchestrator due diligence

Details and sources are in `ADR-0001-build-vs-adopt.md`.
- **FirstMate:** `kunchenguid/firstmate`, MIT, verified by the GitHub API on 2026-10-16. It is an agent "distro" for *coding* agents (worktrees, tmux panes, PRs), not an evidence system.
- **"OrchestraOS":** ambiguous. The only matching project is `Tulum-DAO/orchestraos` (Apache-2.0, created 2026-09-17, self-described "not a finished product"). It is a server harness for coding-agent fleets with approval cards. Other repos with that name are unrelated or unlicensed.
- **Code reuse:** no code was copied, nothing was installed, no setup scripts were run.
- **Concepts re-implemented:** isolated worktrees, escalate-only-real-decisions, restart-proof on-disk state, approve/deny records, `doctor`, append-only state.

## 3. Safety preflight

- **Machine, 11:16 EDT:** load 0.5 on 16 CPUs; 46 GB RAM available; 68 GB disk free (93% used); 2× RTX 2080 Ti mostly idle; no crontab; no SLURM on the host.
- **Owner jobs:** no protected processes at the start.
- **E3:** finished at 07:14 (log line `end 2026-10-16T07:14:01`). This work did not re-run or alter it.
- **Owner process during the build:** at 11:30 the owner's `test_baselines.py` (PID 1104828, launched from another session's working directory) started. `researchos status` detected it. Nothing in this work signals, traces or attaches to processes. The detection reads `/proc/<pid>/cmdline`.
- **No destructive Git actions:** no force-push, no history rewrite, no merge, no branch deletion, no visibility change. The feature branch's upstream was unset at creation, so a bare `git push` cannot target `research-2026-10-13`.

## 4. Demo: the E3 question

**Question:** why does TIa_composition (R² 0.471) beat edgetypes_composition (0.423) on the source-disjoint unseen-formula external test, and what can be inferred without training?

**Run:** `python -m researchos task run T-a5cea8c5230c`. The record is at `demo/T-a5cea8c5230c_decision.md`.

**Inputs verified:**
- `docs/research/out/E3/external/summary.csv`: sha256 `a13dfa7e…`, commit `4a8cb1d`.
- `paired.csv`: sha256 `b4232c7d…`, commit `4a8cb1d`.
- Both match `SHA256SUMS`, and the committed blob equals the working file.
- `info.json` and `predictions.csv.gz` are listed in the manifest but deliberately kept outside Git. They are reported as unverifiable here, not as passes.

**Findings:**
1. All 14 archived external contrasts re-derive exactly from the per-set R². They also match the earlier independent recomputation from predictions (max |diff| 0).
2. **The gap is learner-confounded.** TIa_composition selected extra trees and edgetypes_composition selected SVR. Every internal cell also mixed families across folds. No interval exists for the TIa_composition − edgetypes_composition contrast anywhere.
3. Across the six contexts, edge types match or beat the full set in five:
   - MP E3a −0.006, E3b −0.001;
   - logP E3a −0.003, E3b −0.020;
   - external known-formula −0.026.
   The exception is external unseen-formula (+0.048).
4. **Claim C-EDGE10 is downgraded from `supported` to `inconclusive`.** C-EDGE10 is "the gain is carried by 10 numbers", as worded in `PUBLICATION_READINESS.md` at `8e5bc2e`. The audit log records the reason. The single evidence row that had supported it is `qualitative_llm_review`, which is not numerical.
5. The deterministic reviewer raised seven objections (three high). A narrowed wording is proposed for the owner.
6. Two steps are waiting for the owner and were not run:
   - **Fold-level interval** (analysis-only): needs approval *and* adding the run archive to `allow_roots`.
   - **Fixed-learner experiment** (compute): the protocol draft, the resource estimate (arm 1 ≈ 15–25 min, grounded in the 42-min E3 external stage) and the input checks are in the record. v0.1 will not launch it.

## 5. Tests: exactly what ran

| command | result |
|---|---|
| `~/.venvs/topoqspr/bin/python -m pytest tests/researchos -q -p no:cacheprovider` | **37 passed** (≈3.5 s), exit 0 |
| same, in a fresh clone with a new venv (`/usr/bin/python3 -m venv`, then `pip install pytest` only) | **37 passed**, exit 0 |
| mutation check: grant `analysis-only` by default | 3 tests failed as intended; reverted |
| mutation check: `abspath` instead of `realpath` (symlink escape) | 1 test failed as intended; reverted |
| `python -m compileall -q researchos tests/researchos` | ok |
| `python -W error -c "import researchos.cli, researchos.report"` | ok |
| `git diff --check` | ok |
| AST unused-import scan | one unused import found and removed |
| `pytest tests --collect-only` (whole repo) | 1151/1159 collected (8 deselected by the existing `network` marker); the TopoQSPR suite itself was **not re-run** |
| fresh-clone CLI walk-through (`doctor` 0, `status` 0, `agents list` 0, `claims import` 0, `task propose` → `T-a5cea8c5230c`, `dry-run` 0, `run` 0 → blocked, re-run of the blocked task 1 "refused", wrong `--confirm` 1 "refused", `evidence verify` ok, `report` 0) | as listed |

**Negative paths covered:**
- unapproved step;
- read-only gate when configured;
- approved compute still blocked;
- wrong or missing confirmation;
- spec change voids approval;
- injected "owner approved / rm -rf / force push / curl | sh" text;
- outside path, `..`, symlink escape, private path;
- hash mismatch, missing file, uncommitted change, conflicting manifest;
- concurrent run (lock, in a separate process);
- crash mid-run then resume (no re-execution, attempt + 1);
- completed task not re-run;
- unknown adapter, scope mismatch, cycle, duplicate task block;
- append-only violations;
- newer store schema;
- unknown config key;
- no `claude` on PATH;
- running outside Git;
- a live `test_baselines.py`-named process left untouched.

**Not run:** ruff and mypy (not installed in the project environment). No load or performance tests. No tests on macOS.

## 6. Acceptance criteria (prompt §6)

| # | criterion | status |
|---|---|---|
| 1 | fresh setup from QUICKSTART in an isolated environment | met (fresh clone + new venv, above) |
| 2 | `doctor`/`status` inspect the real project safely | met |
| 3 | an E3 artifact ingested with verified checksum and commit | met (external summary/paired, commit `4a8cb1d`) |
| 4 | reproducible, source-linked report without training | met (task id is deterministic; same id in the fresh clone) |
| 5 | unauthorized execution fails closed, under test | met |
| 6 | ledger distinguishes recomputation from LLM opinion | met (`independence`, `is_numerical`) |
| 7 | intentional safe failure plus resume test | met |
| 8 | isolated on the feature branch; research branch and jobs unaffected | met |
| 9 | gaps, vulnerabilities, licences, unimplemented items documented | met (`LIMITATIONS.md`, `SECURITY_AND_SAFETY.md`, ADR) |
| 10 | code visible on GitHub | see §7 |

## 7. Git

The push outcome is reported in the session summary. Only `feature/researchos-orchestrator` is pushed. Nothing goes to `main`, nothing is merged, nothing is force-pushed.

**Merging caution:** this branch was cut from `research-2026-10-13`, not from `main`. Merging it into `main` would also bring in the research documents. If ResearchOS should land on `main` without them, it should be cherry-picked onto a branch from `main`. That is the owner's decision.

## 8. Corrections to the request (counter-checked)

- **The "ten degree-pair features explain all gain" caution is correct, and the overstatement was mine.** `PUBLICATION_READINESS.md` at `8e5bc2e` lists "The gain is carried by 10 numbers" as a supported claim. Report 08 (N2) and report 07 §7 do mention 0.42 vs 0.47, but the readiness table does not. ResearchOS was told not to edit research documents, so the correction is left to the owner. The proposed wording is in the demo record.
- **"OrchestraOS" is not one product.** See the ADR. The request's description fits the Tulum-DAO project best, but that project is unfinished and is for coding agents.
- **Budget accounting for tokens and HPC** has nothing to measure in v0.1 (no model calls, no HPC). It was not built, rather than faked.
- **"Independent replication":** the E3 "independent recomputation" was a separate script over the *same* predictions, written by the same assistant. ResearchOS records it as `independent_recomputation` (different code, same inputs), not as independent data or an alternative implementation.
