# ResearchOS v0.1 architecture

```
task file (.md)                        owner
  │  only the ```researchos-task JSON block is read;    │ approve <task> --scope … --confirm <task>
  │  prose = untrusted notes (screened, never run)      ▼
  ▼                                               approvals (bound to spec hash)
coordinator ──► TaskSpec (validated DAG of registered adapters) ──► store.tasks
                                                     │
                       dry-run: allowlist + hashes + approvals + health + protected processes
                                                     │
executor (lock per task) ── read-only steps ──► adapters (in-process, no shell)
        │                └─ analysis-only ──► only if approved for this exact spec
        │                └─ compute ────────► never launched in v0.1 (described only)
        ▼
step_results (keyed by input hash) ──► evidence ledger (append-only) ──► claims (status)
        │                                                              ▲
        └──► report.md (commit, sha256, findings, objections, approvals, limitations) ─┘
audit log (append-only) records every transition, approval, evidence row and claim change
```

## Modules

| module | responsibility |
|---|---|
| `schema.py` | Enums and dataclasses: `Claim`, `Evidence`, `TaskSpec`, `Step`; independence levels; canonical hashing |
| `statemachine.py` | Allowed task transitions: proposed → approved → queued → running → verifying → completed / failed / blocked. Invalid transitions raise |
| `store.py` | SQLite (WAL, `BEGIN IMMEDIATE` transactions), schema migrations via `user_version`; refuses a newer schema. Triggers make `evidence` and `audit` append-only |
| `config.py` | Repo discovery, state dir (`<repo>/.researchos`, git-ignored, or `$RESEARCHOS_HOME`), allowlist roots (default: the repo only), private paths (default `~/topoqspr_private`), protected process patterns, resource ceilings. Unknown config keys are rejected |
| `safety.py` | Realpath allowlist (blocks `..` and symlink escapes, and denies private paths), screening of untrusted text, redaction of `$HOME`, secret-named env values, token patterns and e-mails |
| `inspectors.py` | Read-only Git through a command allowlist (no push/commit/reset; `GIT_OPTIONAL_LOCKS=0`), `/proc` process discovery (excludes its own ancestors), CPU/RAM/disk/GPU health, tool detection |
| `evidence.py` | sha256 of artifacts, `SHA256SUMS` parsing (malformed or conflicting lines are errors), Git provenance (last commit, and committed blob = working file) |
| `adapters.py` | The only executable units. Each declares a scope and its path parameters |
| `executor.py` | Approval gate, dry-run planner, locked idempotent runner, crash detection and resume, claim downgrades |
| `agents/` | Role table and the coordinator's task-file parser |
| `report.py` | Markdown decision record |
| `cli.py` | `python -m researchos …` |

## Key design choices

- **Scopes, not trust levels.**
  - `read-only` adapters read allowlisted, committed artifacts in-process and run without approval. The owner can require approval with `read_only_requires_approval: true`.
  - `analysis-only` adapters need an approval recorded for the *same spec hash*: editing a task voids its approvals.
  - `compute` is never launched. A dry run gives the protocol draft, the resource estimate and the input checks.
- **Idempotency by input hash.** A step's identity covers its adapter, parameters, the sha256 of
  every input file, and the hashes of its upstream steps. A completed step with the same input hash is
  reused, never silently re-run. Changed data produces a new hash and a fresh result.
- **Crash safety.** One `flock` per task. If a run finds its task in `running`/`verifying` while holding the
  lock, the previous runner died. The task goes to `failed` ("interrupted"), then to `queued` with
  `attempt + 1`, and completed steps are reused.
- **Claims change only downwards automatically.** Automation can move a claim to
  `inconclusive`/`under_test`. A recommended refutation is stored as `inconclusive` with the
  recommendation in the audit log. Only the owner sets `supported` or `refuted`.
- **Evidence says how independent it is.** See `EVIDENCE_SCHEMA.md`. `qualitative_llm_review` is
  stored but flagged non-numerical.
- **Heavy outputs stay outside Git.** The ledger stores hashes and paths. Files listed in a manifest
  but deliberately absent from Git (predictions, OPERA rows) are reported as "absent by design", not
  as passes.
