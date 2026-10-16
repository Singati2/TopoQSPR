# Agent contracts (v0.1)

`python -m researchos agents list` prints this table from code (`researchos/agents/__init__.py`).
Only the first three roles exist. The others are named interfaces with no behaviour.

| Role | Implemented | Inputs | Outputs | May | May not |
|---|---|---|---|---|---|
| **Coordinator** | yes (deterministic) | task file; upstream step results | `TaskSpec` DAG; claim classification | propose tasks; downgrade a claim to `inconclusive`/`under_test` | approve anything; set `supported`/`refuted`; run compute |
| **Verifier / reviewer** | yes (deterministic) | allowlisted committed CSVs; earlier recomputation files | consistency checks, learner audit, rule-based objections, fold-level intervals (if approved) | read; compute statistics from archived outputs | refit models; read outside the allowlist; contact a model or network |
| **Repository steward** | yes (deterministic) | Git (read-only argv allowlist), `/proc`, disk | artifact hashes and provenance, branch state, protected processes, health | read | commit, push, merge, delete, signal processes |
| Mathematics | no | — | — | — | — |
| Statistics | no | — | — | — | — |
| Experiment | describe-only | compute step params | protocol draft, resource estimate, input checks (dry-run) | describe | launch anything |
| Publication | no | — | — | — | — |
| LLM reviewer (Claude Code) | stub | — | "not enabled" | detect the `claude` CLI | send repository content anywhere |

## Adapter contract

```python
@adapter(name, scope, path_params=[...], role=..., requires=[tools])
def fn(ctx, params: dict, upstream: dict[str, dict]) -> dict:
    return {"status": "pass|fail|flag|not_run|error", "findings": [str], "data": {...},
            "evidence": [Evidence kwargs without claim_id]}   # evidence is attached to the task's claims
```

- Path parameters are allowlisted and hashed by the executor *before* the adapter runs.
- An exception inside an adapter is recorded as `error`. It is never a pass.
- An adapter that names a tool in `requires` is skipped (`not_run`, "executor unavailable") when
  that tool is absent.
- A new adapter needs code review and a test. A task file cannot define one.

## Escalation

The executor stops at, and lists in the report under "Awaiting the owner":
- any step whose scope has no approval for the current spec;
- every compute step;
- protected processes running while compute is requested;
- an input outside the allowlist.
