# Evidence schema (store v1)

SQLite at `<repo>/.researchos/state.db` (git-ignored; override with `$RESEARCHOS_HOME`).
Migrations are keyed by `PRAGMA user_version`. A store newer than the code is refused.

## claims

| field | meaning |
|---|---|
| `claim_id` | stable id, e.g. `C-MP-PRIMARY` |
| `text`, `scope` | exact wording and its scope (endpoint, design, data) |
| `status` | `proposed` · `under_test` · `supported` · `refuted` · `inconclusive` · `blocked` · `superseded` |
| `stage` | `exploratory` · `confirmatory` · `post_hoc` · `not_evaluated`. "Not evaluated" never becomes "failed" |
| `limitations`, `manuscript_refs` | lists |

## evidence (append-only)

`evidence_id` = `E-` + the first 12 hex digits of the sha256 of the canonical record (identical evidence deduplicates).

| field | meaning |
|---|---|
| `claim_id`, `kind`, `summary` | what it is about |
| `supports` | `supports` · `contradicts` · `neutral` |
| `independence` | see below |
| `reviewer`, `method` | who or what produced it, and how (script path, rule, person, model) |
| `artifacts` | `[{path, sha256, commit}]` |
| `metrics` | numbers with intervals |
| `commit`, `command`, `seed`, `protocol_hash`, `run_id` | provenance |
| `dataset`, `license_note`, `limitations` | data terms and caveats |

### Independence levels (weakest first)

| level | numerical? | meaning |
|---|---|---|
| `qualitative_llm_review` | **no** | a model read something and gave an opinion. It gets an automatic limitation |
| `owner_assertion` | **no** | the owner states it |
| `same_code_rerun` | yes | identical code re-executed |
| `consistency_check` | yes | archived files agree with each other (arithmetic only) |
| `independent_recomputation` | yes | different code over the same inputs |
| `alternative_implementation` | yes | an independent implementation of the method |
| `source_disjoint_data` | yes | no shared records; independence at the reference level **not** shown |
| `independent_data` | yes | independence shown at the reference level |
| `external_source` | yes | a published primary source |

These labels never upgrade themselves. "Source-disjoint" stays `source_disjoint_data`,
"exploratory" stays `exploratory`, and a second model's "verified" is `qualitative_llm_review`.

## tasks, step_results, approvals, audit

- `tasks`: `task_id` = `T-` + 12 hex digits of the sha256 of the canonical spec (prose excluded), plus spec JSON, untrusted notes, state and attempt.
- `step_results`: keyed by (task, step, input hash), with status and result JSON.
- `approvals`: (task, scope), plus spec hash, approver, method and time.
- `audit` (append-only): every transition, approval, evidence row and claim status change, with actor and reason.
