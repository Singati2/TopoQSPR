# Security and safety (v0.1)

## Guarantees, and the test that covers each

| Guarantee | Mechanism | Test (`tests/researchos/test_researchos.py`) |
|---|---|---|
| An unapproved analysis step never runs | scope gate, approval bound to spec hash | `test_unapproved_analysis_and_compute_never_execute`, `test_approval_bound_to_spec_hash` |
| Compute is never launched, even if approved | executor never calls compute adapters; the adapter body also raises | `test_compute_not_launched_even_when_approved` |
| Approval needs explicit confirmation | task id typed at a TTY, or `--confirm <task_id>` | `test_approval_needs_exact_confirmation` |
| Text in documents, logs or task prose cannot approve or run anything | only the JSON block is parsed; adapters are referenced by registered name; prose is screened and flagged | `test_untrusted_text_cannot_grant_approval` |
| Paths outside the allowlist, `..`, symlink escapes and private paths are blocked | realpath check plus deny list | `test_path_outside_allowlist_and_symlink_escape_blocked`, `test_private_path_denied_even_if_allowlisted` |
| A changed or missing input blocks downstream steps | sha256 against manifest; committed blob = working file | `test_hash_mismatch_blocks_downstream`, `test_missing_file_reported`, `test_uncommitted_change_detected`, `test_manifest_conflicting_hashes_rejected` |
| No double execution | per-task `flock`; input-hash reuse | `test_duplicate_run_refused_while_locked`, `test_crash_then_resume_without_double_execution`, `test_completed_task_is_not_rerun` |
| Live jobs are observed, never touched | `/proc` reads only; no signals, ptrace or attach | `test_live_protected_process_is_observed_not_touched` |
| Git is read-only | argv allowlist; `hash-object -w` and `remote add` refused | `test_git_wrapper_is_read_only` |
| Evidence and audit cannot be edited or deleted | SQLite triggers | `test_evidence_and_audit_are_append_only` |
| An unknown state fails closed | newer schema refused; unknown config keys refused; outside Git → exit 2 | `test_newer_store_schema_refused`, `test_config_rejects_unknown_keys`, `test_cli_outside_git_fails_closed` |
| No secrets or private paths in reports | redaction of `$HOME`, secret-named env values, token patterns, e-mails | `test_redaction` |
| Works with no LLM and no network | deterministic adapters; the LLM adapter reports "unavailable" | `test_unavailable_executor_fails_closed`, `test_cli_offline_without_claude` |

## Threat model and known limits

- **Not a sandbox.** The gate protects against accidental execution and against instructions that
  arrive through content (prompt injection). It does not stop a person, or an agent with a shell, from
  editing the SQLite file or running Python directly. Approvals record the OS user name; they are not
  cryptographically signed. A signed-approval design is listed in `LIMITATIONS.md`.
- **Screening is heuristic.** The untrusted-text scanner flags common patterns (`rm -rf`, force
  push, pipe-to-shell, "ignore previous instructions", "owner approved"). It is only there to make
  attempts visible. Safety does not depend on it, because nothing in prose is ever executed.
- **Process matching is by substring** of `/proc/<pid>/cmdline`. A protected job that runs under a
  different name is not recognised. Default patterns: `test_baselines.py`, `baselines.py`,
  `bp_tests.py`, `box_sensitivity.py`, `e3_isomer.py`, `bin/topoqspr `, `-m topoqspr`,
  `joblib.externals.loky`.
- **Data handling.**
  - Reports contain hashes and aggregate numbers only.
  - OPERA per-molecule predictions and the owner's private PDFs are never read into reports. `~/topoqspr_private` is on the deny list by default.
  - One committed manifest, `docs/research/out/E3/mp_bradley_E3a/SHA256SUMS`, lists absolute paths that include the user's home directory. That is a minor privacy leak in already-pushed history; ResearchOS redacts it from its own output and verifies the same-named files beside the manifest.
