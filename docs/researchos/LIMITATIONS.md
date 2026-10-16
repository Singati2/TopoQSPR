# Limitations and not-implemented items (v0.1)

**Not implemented (no code, despite appearing in the original request):**
- web dashboard;
- SLURM/FAU HPC adapter (none written; `doctor` only reports whether `sbatch` exists, and it does not on this workstation);
- the mathematics, statistics and publication agents;
- the Claude Code / LLM adapter (stub only, never called);
- budget accounting for tokens and API spend (there is none to account for);
- compute execution;
- signed approvals;
- a multi-user approval UI;
- `researchos task approve` from a phone or web.

**Known limits:**
- The approval gate is a guard against accidents and injected instructions, not a security boundary
  against someone with local shell access (see SECURITY_AND_SAFETY.md).
- The verifier re-derives archived numbers. It cannot detect an error that is consistently present in
  both the archived predictions and the summaries computed from them. That needs a re-run or an
  alternative implementation.
- `fold_level_contrast` uses tabulated t quantiles for 4, 9, 14 and 24 degrees of freedom only (it
  errors otherwise). It bounds the *selected-model* gap; it does not remove learner confounding.
- `learner_family_audit` expects the E3 CSV layouts (`test,set,…,family,r2` and
  `endpoint,estimand,set,…,r2_fold_mean,family_counts`).
- Process protection matches command-line substrings only.
- Static analysis used: `compileall`, `-W error` import, `git diff --check`, an AST unused-import
  scan. No ruff or mypy is installed in the project environment, so none was run.
- The rule-based reviewer only raises the objections it was written to raise. It is a checklist, not
  a reviewer.

**Next milestones**, in priority order. Scope is given as complexity only; there is no cost basis to quote:
1. Fixed-learner protocol integration: the owner freezes the draft from the worked example, then a
   `compute` executor with an explicit owner token, resource ceiling and protected-process wait.
   Medium.
2. Signed or second-factor approvals and an approval inbox (`researchos approvals list`). Small to medium.
3. An LLM reviewer adapter that sends only owner-approved, redacted excerpts and records the output
   as `qualitative_llm_review`. Medium.
4. A SLURM adapter, dry-run by default, written against the real FAU cluster documentation once
   access exists. Medium. Its configuration must not be guessed.
5. A read-only dashboard over the same store. Medium. Only after the CLI has been used in practice.
