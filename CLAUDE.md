# Instructions for Claude Code / Codex — TopoQSPR

The repository contains an independent, source-grounded adversarial audit dated 2026-10-12:

**[docs/INDEPENDENT_ADVERSARIAL_AUDIT_2026-10-12.md](docs/INDEPENDENT_ADVERSARIAL_AUDIT_2026-10-12.md)**

Read it **before** changing mathematical indices, cohort construction, nested cross-validation, statistics, experiment result tables, or the manuscript. Also read docs/EXPERIMENT_PLAN.md **including its amendments**, docs/CHECKPOINT.md, docs/FEYNMAN_REVIEW.md and docs/INDEX_VERIFICATION.md.

## New exploratory research proposal (countercheck first; 2026-10-13)

Read **[docs/POSITIVE_RESULT_FALSIFICATION_PLAN_2026-10-13.md](docs/POSITIVE_RESULT_FALSIFICATION_PLAN_2026-10-13.md)** when investigating whether graph-theoretic descriptors can show *genuine incremental predictive value* beyond molecular composition. **Do not accept its recommendations on authority.** Start by independently verifying or falsifying every premise (including mathematical small-graph counterexamples, prior-art novelty, data feasibility, chemical identity, leakage risks and metric validity) and write a recommendation-by-recommendation **VALID / VALID WITH CHANGES / INVALID / UNRESOLVED** ledger **before expensive runs**. Follow its gate-based protocol; retain negative results, freeze any new exploratory design before running it, and do not confuse a positive hypothesis with an observed result.

## Required workflow

1. Check the exact current Git SHA, read the frozen plan tag, and examine both original and amended Experiment B run records. **The main selected-model test is not the Ridge-only or SVR-only control.**
2. Investigate each audit item with a minimal counterexample or independent test. Challenge the audit too; mark a concern **not a bug** when the evidence proves it false. Never claim Codex independently checked anything unless Codex actually ran.
3. Add reproducible regression tests before modifying scientific calculations; run the relevant tests and, if feasible, the full suite. Never alter a published number silently.
4. Maintain a provenance ledger: issue → tested claim → file:line → test → code diff → rerun metrics → manuscript claim.
5. Do not change the original pre-specified four-endpoint confirmatory question retrospectively. Mark added composition-adjusted and chemical-disjoint analyses as exploratory unless an appropriate new protocol is prospectively frozen.
6. Keep known Ridge bootstrap failures, SVR amendments, four verified extrapolating molecules, and other negative findings visible. The paper's credibility depends on honest reporting.
7. Check preregistration language: a frozen Git tag is not necessarily a prospectively registered OSF protocol. Require dates and links for any 'preregistered' claim.
8. For each completed change, commit with the issue ID (A1–A9), update the audit's disposition and link tests and results. Do not claim a full audit or clean rerun on incomplete evidence.

## Immediate priority

A1 (selection-vs-control distinction), A6 (preregistration wording), A5 (102-vs-104 registry provenance), A3 (group permutation validity), A4 (degenerate corrected-t reporting), A2 (exploratory incremental gain beyond composition). See the linked report for precise source files, hypotheses, counterexamples and acceptance criteria.
