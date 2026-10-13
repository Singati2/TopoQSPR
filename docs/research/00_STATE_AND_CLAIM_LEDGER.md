# 00 — State and claim ledger (Gate A), 2026-10-13

## Repository state (recorded 2026-10-13T13:36Z, refreshed 14:40Z)

| ref | SHA | note |
|---|---|---|
| remote `main` | `a0fd2a6` | `fecb800` + 7 docs-only commits (audit prompts, plans, this mission) |
| remote `v1.0-dev` | `fecb800` | |
| remote `audit-tests-20261012` (PR #2, open) | `5fc799e` | tests + permutation simulation; `test_a4` conflicts with the F2 fix (`pr2_test_a4.patch`) |
| tag `plan-freeze-v1` | `38eda80` | frozen plan (102 indices) |
| local `audit-2026-10-12` | `b36d591` | 10 commits on `fecb800`, **not pushed** |
| local `research-2026-10-13` (this work) | from `b36d591` | **not pushed**; working tree clean at each commit |

Package: `src/topoqspr/` (v1.0); the root `main.py`/`ml.py` is the v0.x prototype. The
E2 results of 2026-10-13 are on the local audit branch only (`docs/audit-2026-10-12/out/E2_*`).
They were re-verified from their own predictions (TEST_LOG rows 14–15). Nothing on GitHub
contains them.

## Paper claims

Claim-by-claim reconciliation with PASS/LIMITED/UNRESOLVED/FAIL:
`docs/audit-2026-10-12/CLAIM_PROVENANCE_LEDGER.md` (14 claims; 5 FAIL as worded, 4 LIMITED,
1 UNRESOLVED, 4 PASS). The confirmatory numbers are now re-derived from raw predictions (04).
Index counts per experiment: frozen confirmatory, Experiment C and E1/E2 use **102**; amended SVR
control, `main_scaffold`, `exploratory_families`, Experiment D and E3 use **104**.

## Mission and Oct-13-plan recommendations: verdicts

| # | Recommendation (source) | Verdict | Evidence / counterexample |
|---|---|---|---|
| R1 | "Unlabeled TI fail to distinguish isomers" must be narrowed to isomorphic graphs (plan §1.1) | **VALID** | 95.0% of Bradley within-formula pairs are non-isomorphic, and all 104 TI distinguish every one of them (05 §3) |
| R2 | Ethanol/DME as collision example (plan §1.2; mission §3) | **VALID** | 104/104 equal; NIST bp 351.5 vs 248–249.5 K |
| R3 | Collision reduction ≠ prediction improvement (plan §1.3) | **VALID** | 03 M4 C3, from our own data |
| R4 | Formula baseline constant within formula; no held-out formula means (plan §1.4) | **VALID**, implemented | E3b groups = formula; no target-derived formula feature |
| R5 | Novelty must be checked; labelled indices not new (plan §1.5; mission §8) | **VALID** | 02: E = ALREADY KNOWN |
| R6 | QM9 is computed, not measured (plan §1.6; mission §7) | **VALID** | QM9 not used. Measured Bradley chosen instead |
| R7 | Drug cohorts may be too sparse for isomer work (plan §1.7) | **VALID → NO-GO on drugs** | 05 §1 |
| R8 | Separate stereo / tautomer / salt artefacts (plan §1.8) | **VALID WITH CHANGES** | Stereo grouped by first block. Tautomers handled only through the standardiser. The salt/parent MP mismatch is a remaining risk (01) |
| R9 | Fit all preprocessing within folds (plan §1.9) | **VALID**, verified | asinh is fit-free; imputation and scaling are inside pipelines; ECFP is per-molecule |
| R10 | Small-data advantage is a hypothesis (plan §1.10, E5) | **VALID**; E5 **NOT RUN** | out of budget; not prioritised |
| R11 | Default framing question (mission §1) | **VALID WITH CHANGES** | The question is not new (Altairi et al., Sci Rep 2026, asked it and answered yes). The niche is the measured-endpoint, composition-controlled, nested, isomer-controlled test (02) |
| R12 | Title "What Does Molecular Topology Add? …" | **CONDITIONAL** | Only if E3 supports a clear statement either way; see 08 |
| R13 | Verify every index independently (mission §3) | **DONE** | 03a (implementation, 104/104), 03b (primary sources, 45/104 read), 03 |
| R14 | Recompute statistics from raw predictions (mission §6) | **DONE** | 04 |
| R15 | Bayes-risk statement E[Var(Y\|T)] (mission §3.3) | **VALID**, textbook | 03 M3 F3 |
| R16 | Cyclomatic check (mission §3.4) | **VALID** | `n_rings` = m − n + 1; size_rings is redundant |
| R17 | Hosoya overflow / float (mission §3) | **VALID concern; no effect** | 03 M5 |
| R18 | OPERA may contain predictions (mission §4) | **INVALID** for the files used | the values are measured PHYSPROP data (01b) |
| R19 | External validation (owner) | **VALID WITH CHANGES** | Bradley contains PHYSPROP, so the test was redesigned to be source-disjoint (06b R3) |
| R20 | Lean 4 formalisation (mission §3) | **NOT DONE** | no toolchain, and the statements are textbook (low added rigour) |
| R21 | GNN baseline (mission §7.9) | **NOT RUN** | not needed for the question; ECFP and RDKit 2D are the competitive baselines |
| R22 | Website destructive tests; clean-environment end-to-end reproduction (mission §5, §9) | **NOT RUN** | live site is v0.x (F7); clean rerun about 5 h |

## Gate A verdict

State recorded. Novelty: the question is known and contested. A niche remains for a measured,
composition-controlled, isomer-controlled, source-disjoint test. Proceed to Gates B–E.
