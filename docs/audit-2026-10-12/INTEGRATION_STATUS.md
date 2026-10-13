# Integration status (2026-10-12, late evening)

## Repository state (verified with `git fetch` + `gh`, 2026-10-12 ~21:40 EDT)

| ref | SHA | content |
|---|---|---|
| remote `main` | `ab766682825c3cb73a5d066313fed5a927222ed2` | `fecb800` + `c6d5d80` (independent audit doc) + `639c413` (CLAUDE.md) + `ab76668` (follow-up prompt) |
| remote `v1.0-dev` | `fecb800d2b10fb8802d889b6d0b4c3a192eef090` | unchanged since the 2026-10-12 trailer rewrite |
| remote `audit-tests-20261012` (PR #2, open) | `5fc799eff50b1f91202256588c2902bf33b3a6bc` | `639c413` + adversarial tests + permutation simulation |
| tag `plan-freeze-v1` | tag `f2ed9e5` → commit `38eda80` | frozen plan |
| local `audit-2026-10-12` (**not pushed**) | `fecb800` + 8 commits (as of `533dcbb`) | this audit |

The local audit branch is **not visible on GitHub**. Nothing from it has been pushed. Its base is
`fecb800`, so it does not contain remote `main`'s three docs commits. Those touch only new files
(`CLAUDE.md`, `docs/INDEPENDENT_ADVERSARIAL_AUDIT_2026-10-12.md`,
`docs/CLAUDE_FOLLOWUP_PROMPT_2026-10-12.md`), so a merge or rebase has no textual conflict
(checked with `git diff --stat`; no overlapping paths).

## PR #2 vs this branch

* `tests/test_adversarial_audit_20261012.py` run against this branch: 5 pass, `test_a4` fails
  by design. It asserts the old p = 0 contract that F2 changed. `pr2_test_a4.patch` makes all 6
  pass (TEST_LOG rows 9–10).
* Statistical justification for the changed contract (F2): with zero observed variance among
  the 25 (or J) correlated fold differences, the corrected SE is 0 by arithmetic, not by
  estimation. A p-value of 0 and a point CI would claim certainty about the generalisation
  difference that the resampling cannot support. NaN with an explicit reason is the honest
  output; the descriptive mean is kept. Downstream consumers: `evaluate.compare_feature_sets`,
  `experiments._ci_of_mean`, `experiments.experiment_d` (all propagate NaN; BH on NaN p is
  handled by `stats._family_size`). No archived cell reached the branch (checked over all
  paired groups of the four Experiment B runs).
* Merge order proposal (needs the user's authorisation): (1) apply `pr2_test_a4.patch` on PR #2,
  (2) merge PR #2, (3) merge or rebase `audit-2026-10-12` onto `main`, (4) run the full suite.

## Findings register

| ID | short | verdict | evidence |
|---|---|---|---|
| F0 | run-record commit IDs absent from GitHub after rewrite | provenance defect | `sha_map_2026-10-12.tsv` |
| F1 | catastrophes: Hosoya Z (TI sets) and RDKit Ipc (rdkit2d); paper names only Z | confirmed (mechanism); manuscript error | `out/catastrophic_main_confirmatory.csv` |
| F2 | corrected t: p = 0, point CI at zero variance | confirmed bug, fixed (`33e8198`) | `tests/test_stats.py` |
| F3 | abstract "1821 with measured endpoints" (959) | manuscript error | `targets.csv` |
| F4 | "preregistered" without registration | manuscript/integrity | `main.tex:38,265`, `supplement.tex:63` |
| F5 | scaffold logP "same size" | manuscript error | RESULT_RECONCILIATION §3a |
| F6 | Exp. D non-significance worded as equivalence | inferential limitation | §6 |
| F7 | live site is v0.x; `/health` hard-codes 1.0.0 | confirmed | `website.py:28,56` |
| F8 | README stale | doc defect | `README.md:11-22` |
| F9 | in-data collisions; all are graph isomorphisms; non-binding ceilings | verified | `out/collision_*` |
| F10 | epirubicin MP 344.53 °C | unverified (suspect value) | Bradley row 16461 |
| F11 | rdkit2d keeps Crippen MolLogP for measured logP | methodological risk | `target_providers.py:290` |
| F12 | controls ≠ selected learner | limitation (documented) | run records |
| F13 | "never imputed" vs in-model median imputation | doc defect | `models.py:11` |
| F14 | dataset licences (OPERA, Biogen as "MIT") | unverified | `expdata.py:219,231` |
| F15 | inexact CSV float parsing | confirmed (tooling) | TEST_LOG row 6 |
| F16 | SVR control from a dirty tree | provenance defect; source unchanged by circumstantial evidence; exact tree unrecoverable | run record, `git diff 1a504a9 ab2c30f -- src/ studies/` empty |
| F17 | permutation fallback for unique-size groups | real in general (0.48 at 0.05); 0 rows in archived runs | `scripts/permutation_block_check.py` |
| A1–A9 | PR #2 audit items | cross-walk in AUDIT_FINDINGS.md | — |

## Blockers needing the user

1. Authorisation to push `audit-2026-10-12` and to merge PR #2 (with the patch).
2. Render account for the v1.0 deployment (F7).
3. Journal guideline pages blocked to automated fetch. They must be checked in a browser.
4. Dataset licence confirmation (F14). Epirubicin value check (F10).
5. Decision on citing the new commit IDs plus the map in the paper (F0); purge of `related_papers/`
   (license risk; history rewrite) only with explicit authorisation and before the release tag.
