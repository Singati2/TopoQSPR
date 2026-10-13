# Claude Code Prompt — TopoQSPR Audit Integration and Publication Readiness

**Prepared:** 2026-10-12  
**Repository:** `Singati2/TopoQSPR`  
**Purpose:** Turn the October 12 adversarial audits into a reproducible, scientifically honest release and journal submission. **This is an instruction prompt, not evidence that the audit branch has been pushed or that any experiments have been rerun.**

## Role and non-negotiable rules

You are acting as an adversarial computational-science auditor with expertise in mathematical graph invariants, cheminformatics, experimental design, statistical inference, software verification, and reproducible research. Try to **falsify** each claim. Never equate a clean test suite with a valid scientific inference.

1. Inspect the **current remote HEAD**, branches, open PRs, diffs, tags, frozen plans, run records, and working-tree status before touching code. Record complete SHAs.
2. Read `CLAUDE.md`, `docs/INDEPENDENT_ADVERSARIAL_AUDIT_2026-10-12.md`, `docs/EXPERIMENT_PLAN.md` (including amendments), `docs/CHECKPOINT.md`, `docs/FEYNMAN_REVIEW.md`, `docs/INDEX_VERIFICATION.md`, `docs/EVALUATION_PROTOCOL.md`, and the paper and supplement. If present locally, also read `docs/audit-2026-10-12/` and its `AUDIT_FINDINGS.md`, `RESULT_RECONCILIATION.md`, `SCIENTIFIC_RESEARCH_PLAN.md`, `TEST_LOG.md`, and `sha_map_2026-10-12.tsv`. **Do not assume the local audit branch is visible on GitHub.**
3. Distinguish **confirmed bug**, **possible inferential limitation**, **documentation/provenance defect**, **unverified claim**, and **new exploratory experiment**. For every confirmed bug give a minimal counterexample, failing-before/passing-after regression test, exact source location, and quantitative effect.
4. Preserve original frozen analyses, negative results, controls, and full audit trail. Do not change archived results, retcon a prospective protocol, delete inconvenient molecules, force-push, rewrite history, or merge PRs without inspecting impacts. For any new experiment, create a dated, immutable amendment **before** execution; label it exploratory unless genuinely prospective.
5. Do not say “verified”, “independent”, “preregistered”, “significant”, “equivalent”, “no effect”, or “published result unaffected” unless evidence specifically supports that wording. A confidence interval containing zero does not prove equivalence or zero gain.

## A. Integrate the audit without destroying provenance

The latest *user-supplied Claude run report* says five commits on local `audit-2026-10-12` document 17 findings and analyses, while remote `main` already has `CLAUDE.md` and an independent audit and PR #2 (`audit-tests-20261012`). **Verify this state; it may have changed.**

- Reconcile local audit work with remote main and PR #2 using non-destructive diffs. If the local branch is still private to a workstation, explicitly report that and request access or publish it only when authorized. Do not invent GitHub links for local-only files.
- Run PR #2’s adversarial tests against the proposed fixed statistics implementation. A reported conflict is `test_a4_degenerate_corrected_ttest_behavior_is_explicit`: the original assertion expects `p=0` and a zero-width CI for `[0.1,0.1,0.1]`; the audit’s proposed corrected contract makes the inferential p-value and CI **undefined** while preserving the descriptive mean. Confirm that the revised behavior is statistically justified, and update the test only after verifying the implementation. The local audit reportedly has `pr2_test_a4.patch`; inspect it rather than assuming its correctness.
- Audit the report’s `sha_map_2026-10-12.tsv` against `run_record.json`, `paper/main.tex`, and `\\runCommit`. Broken historical IDs must be fixed using a transparent old-to-new mapping or immutable archival reference—not silent substitution with current HEAD.
- One amended SVR control run reportedly has `dirty=true`, with no captured working-tree diff. Distinguish circumstantial evidence about untouched `src/` and `studies/` from exact reproducibility. Require clean pinned rerun, or mark this control’s precise source tree unrecoverable.
- Confirm the original confirmatory analysis used **102** indices while software/amended analyses may use **104**, and the main selected-model test differs from Ridge-only or SVR-only controls. Make source-to-table version labels unambiguous.

## B. Test the most consequential statistical claims

1. **Beyond composition:** A reported exploratory asinh-sensitivity study tested indices + element composition against composition alone on identical folds and found no 95% interval excluding zero for four endpoints, either random or scaffold splits. Reported melting-point paired gains: random `+0.05 [-0.01, +0.12]`; scaffold `+0.09 [-0.01, +0.19]`. Independently reproduce fold definitions, train-only transforms, model selection, paired contrasts, corrected-CV uncertainty and multiplicity. Never write “indices add nothing” as a proven equivalence conclusion. Identify whether asinh was applied inside train/test preprocessing without target leakage, and distinguish descriptor transform from index definition.
2. **Original logP headline:** Confirm the frozen size-baseline contrast `+0.15 [0.08, 0.23]`, and reported scaffold contrast `+0.10 [-0.05, 0.25]` using the original **102** indices. Fix any prose saying the latter is “the same size but wider”: it has a smaller point estimate and uncertainty includes zero.
3. **Catastrophic extrapolation:** Produce per-molecule, per-fold attribution of extreme errors: Hosoya Z for indicated large drugs including lanatoside C, atracurium and digoxin/digitoxin as applicable, and RDKit `Ipc` for the RDKit descriptor set. Confirm the exact endpoint attribution before altering any prose. Show training-range violations and the learner/prediction paths; keep all frozen primary-cohort molecules.
4. **Permutation exchangeability:** The reported group-size fallback can shuffle observations within a unique-size correlated cluster, leading to severe anti-conservative Type I error under some synthetic nulls; the supplied report says a simulation rejected `96/200` times at nominal `0.05`. It also reports **zero fallback rows across the archived main/SVR cells**, all multi-member groups being same-size pairs. Independently verify both claims from real archived run-stage group arrays; assess all other exchangeability assumptions. Add a fail-closed warning or group-respecting alternative for future runs, and test calibration under several nulls. Do not claim the *entire* permutation procedure is universally valid simply because this specific fallback was not triggered.
5. **Degenerate corrected t-test:** Establish and document a defensible contract for constant nonzero paired fold differences; no p=0/point CI presented as population certainty. Check all downstream consumers and whether any published cells actually use the edge case.
6. **Representation collisions:** Validate reported examples (including ephedrine/pseudoephedrine) by recomputing all indicated descriptors, compositions, stereochemistry, target values and data provenance. Re-derive any claimed performance ceiling mathematically; verify direction of the inequality (an upper ceiling is `R^2 <= bound`, not `R^2 >= bound`), distinguish irreducible label collisions from finite-sample learning error, and quantify how common collisions are.

## C. Run the missing, predeclared exploratory benchmark only after amendment

Write a dated amendment specifying the exact cohorts and endpoints, IDs and exclusions, matched outer folds, grouped inner CV, seeds, learner selection, transforms, hyperparameter budget, corrected uncertainty, family of hypotheses, multiplicity correction, compute budget and stopping rule. **Do not optimize the design using observed results.**

- Paired **composition-only** vs **composition + 102/104 explicitly versioned indices**, plus **ECFP/Morgan fingerprint** and ordinary RDKit 2D baselines, with comparable optimization and held-out predictions.
- Verify one chemically disjoint **fingerprint-similarity cluster split** in addition to random and Bemis–Murcko scaffold splits. Prevent stereoisomer, parent/salt, tautomer, and duplicate leakage where the intended scientific estimand demands it; explicitly document remaining similarities.
- Optional isomer-controlled study only if the data actually support defensible within-composition comparisons. Report effective sample size and uncertainty.
- Archive raw fold predictions and splits, environment lock, descriptor manifests/hashes, failure handling, CPU/GPU requirements, wall time, and the exact entrypoint. Run once under the frozen exploratory amendment; unfavorable outcomes remain results.
- If the interval supports a narrow melting-point gain, describe **only that**. If no convincing gain remains, frame an informative reproducible negative benchmark. If fingerprints outperform indices, say so.

## D. Correct scientific writing and publication integrity

Check and correct (with evidence, not a global blind search/replace):

- “Preregistered” vs Git-frozen/prespecified plan: the report says no prospective OSF registration was filed. Preserve dates and never retroactively call this preregistration.
- The abstract reportedly says **1,821 drugs with measured endpoints**, while **959** reportedly have at least one endpoint measurement. Recalculate independently: 1,821 is the candidate/parent cohort, not necessarily the measured subset.
- Distinguish non-significant differences from **equivalence**; supply an explicit equivalence margin and suitable interval test before making any equivalence claim.
- Selected learner versus Ridge/SVR family-specific controls; 102/104 descriptors; size-only versus composition; random versus scaffold inference; measured versus estimated endpoints; hypothesis family and BH q-values.
- ChemGraphX coauthor overlap and related work: describe the relationship fairly, check licences, attribution, originality, and journal disclosure requirements. Never purge a tracked reference PDF or rewrite Git history without identifying the license risk and obtaining explicit authorization.
- Check current *official* author guidelines, open access fees, page/word limits and scope for **Journal of Cheminformatics**, **Digital Discovery**, and **SoftwareX**. Prefer a scientifically honest research article if the completed benchmark is substantive; otherwise consider a constrained software paper. Do not assert journal acceptance or acceptance probability.

## E. Software/release hardening and reproducibility

- Compare all paper tables, figure labels, statistics, source citations, CSV data, endpoint coverage and machine-readable run records through a claim-provenance ledger.
- Run focused regression tests, the entire test suite, and end-to-end reconstruction that is feasible on the available machine; report exact commands, duration, exit status, skipped tests, environment, and differing hashes. The supplied run report claims 1,103 tests before and 1,104 after the statistical fix; treat this as **reported**, not independently confirmed.
- Check deployed site versus repo: the supplied audit says the public Render deployment still serves v0.x while `/health` hard-codes `1.0.0`. Verify independently and either correct version reporting or plan an authorized v1.0 deployment; do not deploy automatically.
- Make future confirmatory runs refuse dirty trees or archive complete diffs with run records. Enforce data-source/chemical-data licensing and citation requirements before a public release.

## Required deliverables

Produce/maintain:

1. `docs/audit-2026-10-12/INTEGRATION_STATUS.md` — branch/PR/commit comparison, findings F1–F17 and A1–A9, disagreements and evidence, blockers.
2. `docs/audit-2026-10-12/CLAIM_PROVENANCE_LEDGER.md` — manuscript claim -> data/cohort -> run SHA/config -> exact test -> CI/effect -> source paths -> permitted wording.
3. `docs/audit-2026-10-12/EXPERIMENT_AMENDMENT.md` — prospective *for the new exploratory run* (not a retroactive preregistration), locked before execution.
4. `docs/audit-2026-10-12/RESULTS_DIFF.md` — old/new values and checksum differences; explicit “none” if no values changed.
5. A concise `PUBLICATION_READINESS.md` with journal strategy, claim scope, unresolved data/licence issues and reproducibility checklist.

Use filenames that do not overwrite already existing audit artifacts without careful merging; adjust destinations as necessary.

In your final response, provide a table with **finding | verdict | minimal test/counterexample | code/manuscript path | fix commit | tests executed | result change | remaining limitation**. Show the actual `git status`, remote branch, PR state, pytest summary, paper build status, and links **only to files actually pushed**. Separate completed work from proposals. **No push, merge, force-push, deployment, history rewrite, or frozen-result replacement without explicit user authorization.**
