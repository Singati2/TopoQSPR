# TopoQSPR — Independent Publication Assessment and Claude Research Priorities

**Review date:** 2026-10-16 (America/New_York)  
**Review basis:** GitHub `Singati2/TopoQSPR`, branch `research-2026-10-13` at `9dce4d3291e1aeb9bbf6fd70d10e5d45b316d11d`, and `feature/researchos-orchestrator` at `5fea2fa0f1788f3f2baf6ca3d8fb279b921d1f24`.  
**Nature of review:** Independent *critical reading* of checked-in reports and commits. **I did not re-run experiments or inspect the private full prediction archives.** The reported "independent recomputation" uses alternate code on the same predictions; it is not independent-data replication.

## Executive decision

**CONDITIONAL GO for one narrowly framed empirical cheminformatics research article; NOT READY TO SUBMIT.**

A defensible core result is that *unlabeled molecular topology adds measured melting-point (MP) predictive information beyond elemental composition under a formula-controlled evaluation, while adding little or no practical value beyond strong labeled chemistry descriptors*. The corresponding controlled `logP` result is practically null. This is a useful endpoint- and representation-dependent evaluation, **not** proof that topological indices generally improve QSPR, a new index, or a new theorem.

The best-fit provisional venue is a **Journal of Cheminformatics research article**, with **Digital Discovery** as a possible alternative if the mechanistic/generalizable insight becomes stronger. These are proposed targets, **not** acceptance predictions. **SoftwareX** may fit a *separate*, mature software paper after release-quality work; it is not the primary destination for this scientific result.

**Highest-value next work:** Resolve the full-index-versus-edge-count learner confounding; then conduct genuinely prospective confirmation on an untouched measured source/endpoint. Avoid adding arbitrary new indices or fishing for positive endpoints.

## Evidence snapshot (as archived, not rerun in this review)

Source: `docs/research/07_RESULT_RECONCILIATION.md`, `PUBLICATION_READINESS.md` and associated E3 small outputs.

| Experiment and comparison | Reported result | Interpretation / qualification |
|---|---:|---|
| MP E3a, known formula / unseen constitution: TIa + composition vs composition | **ΔR² = +0.085**, 95% CI **[+0.065, +0.106]**; MAE **4.8 °C lower**, 15/15 folds positive | Strongest archived positive contrast; exploratory, data-driven endpoint choice |
| MP E3b, held-out formula: same comparison | **+0.106 [0.087, 0.125]** | Supports transport across formula groups in same measured dataset |
| MP source-disjoint external, known formula | **+0.072 [0.019, 0.138]** | Single-fit external contrast; bootstrap covers test sampling |
| MP source-disjoint external, unseen formula | **+0.240 [0.203, 0.281]** | **Do not oversell**: composition alone R² ≈ 0.23 vs topology+composition ≈ 0.47; ECFP ≈ 0.57 and RDKit 2D ≈ 0.67 |
| logP E3a, known formula: TIa + composition vs composition | **+0.016 [0.001, 0.032]**, BH-adjusted q = **0.056** | Not significant under the stated multiplicity correction; practically small |
| logP E3b, held-out formula | **−0.001 [−0.023, +0.020]** | Practically null under the reported ±0.05 equivalence bound |
| MP E3b, TI added to ECFP | **−0.015 [−0.029, −0.002]**, q = **0.048** | Additional TIs may hurt this stronger representation |
| Within-formula MP, centered R² | TIa + composition ≈ **0.18–0.19**, labeled representations ≈ **0.52–0.59** | **Post hoc** centered metric. The original pre-declared within-formula R² is negative for the TI representations. Report both |
| E3a MP negative control, within-formula target permutation | **−0.028 [−0.052, −0.004]** | Did not produce a spurious positive gain; a control, not a fresh independent test |
| MP sensitivity, replicate filter disabled | **+0.088 [0.052, 0.123]** | Supports robustness to this filter, on the reported seed-0 sensitivity run |

Reported MP development sample: **6,725 molecular parents / 758 formula groups**; logP: **3,994 parents / 495 formula groups**. Three seeds × five outer folds are documented for the nested evaluation. E3 execution, 42-contrast BH family, and same-prediction recomputation were reported completed on 2026-10-16. Before manuscript submission, independently inspect exact scripts, fold/group assignments, archives, hashes and statistical estimands rather than trusting this table.

## Contribution ledger: what is new, what is not

### Plausible and publishable *if carefully scoped*

1. **N1 (primary):** A controlled, large measured-MP comparison of element counts versus unlabeled molecular topology, with formula-aware splits, paired uncertainty, controls, and source-disjoint external evaluation. The design is a stronger contribution than claiming novel descriptors.
2. **N3 (supporting):** Isomer discrimination by a descriptor family is **not** the same as usefulness for measured-property prediction: distance-based indices can distinguish more graphs while predicting MP worse under the current model-selection pipeline. Check fixed-learner alternatives before attributing this solely to the information content of the indices.
3. **N4 / N7 (scope):** Topology has a measurable contribution relative to the deliberately weak but scientifically interpretable *elemental-composition* reference for MP, **not** to ECFP or RDKit 2D; the `logP` contrast is effectively null. This constrains any generalization across endpoints.

### Known, incremental, or presently unresolved

- Degree-based topological indices as linear functions of degree-pair edge counts are **already known** (see the project's prior-art matrix and cited Rada/Gutman literature); **do not claim a new theorem**.
- Low rank/redundancy, descriptor collision limits, and superiority of richer chemical descriptors are known in substance. The registry-specific facts (**104 indices, rank 38; 81 reducible within the specified degree conditions**) are useful checks, not a headline mathematical discovery.
- **N2 remains unresolved in its strongest form.** In internal cells and external known-formula evaluations, nine degree-pair counts plus maximum degree closely match the full TI set. But for the **external unseen-formula** test, the archived absolute R² is **0.471 (full TIa+composition) vs 0.423 (edge types+composition)**. Different selected learners (**extra trees vs SVR**) and missing paired uncertainty prevent the blanket claim that "10 numbers carry all predictive gain." The October 16 research-branch correction rightly narrows this.
- The companion `Singati2/topological-index-orthogonality` manuscript overlaps strongly in the edge-count span, redundancy and descriptor-ceiling theme. Coordinate scope, authorship and citations **before** pursuing parallel submissions.
- No claim of disease-specific predictive validity follows from the cohort search or E2 work.

## Priority work plan for Claude — do not blindly follow this advice

**Meta-instruction: independently challenge every recommendation in this file.** Check the archived protocol, actual code, source literature, data provenance, outcome definitions, and whether a proposed experiment really separates alternative explanations. Mark each item `VALID`, `MODIFY`, or `REJECT` with evidence and rationale. Never optimize toward a desired positive result; preserve and publish negative and inconclusive findings.

### P0 — Reconcile the full-index vs edge-count gap (highest immediate priority)

**Question:** Does the +0.048 absolute-R² external unseen-formula difference reflect extra non-degree information, learner-family selection, or test variation?

1. **Analysis first, no new training:** recover all fold-level paired scores and the external predictions from authorized local archives; verify SHA256 provenance; calculate a clearly labeled paired contrast and appropriate formula-cluster uncertainty where raw predictions and group IDs are available. Note this *does not* eliminate learner confounding.
2. Inspect the exact selected learners and hyperparameters for each feature set, fold and source. Report absolute R², MAE, sample sizes, selected family and feature preprocessing.
3. Propose a **fixed-learner controlled experiment** using the *same* learner families (at minimum SVR and extra trees), identical train/test groups and preprocessing rules, nested tuning **within training only**, and paired inference over a defensible grouping unit. Permit fair within-learner hyperparameter tuning, but never select by external test performance.
4. **Freeze and record** estimand, dataset IDs/hashes, grouping, tuning budget, seed plan, practical threshold, correction family, analysis code and stop rule **before fitting**; obtain owner approval and resource check before computation. The result can be positive, zero or negative.
5. The decision must distinguish **(a)** equivalent predictive performance within ±specified margin, **(b)** evidence that other descriptors help, and **(c)** unresolved uncertainty. Do not force a binary conclusion.

### P1 — Obtain genuinely prospective confirmation (not just code recomputation)

The E3 hypothesis and melting-point endpoint were selected after inspecting earlier E1/E2 results. A protocol committed before E3 fitting is transparent and valuable, but **is not third-party preregistration nor independent replication**.

- Identify an **untouched, appropriately licensed measured** endpoint/source with credible molecule identity, units, salt/tautomer/parent handling, and overlap tracing; an independent MP collection is ideal. A packing-sensitive second endpoint is another option, but feasibility and prior exposure matter more than an attractive name.
- Perform a **target-blind source and overlap census first**. Record selection criteria before reading targets or scores.
- Freeze primary contrast `TIa+composition vs composition`, formula-group splitting, external holdout, metrics, minimum effect, sample/exclusion rules, seed/power justification, and multiplicity plan before fitting.
- Retain RDKit 2D and ECFP as full-strength competing references. Test whether any topological advantage remains when chemically labeled information is included; a null result is informative.
- **No invented "independent validation"**: document if reference-level source overlap is unknown. In particular, do not assume OPERA has demonstrated independent underlying measurements merely because files come from different sources.

### P2 — Make the statistical reporting reviewer-resistant

- Report the **negative pre-declared** within-formula R² beside the **positive post hoc centered** variant, with precise estimands and reasons they differ.
- Keep the E3 analysis **exploratory**; distinguish pre-fit dated amendments, exploratory secondary contrasts, sensitivity checks and genuinely new confirmatory tests.
- Show absolute performance next to every ΔR², with uncertainty and practical equivalence margins; do not interpret a nonsignificant p/q alone as proof of equivalence.
- Audit repeated-CV dependence, corrected-t assumptions, formula-group bootstrap coverage, label-source clustering and multiple-comparison families. Where assumptions are fragile, give robust sensitivity analyses with limitations.
- Probe artifact risks: duplicate structures, isomorphic unlabeled graphs with distinct atom labels, cross-source overlap, salts, tautomer handling, zwitterions/decomposition temperatures and endpoint measurement noise.
- Distinguish a property of **feature information** from a property of the **feature × learner × training budget** combination; the present nested-CV estimand is best-in-class performance per representation.
- Include fixed-learner ablations and available failures, not only averaged winning models.

### P3 — Lock the novelty boundaries and manuscript architecture

Recommended working title: **"What Does Unlabeled Molecular Topology Add? A Formula-Controlled Evaluation of Melting-Point and Lipophilicity Prediction."**

A strong paper arc:
1. **Question:** What does *unlabeled topology* add over *elemental composition* in measured QSPR, when evaluation is grouped by formula?
2. **Design:** Explain why composition is a meaningful isolation control, but plainly acknowledge that it is weaker than chemistry-rich fingerprints/descriptors.
3. **Evidence:** MP gain under E3a/E3b/source-disjoint test; `logP` practical null; negative control; sensitivity and external caveats; absolute baselines and uncertainty.
4. **Mechanism / limits:** Degree-pair versus full-index attribution **only within evidence**, discrimination vs prediction, information-losing unlabeled graph collisions, inability to capture labeled chemistry/packing.
5. **Interpretation:** Clearly distinguish mathematical prior art, supported measurements, unexpected observations, unresolved questions and post hoc analyses.
6. **Reproducibility:** Data provenance, licences, versions, pinned environments, archived hashes, exact splits, rerun scripts, tables, and a claims-to-evidence ledger.

**Do not submit parallel overlapping mathematics and TopoQSPR papers** without a defensible division of contributions and explicit cross-reference. Ask the coauthors to decide: merge some findings, or keep this paper empirical and cite the orthogonality work.

### P4 — Repository/release/licensing gate

- Confirm the copyrighted article page-photo PDF is gone from branch tips; it **remains in Git history** according to `11_HISTORY_CLEANUP_PLAN.md`. **No force push, history rewrite, branch deletion, visibility change or purge without the owner's explicit approval.** Protect the private research copy.
- Resolve the paper/data licensing questions, including ChEMBL attribution, source licences and OPERA redistribution terms. Do **not** commit private per-molecule OPERA data, publisher PDFs or sensitive local files.
- Finish software release claims separately; the existing website must not be represented as a deployed, verified v1.0.
- Maintain a dated release candidate and reproducibility checklist only after the scientific claim and legal issues are resolved.

### P5 — ResearchOS: a useful assistant, not experimental evidence

The branch `feature/researchos-orchestrator` reports a v0.1 orchestrator with **37 tests passed** and a worked audit that caught the N2 overclaim. This supports traceability, not independent scientific replication. The LLM adapter, execution adapters for FAU HPC/SLURM, agent network and web UI are not implemented. Do not divert the publication effort into expanding orchestration features until P0/P1 are resolved.

## Required deliverables from Claude

Create **new, reviewable** files on a non-`main` branch (do not overwrite authoritative reports without an explained change record):

1. `docs/research/13_RECOMMENDATION_COUNTERCHECK.md`: each recommendation above marked `VALID/MODIFY/REJECT`, with primary-source citations, code evidence and any errors in this document.
2. `docs/research/14_FIXED_LEARNER_PROTOCOL.md`: falsifiable P0 plan, exact estimand, groups, model comparisons, paired intervals, budget, potential outcomes, frozen stop rules and explicit owner-approval checkpoint.
3. `docs/research/15_CONFIRMATORY_REPLICATION_FEASIBILITY.md`: target-blind candidate datasets, overlap/licensing assessment, required sample size / achievable power, and go/no-go rationale.
4. `docs/research/16_PAPER_CLAIM_LEDGER_AND_OUTLINE.md`: exact allowable/forbidden wording, novelty overlap ledger, abstract draft and paper outline (state all uncompleted experiments as pending).
5. `docs/research/17_VERIFICATION_REPORT.md`: what was actually read, rerun, independently recomputed, checked from hashes, failed or not attempted; machine/branch/commit identities.

**Operating guardrails:** Do not stop or disturb running experiments (including any separately launched `test_baselines.py`). Do not secretly retrain, cherry-pick results, change frozen protocols mid-run, create unauthorized GPU/HPC load, change Git history, merge to `main` or push there. Use isolated branches/worktrees and ordinary non-force commits. Before costly computation or destructive changes, obtain explicit owner approval. Push review files to the authorized research branch and provide full links/SHAs, distinguishing local work from work actually on GitHub.

## Final reviewer judgment

TopoQSPR now has a credible **moderately novel, carefully constrained empirical story** and materially better reproducibility than its earlier exploratory version. It **does not yet establish a general mathematical contribution or topological-index superiority over modern chemistry descriptors**. Its publication strength will be decided by rigorous fixed-learner attribution, genuine confirmation and transparent writing of both positive MP and practically null `logP` findings. If P0/P1 weaken the headline, retain and publish the scientifically informative negative answer rather than tuning toward a positive one.

**References inside this repository:**
- `docs/research/07_RESULT_RECONCILIATION.md`
- `docs/research/08_PAPER_NOVELTY_AND_SUBMISSION_DECISION.md`
- `docs/research/09_REVIEWER_TWO_ADVERSARIAL_REPORT.md`
- `docs/research/PUBLICATION_READINESS.md` (updated at commit `9dce4d3`)
- `docs/research/02b_PRIOR_ART_COMPLETION.md`
- `docs/research/10_LICENSING_RISK_REVIEW.md`
- `docs/research/11_HISTORY_CLEANUP_PLAN.md`
- `feature/researchos-orchestrator:docs/researchos/IMPLEMENTATION_REPORT.md`
