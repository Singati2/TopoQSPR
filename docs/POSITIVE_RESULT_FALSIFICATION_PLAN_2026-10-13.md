# TopoQSPR — Adversarial plan to test genuinely positive graph-descriptor contributions

**Date:** 2026-10-13  
**Status:** PROPOSAL / NOT APPROVED, NOT RUN, NOT A FINDING  
**Audience:** Claude Code, Codex, co-authors, independent reviewers  
**Entry point:** read `CLAUDE.md`, `docs/INDEPENDENT_ADVERSARIAL_AUDIT_2026-10-12.md`, `docs/EXPERIMENT_PLAN.md` (all amendments), `docs/CHECKPOINT.md`, and the October 13 local audit if/when it has been pushed.

> **Mandatory instruction to Claude:** Do not agree with this document by default. **Countercheck every recommendation, hypothesis, dataset suggestion, mathematical claim, novelty claim, metric, and proposed experimental design.** Look for a counterexample or an experiment that would refute it. Mark each recommendation **VALID / VALID WITH CHANGES / INVALID / UNRESOLVED**, with exact supporting evidence. A proposed path to a *positive* result is not permission to engineer significance or suppress negative evidence.

## 0. Background — evidence that must remain intact

The existing four-endpoint selected-model confirmatory experiment on approved drugs (351–711 measured compounds per endpoint) showed only a logP improvement over the two-count *size* baseline after BH adjustment. It did **not** show superiority over element composition. The October 13 audit report supplied by the project team describes an exploratory E2 similarity-cluster analysis comparing **topological indices + 11-count element composition** against composition:

| Endpoint | E2 paired ΔR² [95% interval] | BH q |
|---|---:|---:|
| Melting point | +0.069 [−0.017, +0.154] | 0.45 |
| logP | +0.017 [−0.080, +0.115] | 0.84 |
| logD | −0.045 [−0.128, +0.039] | 0.57 |
| Solubility | +0.011 [−0.104, +0.126] | 0.84 |

**These E2 figures come from the user-supplied local audit log, not independently retrieved remote results.** At the time this proposal was authored, the local audit branch and its 10 commits had not been pushed to GitHub `main`; verify whether that has changed. Do not silently rewrite the frozen confirmatory experiment, conflate the 102-index confirmatory registry with the later 104-index software registry, or call a nonsignificant result evidence of equivalence. Preserve full CIs, corrected testing, all negative results, Ridge extrapolation failures and amended SVR *family-specific* controls. No experiment below has been executed or independently validated as part of this document.

## 1. REQUIRED FIRST STEP — falsify or qualify this plan before any expensive computation

Create `docs/research/RECOMMENDATION_COUNTERCHECK_2026-10-13.md` (or an equivalently named versioned report) containing:

| ID | Proposed recommendation | Is mathematical premise true? | Does existing literature already do it? | Is data valid/available/licensed? | Minimal counterexample or falsification test | Decision (VALID / VALID WITH CHANGES / INVALID / UNRESOLVED) | Evidence and actionable correction |
|---|---|---|---|---|---|---|---|

**Check at least all of the following, rather than assuming them:**

1. Is "unlabeled TI fails to distinguish isomers" universally true? **No.** Many constitutional isomers have *different unlabeled graphs* and can have different classical TI values. The narrower correct assertion is that **if their unlabeled graphs are isomorphic, all exact unlabeled-graph invariants coincide**, even when atom labels differ. Construct both a supportive example and a counterexample to an overbroad claim.
2. Example: ethanol (heavy-atom path C–C–O) and dimethyl ether (C–O–C) both have molecular formula C2H6O and the same **unlabeled** three-vertex path P3, but different vertex-label placements and properties. Independently verify chemistry, identity, standard conditions, graph representation (heavy atom vs explicit H), and a reliable property source. Show that labeling can distinguish them, **not** that it guarantees better prediction.
3. **Collision reduction is not prediction improvement.** A descriptor can distinguish more molecules yet predict worse because of noise, sample size, correlation, or model regularization. Construct at least one synthetic example illustrating this difference.
4. An 11-count elemental-composition vector is constant within a fixed molecular formula, but different formula groups can have different baselines; test what happens when a formula occurs only in test folds. Do not use target-derived formula means for held-out formulas.
5. **Literature novelty must be investigated:** atom-/bond-weighted graph invariants, heteroatom-sensitive Zagreb/Wiener/Randić variants, functional-group descriptors, edge labeling, graph kernels, hashed fingerprints, chemical graph neural networks, and isomer-discrimination studies already exist. Check 2026 prior art. No claim of inventing chemical labeling, isomer studies, or a new TI family without precise novelty evidence.
6. **QM9 suitability and limitations:** the primary publication documents ~133,885 small CHONF molecules and **6,095 C7H10O2 constitutional isomers** with computed quantum-chemical properties. Verify exact archive, identifiers, target names, units, exclusions and data license. **QM9 is not measured melting point, aqueous solubility, experimental logP or a representative approved-drug dataset.** It can be a mechanistic/computational benchmark, not direct experimental replication of TopoQSPR's four endpoints. If this limitation invalidates a proposed comparison, reject or redesign it.
7. Check whether approved-drug cohorts actually include sufficiently large molecular-formula isomer groups with matched measured endpoints; calculate counts **before looking at target outcomes**. It may be infeasible on the current cohorts. If infeasible, report **NO-GO**; identify a suitable external *measured-property* dataset with clear source/assay/temperature/units, not a cherry-picked subset.
8. Distinguish genuine connectivity effects from stereoisomer, tautomer, salt, protonation, stereochemical representation and dataset-source artifacts. Canonicalize the correct chemical identity and keep legitimate distinctions based on the scientific estimand.
9. Check leakage from RDKit descriptors or learned representation preprocessing, similarity clustering, feature selection, scaling, hyperparameter tuning and duplicate molecules; perform all fitted operations within training folds only.
10. Verify whether proposed sample-efficiency gains are plausible and identifiable. A favorable small-n curve is a *hypothesis*, not a promise.
11. **Countercheck this entire plan itself** against the current pushed code, the newer local audit, existing experiment registry and original paper; state what is already done or redundant.

If important premises fail, **stop and report with recommendations to replace or drop that experiment**. Do not execute a four-hour run just because it is listed here.

## 2. Proposed E3 — molecular-formula-controlled isomer benchmark (exploratory)

**Scientific question:** At fixed elemental composition, does connectivity-sensitive information improve held-out prediction over composition alone, and how does unlabeled topology compare with chemically labeled descriptors and standard fingerprints?

**Datasets (conditional, not pre-approved):**
- Feasibility audit of TopoQSPR's current four measured-endpoint sources; table per endpoint of identical-formula groups with ≥2, ≥3, ≥5 distinct standardized constitutional structures, unique source records and properties, missingness, intra-formula target variance, and cluster sizes. Establish a **target-blind** inclusion rule.
- QM9 C7H10O2 subset as a **separate computed-property** diagnostic if data rights and preparation are acceptable. Cite Ramakrishnan et al., *Scientific Data* (2014), https://doi.org/10.1038/sdata.2014.22 and data citation https://doi.org/10.6084/m9.figshare.978904. Choose only actual QM9 target columns; do not mislabel DFT calculations as experimental assays.
- Identify a genuinely measured external property dataset if the drug endpoint sample is inadequate; log source version, exact target definitions, conditions and permissions.

**Baselines and candidate descriptors — identical molecules, splits, model-search budget and outer folds:**
1. Null/target mean and formula-only 11-element composition; molecular size separately.
2. Existing **unlabeled** TI registry, with feature-version and non-finite-value policies preserved.
3. Formula + existing unlabeled TI.
4. Chemically **labeled** indices (only after E4's correctness/novelty audit), plus formula.
5. Standard atom/bond-aware RDKit 2D and Morgan/ECFP descriptors, alone and alongside formula.
6. A complexity-matched compact baseline (for fair attribution of gains to chemistry vs descriptor count).
7. Optional GNN/3D baselines **only** if data, resources and novelty question warrant them; name any unavailable baselines.

**Evaluation:**
- Primary experimental design to be **locked in a dated amendment before observing performance** after the feasibility check. Separate (a) held-out constitutional isomers within formulas represented in training, (b) unseen formula groups for transfer, and (c) dissimilarity-cluster / scaffold transfer where meaningful. These test different estimands; a good result on one cannot be substituted for another.
- Ensure one chemical entity or near-duplicate cannot occur across training/test; avoid group leakage through isomer pairs generated from the same molecule. Group/cluster split **before** any target-driven selection or fitted preprocessing.
- For same-formula contrasts, formula-only features are constant *within a formula*; report MAE/RMSE and target variation alongside R² (which can be unstable/undefined when evaluation target variance is zero). Pairwise target differences are dependent if pairs share a molecule: do not treat all pairs as independent observations.
- Define a primary metric/direction, primary endpoint or a multiplicity-adjusted family, smallest meaningful effect and an analysis plan *before* the confirmatory follow-up; run seed/fold integrity checks. Use cluster-aware uncertainty at a suitable independent unit (molecule, formula group or chemical series, depending on split) and report the actual independent sample count. Explain if classical corrected repeated-CV t intervals are not justified for the proposed grouping.
- Report per-formula sample sizes, value ranges, descriptor collisions, negative outcomes and failures. When formula-held-out testing gives an unseen formula, **never** calculate a target-derived formula mean using held-out labels.

**Falsification:** The claim that chemical labeling adds predictive value is false in a tested setting if properly held-out improvements are absent or underperform stronger baselines; report that without selectively changing the dataset, split threshold, endpoint or descriptor definitions.

## 3. Proposed E4 — chemically labeled topological descriptors, with proofs and counterexamples (exploratory)

**Hypothesis, NOT established novelty:** Some atom-/bond-aware, chemically meaningful graph features convey property-relevant information that an unlabeled graph invariant provably cannot encode in particular collision classes.

Before implementing, do a systematic comparison with prior atom-weighted/heteroatom-sensitive graph index literature and RDKit/ECFP. Define exactly what is new **relative to the closest mathematical and cheminformatics construction**. If a proposed formula is a renaming or special case of an existing one, **mark NOT NOVEL and do not claim invention**.

Minimum independent correctness gates:
- Precisely define graph type: heavy atoms vs explicit H; elements; formal charges; bond order/aromaticity; isotope, stereochemistry and tautomer handling; missing and disconnected graphs.
- Give explicit equations, graph-isomorphism invariance proofs **for chemically labeled graph isomorphisms**, and note any counterexamples to injectivity; no descriptor is automatically complete.
- Work small graphs by hand: P2, P3 (C–C–O vs C–O–C), paths, cycles, branched trees, disconnected cases, aromatic heterocycles and stereo examples as appropriate. Compare to a separately coded oracle/enumeration and RDKit invariants.
- Construct two pairs: (A) **unlabeled-isomorphic but chemically distinct** with a labeling-aware descriptor separation where justified; (B) two nonisomorphic labeled chemical graphs with a **collision under the new descriptor**, unless a mathematically sound injectivity theorem and scope are actually proven. Detect numeric and algorithmic collisions separately; `np.isclose` alone is not a proof of symbolic equality.
- Check complexity and overflow, adversarial high-degree/high-size molecules, atom-renumbering invariance, salt mixtures, invalid SMILES, finite-range scaling, and train-only normalization; add reproducible tests.
- **Mathematical gain and empirical gain are separate claims.** A valid graph theorem does not imply a prediction advantage or novel index.

Deliver exact statement → proof or finite witness → independent check → unit test → measured predictive result, with limitations. Use Lean 4 if justified and available; never call a Python test a formal proof or imply Lean was run if it was not.

## 4. Proposed E5 — low-data sample-efficiency study (optional; exploratory)

Question: Do simple, validated chemical graph descriptors retain useful inductive bias when training data are scarce?

- Compare composition, unlabeled TI, labeled TI, RDKit 2D and ECFP using **predefined nested subset sizes and seeds**; sample only from each training fold, leaving the outer test fold untouched.
- Equalize learner tuning and data budgets; include constant-capacity and regularized baselines. Repeat enough times for uncertainty and report failed fits and coverage.
- Predefine how to evaluate area under a learning curve and/or performance at a specific training size. Adjust for multiplicity if several sizes/endpoints are tested.
- **Counterexample requirement:** identify at least one endpoint/data regime where extra graph features hurt; report when increased complexity *degrades* generalization.
- Reject a "small-data advantage" claim if it disappears against a fair descriptor-count- and learner-matched comparison.

## 5. Mandatory scientific safeguards

1. **Do not optimize for statistical significance.** No label peeking, cherry-picking isomer groups, multiple unreported endpoints, hidden failure exclusions, selective seeds, threshold shopping or post hoc relabeling as confirmatory.
2. Keep the original study and frozen protocol untouched. New analyses begin as **exploratory**, and any later confirmatory study needs a genuinely prospective protocol and disjoint held-out confirmatory evidence; a new Git tag after outcomes have been inspected is not prospective confirmation.
3. Compare **incremental information beyond composition**, not just beyond molecule size. If classical unlabeled TI cannot resolve a pair due to a graph-automorphism/isomorphism collision, report this formal limitation without overclaiming all isomers.
4. Baseline against competitive RDKit 2D and Morgan/ECFP fingerprints, not just null/size/composition. The claim "competitive" needs supported intervals, not a single higher point estimate.
5. Check ChEMBL, QM9 and other data licences/redistribution separately. Maintain parent-structure, source ID, target unit, version, hash and duplicate/cluster provenance.
6. Do not rewrite, drop, or hide the existing negative E1/E2 or frozen Experiment B results; discuss scientific limitations openly. Do not assert zero information from a nonsignificant difference or equivalence without a powered, margin-defined equivalence test.
7. Never say "the code is verified" simply because pytest passed. Independently crosscheck chemical identities, graph formulas, numerical results, split integrity and statistical inference.

## 6. Execution sequence / hard gates

**Gate 0 (first; no expensive runs):** Inspect current GitHub HEAD, unpushed `audit-2026-10-12` work and open PR #2. Confirm audit integration status; avoid conflicting experimental code. Write the recommendation countercheck and prior-art/source audit. Every recommendation gets a VALID / VALID WITH CHANGES / INVALID / UNRESOLVED verdict, a counterexample where applicable, source, and explicit go/no-go.

**Gate 1 (data feasibility):** Execute target-blind isomer group census, chemistry standardization audit, subset hash, duplicate leakage checks and QM9 measured-vs-computed distinction. If data are inadequate, stop and publish NO-GO.

**Gate 2 (protocol):** Write `docs/research/ISOMER_EXPERIMENT_AMENDMENT.md` with frozen hypotheses, data/splits, metrics, nulls, learners, primary comparisons, effect threshold, CIs and multiplicity. Commit it *before* the substantive run, record the commit SHA; label first research run exploratory.

**Gate 3 (minimal smoke/counterexamples):** Hand-checked graph witnesses; independent implementation; labels-permuted or target-shuffled negative control within appropriate strata; no group leakage; reproducible from a clean tree. Smoke outputs must not be used to choose a favorable confirmatory hypothesis.

**Gate 4 (main study):** Only after gates pass, run E3 then justified E4 and optional E5. Preserve exact configs, environment locks, data checksums, seeds, per-fold predictions, failures, model fits and hashes. Do not modify the original frozen experiment, results, or manuscript tables silently.

**Gate 5 (adversarial re-check):** Independently recompute the headline metric and adjusted inferential values from archived predictions, audit one formula/chemical series end-to-end, identify plausible confounding, compare against all baselines and try to reproduce the result on an independent measured-property dataset.

**Gate 6 (publication):** If supported, write carefully scoped findings, proof limits and counterexamples; otherwise publish a useful negative result. Do not call it publication-ready until provenance, deployment/licences and co-author/venue checks are resolved.

## 7. Required report back to the project owner

Produce Markdown:
- `docs/research/RECOMMENDATION_COUNTERCHECK_2026-10-13.md`: per-recommendation validation/rejection ledger with **counterexamples** and citations;
- `docs/research/ISOMER_DATA_FEASIBILITY.md`: endpoint × formula-group census, licenses, source identifiers, target-blind eligibility, NO-GO possibilities;
- `docs/research/ISOMER_EXPERIMENT_AMENDMENT.md` (only after feasibility), locked before computation;
- `docs/research/MATH_AND_NOVELTY_AUDIT.md`: claims, exact definitions/proofs, collisions, prior art, independent small-graph checks;
- `docs/research/RESULT_RECONCILIATION.md`: if a run occurs, all per-fold results, CIs, adjusted tests, hash/provenance and null controls, with confirmed vs exploratory status;
- `docs/research/PUBLICATION_DECISION.md`: defensible positive/negative claims, limitations, and what evidence would change the verdict.

**Report concise final answer:** Which recommendations survived adversarial checking? Which were **invalid** and why? What counterexamples were found? What can be done now without risking the frozen study? Which experiments did you *actually* run? Which remain unrun, infeasible, or unverified?

## References / verification starting points

- Ramakrishnan, R. et al. (2014), *Quantum chemistry structures and properties of 134 kilo molecules*, **Scientific Data** 1, 140022. https://doi.org/10.1038/sdata.2014.22 (6,095 C7H10O2 isomers and computed properties; not experimental drug endpoints).
- Original QM9 data citation: https://doi.org/10.6084/m9.figshare.978904 (verify archive version, content and reuse rights before downloading).
- Existing TopoQSPR evaluation/audit: `docs/EXPERIMENT_PLAN.md`, `docs/INDEPENDENT_ADVERSARIAL_AUDIT_2026-10-12.md`, `docs/FEYNMAN_REVIEW.md`, `docs/INDEX_VERIFICATION.md` and current frozen result manifests.

---

**Authoring instruction:** This document records *questions and proposed falsification tests*, not experimentally supported positive conclusions. The highest-priority deliverable is the **independent countercheck**, not a favorable `p`-value.
