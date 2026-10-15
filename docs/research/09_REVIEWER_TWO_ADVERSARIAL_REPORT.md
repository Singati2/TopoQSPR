# 09 — Reviewer 2: adversarial report on the candidate paper (2026-10-14)

This is written as a hostile but fair referee would write it, against the claims in 08 §3. Each
point gives the objection, then its severity, then the response available from our own evidence.

## Major

**R2-1. "Exploratory amendment" is a garden of forking paths.** E3 was designed after E1/E2 had
been seen. The endpoint (MP), dataset (Bradley) and baseline (composition) were chosen knowing that
MP was the endpoint with the best E2 point estimate.
*Severity: high.* *Response:* the protocol (06, 06b) was committed at `fc69883` before any
real-target fit. The primary test is single and declared. Secondary contrasts are BH-controlled.
The negative control is null. But the choice of MP was informed by E2. That makes the result
**hypothesis-generating plus internally replicated** (E3a → E3b → external), not confirmatory.
State it so. A confirmatory run on an untouched endpoint or source would remove the objection
(08 condition 5).

**R2-2. The composition baseline is a straw man.** Eleven element counts cannot represent anything
structural. Beating them shows only that "some structure helps", which nobody doubts. On unseen
external formulas composition reaches R² 0.23, barely above size (0.22), so the +0.24 is mostly
composition's failure.
*Severity: high for framing, low for validity.* *Response:* the composition reference is the
right control for the specific question "does unlabeled topology add anything to *what the
molecule is made of*?". It isolates the graph. The paper must present the ECFP and RDKit 2D
comparisons with equal weight. They show that TIs are dominated. Report absolute R² next to every Δ.

**R2-3. Learner confounding.** Each feature set picks its own learner by inner CV: mostly SVR
for TIa_composition, extra trees or SVR for composition, elastic net for ECFP_TIa (E3b) but extra
trees for ECFP. ΔR² therefore mixes
representation with learner choice.
*Severity: medium.* *Response:* this is the pre-declared estimand, best-in-class per
representation, and it is the practically relevant comparison. A fixed-learner breakdown (same
learner, two representations) is not in the archive. It needs a refit and would be a labelled
supplementary analysis.

**R2-4. Only one endpoint under the controlled design.** logP is pending. E2's drug-cohort results
for logP, logD and logS showed no detectable gain.
*Severity: high for generality.* *Response:* restrict every claim to MP until logP E3 is reported.
Present E2's nulls.

**R2-5. Overlap with the authors' own companion manuscript.** The m_ij reduction, the information
ceiling and "RDKit 2D > topology" appear in `topological-index-orthogonality`.
*Severity: high (editorial).* *Response:* cite it and disclose it in the cover letter. Confine
this paper to N1–N4 (08), or merge the two.

**R2-6. MP label quality.** Bradley mixes sources. Amino acids and zwitterions decompose. Salt or
parent mismatches occur (1.8% salts in OPERA MP). Replicate disagreement is filtered by a rule that
could be tuned.
*Severity: medium.* *Response:* the rule was pre-specified. S-noflt (filter off) gives +0.088
[0.052, 0.123]. The failure analysis names the worst formulas. A sensitivity run excluding
amino-acid-like formulas would be post hoc and must be labelled so.

## Minor

- **External uncertainty.** These are single-fit intervals. The bootstrap covers test sampling
  only. Say this explicitly. Report only an independent-bootstrap robustness check (CI ends within
  0.007), not a significance claim beyond that.
- **"Source-disjoint" ≠ independent.** OPERA has no reference field, so overlap at the reference
  level cannot be excluded. Keep the wording.
- **Seeds.** Three outer seeds (from the smoke rule) give 15 folds. The corrected t accounts for
  overlap, but intervals are wide relative to 5 seeds. Acceptable; report it.
- **Within-formula metric.** The pre-declared within-formula R² is negative for all TI sets. The
  centred metric that makes the story look better is post hoc. Report both, with the post hoc label.
- **asinh transform.** It is applied to TIs and RDKit 2D. It is
  fit-free and monotone, so it cannot leak. It can still change learner behaviour, so justify it
  with E1 (raw-TI linear-learner failures).
- **Why 104 indices?** 81 reduce to 10 numbers and the registry rank is 38. Readers will ask what
  computing the rest bought. Answer: almost nothing for MP, and that is itself a result (N2).
- **Distance-set scaling.** distance_composition's weaker performance could reflect feature
  scaling or learner fit rather than information. Partly answered: the same pipeline is used for all
  sets, and the external unseen result agrees.
- **Reproducibility.** The MP cells took about 34 h wall time on 14 workers (roughly 20 CPU-days). Provide the archived predictions and hashes,
  plus a reduced-seed reproduction script.
- **Software claims.** The live website is v0.x (audit F7). Do not describe it as the v1.0 tool.

## Verdict (as Reviewer 2)

**Major revision.** The core empirical finding (N1) is carefully supported, internally replicated
and survives a negative control. The headline must be "unlabeled topology adds real MP information
beyond composition, all of it in 10 degree-pair counts, none of it beyond standard descriptors". It
must not be "topological indices predict melting point". R2-1, R2-4 and R2-5 have to be resolved,
or explicitly scoped, before acceptance.
