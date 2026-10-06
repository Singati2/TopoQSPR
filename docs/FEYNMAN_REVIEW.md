# Feynman review of the TopoQSPR article (2026-10-06)

Method: every claim the article makes is restated in plain words, as to a sharp outsider; then the
question "what did we actually measure that makes this true, and what would make it false?" is
answered from the result files. Claims that did not survive were changed in the manuscript the same
day; the changes are listed with each claim. Reviewer: the lead author (no agent; the agent usage
limit lifts 2026-10-07, when an independent pass can repeat this).

## Claims that survived as written

| Claim (plain words) | What makes it true | What would make it false | Verdict |
|---|---|---|---|
| Knowing only how many atoms and bonds a drug has already predicts its properties somewhat; adding 102 graph indices helps for log P and not clearly elsewhere. | `confirmatory.csv`: ΔR² vs the 2-feature size baseline on identical folds, corrected repeated-CV t, BH over the 4-member family; log P +0.15 [0.08, 0.23], q = 0.002; the other three not significant after correction. | A size baseline that is too weak (then anything beats it) — but the baseline is the same two features the genre never adjusts for, and composition (11 counts) beats the indices anyway. | holds |
| Counting the elements in a molecule predicts better than all the indices. | Pooled and in-domain R² of composition > TI on all 4 endpoints, same learners and folds. | A paired test with an interval; not run (every paired test in the frozen run is against size). | holds as a descriptive statement; the text now says so explicitly |
| Published single-index studies reproduce on their own data. | `results/expD`: Godlin's 9 fits and 504 printed index values reproduce to 3 decimals; Almulla and Kabeer reproduce where the printed graph is the resolved structure (16/20 and 15/25 compounds), and where not, every index differs at once. | A discrepancy that could not be traced to a different graph — none found. | holds; the abstract said "reproduce exactly", now "reproduce (one of them to three decimals)" |
| Under honest cross-validation those studies' indices do no better than molecular size. | Nested CV of the genre's own procedure vs the same procedure on size: TI ≤ size in 17/28 fits; the paired interval excludes zero in 0/28. | A fit with an interval above zero — there is none, but n ≤ 25 makes every interval wide, so "no better" is "not distinguishable". | holds with the softer wording now used ("statistically indistinguishable from", not "collapse to") |
| The v0.x protocol (tune on everything, report the best family) inflates R². | `results/expC`: training fit up to 0.91 vs nested 0.31 (log P, random forest); ordinary CV only +0.01 to +0.05 above nested. | — | holds; the training-fit number is the headline, the CV optimism is honestly small and said so |
| A free-text disease query retrieves 4–16× more structures than the controlled rules. | `results/expA`: counts and Jaccard per disease. | — | holds |
| Every printed value we checked from seven other groups' papers reproduces, or the discrepancy is in the paper. | `tests/test_literature_values.py`, 199 cases, each discrepancy classified against the article text. | A discrepancy we attributed to the paper that is really ours — the Rao Theorem 2.2 and 2.4 cases were re-derived by hand and checked against the Europe PMC text. | holds |
| The ChemGraphX live service and TopoQSPR give the same numbers on what both compute. | `tools/compare_sites.py` on the same graph: 21/21 indices, 13/13 entropies. | — | holds |

## Claims that were wrong or overstated, and what changed

1. **"Three cardiac glycosides (digoxin, digitoxin, acetyldigitoxin) broke the linear models."**
   Re-checked with `tools/check_extrapolation.py` (scale-free rule: a residual over 10 SD of the
   target). Four molecules, not three, and one name was wrong: lanatoside C (not acetyldigitoxin)
   for log P, digoxin and digitoxin for solubility, and atracurium (not a glycoside) for melting
   point. Each is the largest or second-largest molecule of its endpoint's set. The manuscript now
   names them from the generated macros; the amendment in `EXPERIMENT_PLAN.md` carries a dated
   correction.
2. **The mechanism was asserted, not shown.** "Unbounded indices such as Hosoya's Z make linear
   models fail" was plausible but unmeasured. Measured: on the training folds that exclude them,
   the flagged molecules' Z lies 12 to 1,593 standard deviations from every other molecule of the
   same endpoint (Szeged and edge-Wiener 12–17 SD). Note the subtlety the first check missed: on the
   whole pooled set these molecules are *not* extreme (|z| ≤ 3) because a still larger molecule sets
   the scale; they are extreme within their endpoint's set, which is what a training fold sees.
   The manuscript sentence is now macro-driven from the check.
3. **"Which endpoint passes depends on the fold design."** This compared a multiplicity-corrected
   q (random folds) with an uncorrected p (scaffold folds). Rewritten: scaffold folds change the
   precision of each estimate (melting point tighter, log P wider); the scaffold run is outside the
   confirmatory family and uncorrected.
4. **"The rest [of the recomputed R²] lie above the published value, as expected when a paper's
   graph is smaller."** No reason to expect that direction. Rewritten: they differ because the
   paper's graphs differ from the resolved structures for those compounds.
5. **"A disease's drug set is defined as much by the retrieval rule as by the disease."** The
   comparison is between two rules of our own (free-text vs controlled); the genre's hand-made
   lists are a third rule we did not measure. Rewritten in both places.
6. **"Evaluates ten regressor families."** The package offers ten; the frozen run used six. Now
   "up to ten".
7. **TI + composition vs composition** was stated as if tested; it is a difference of pooled values
   without an interval. Now labelled descriptive, with the reason the paired test was not run.

## Things an outsider would still ask, and the honest answers

* *Why 4 endpoints?* They are the measured endpoints with ≥ 200 drugs in the licence-checked sources
  after full-InChIKey matching (plan §5). More would be computed values, which the article argues
  against using as targets.
* *Is 1,821 drugs a lot?* For this genre (median 20) yes; for QSPR in general no. Scaffold-grouped
  folds show the estimates are not inflated by near-duplicates.
* *Could the indices matter for other properties?* Unknown; the article claims only the four
  tested, and says so in the limitations.
* *Why trust the applicability-domain explanation and not just drop the molecules?* The plan
  pre-specified all-molecule metrics as primary and in-domain as secondary; both are reported, and
  the explosion is reported as a result rather than removed.
* *What would change the headline?* A sixth family or a different fold seed could move log P's ΔR²
  by a few hundredths; it cannot move composition below the indices, and it cannot make the
  replicated studies beat size, because those intervals are centred near zero with n ≤ 25.
