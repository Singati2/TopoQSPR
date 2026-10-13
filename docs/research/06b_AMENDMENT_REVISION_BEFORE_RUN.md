# 06b — Revision of AMD-2026-10-13-E3 before any real-target fit

**Status:** committed before any model in E3 was fitted to real target values. Only the
smoke test (fully permuted targets, numbers unread) has run. This revision supersedes 06
where they differ. **Everything in E3 is exploratory.**

## R1. Wording correction: "target-blind" (owner's point 4)

The inclusion rule is **not** fully target-blind. `expdata.get_target` applies the package's
pre-specified replicate-agreement filter, which drops a parent whose replicate measurements
disagree by more than 10 °C (MP) or 0.5 log units (logP). That uses the measured values'
**spread**, not their level and not any relation to features. The correct description is:
"selection independent of target level and of any feature–target relation; a pre-specified
replicate-agreement filter (frozen in the package before this study) removes 299 of 19,799 MP
parents". **Why the filter is appropriate:** replicate disagreement beyond 10 °C mostly
signals identity errors (polymorph or salt confusion, decomposition reported as melting) or
transcription errors. Such parents contribute label noise of a size comparable to the effects
under study. The filter does not condition on the feature–target relation, so it cannot
favour one feature set. Census impact: 769 formula groups (6748 molecules) before the filter,
758 groups (6587 molecules) after. The 11 groups lost are listed in
`out/isomer_census_formula_groups.csv`. **Sensitivity (S-noflt):** the primary contrast is
also run with disagreeing parents retained (median value), seed 0 only.

## R2. Feature sets added (information-content attribution)

The target-blind information-content analysis (`out/family_information_content.csv`)
showed three things. The 11 distance indices alone distinguish every pair that any unlabeled
invariant can (33,592 / 35,369). The degree and reverse-degree families equal the 9
degree-pair edge counts plus Δ (30,902). All 104 indices add nothing beyond distance in
discrimination. Added sets: `distance_composition` (asinh of the 11 distance indices plus
composition) and `edgetypes_composition` (9 edge counts, Δ and composition). Added secondary
contrasts: S6 `distance_composition − composition`, S7 `edgetypes_composition − composition`.
These attribute any TI gain to a family. **Distinguishing more isomers is not evidence of
better prediction**; only held-out paired contrasts count (owner's point 8).

## R3. External validation, source-disjoint (owner's point 5)

Bradley contains 9,687 PHYSPROP-sourced records (8,811 of 19,500 kept parents have ≥ 1).
OPERA MP *is* PHYSPROP, so a naive Bradley → OPERA test is partly circular. E3-ext is:

* **Training:** Bradley parents with **no PHYSPROP-sourced record**, in formula groups with
  ≥ 5 constitutions *within that subset*.
* **External test:** OPERA `mp_opera` parents (package rules) whose constitution does not occur
  in the training set and whose formula is a training formula (known-formula external test,
  the E3a estimand on new molecules from a different measurement source). A second test set uses
  OPERA parents with formulas absent from training (E3b-like).
* **Fitting:** each set is fitted once on the full training set, with in-fold family and
  hyperparameter selection by grouped inner 3-fold CV (groups = constitution) using
  `evaluate.select_and_refit`, then predicts the test set.
* **Uncertainty:** paired ΔMAE and ΔR² on the same test molecules, 95% percentile intervals from
  a formula-cluster bootstrap (B = 2000, `default_rng(2026)`). Contrasts S1–S7. Secondary.
* **Caveat retained:** PHYSPROP and Bradley's other sources (Alfa Aesar catalogue,
  Karthikeyan, academic sites) may cite common primary literature. Source-disjoint is
  stronger than molecule-disjoint, but it is not proof of independent measurements.

## R4. Independent units and uncertainty

* Units: the standardised parent. Constitutions (InChIKey first block) group stereo variants;
  formula groups are the cluster unit for E3b and the bootstrap.
* Internal CV (E3a, E3b): corrected repeated-CV t over outer folds. For E3b, with folds of whole
  formula groups, the interval is approximate (a heuristic for correlated folds). It is reported
  as such, alongside the number of formula groups.

## R5. Multiplicity

The primary test (MP E3a, TIa_composition − composition) is a single test at α = 0.05. All
secondary contrasts across all cells actually run (internal and external, S1–S7) form one BH
family at q = 0.05. Negative control and sensitivity are reported, not tested.

## R6. Failure analysis (descriptive, declared now)

For MP E3a: (i) the 1,777 within-formula pairs that every unlabeled invariant must confuse
(isomorphic graphs). Report the paired absolute error difference of `TIa_composition` and
`ECFP` on those molecules versus the rest. (ii) The 10 formula groups with the largest
`TIa_composition` MAE. (iii) Molecules with absolute error > 100 °C. No exclusion follows
from any of these.

## R7. Run order and budget

MP E3a → MP E3b → E3-ext → negative control → S-noflt → logP E3a → logP E3b. Seeds rule as in 06
(0–4 if the smoke projection ≤ 10 h, else 0–2). Cells beyond the budget are reported NOT RUN.
