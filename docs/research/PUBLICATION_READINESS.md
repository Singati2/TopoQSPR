# Publication readiness: TopoQSPR research track (2026-10-14; final E3 numbers 2026-10-16)

Supersedes `docs/audit-2026-10-12/PUBLICATION_READINESS.md` for the scientific claims. That
document still governs the software and provenance items. Details: 07 (results), 08 (novelty and
decision), 09 (adversarial review), 10/11 (licensing and history).

## Verdict: **CONDITIONAL GO** (evaluation-led research article; positive for MP, null for logP)

Not ready to submit today. The science for one endpoint is done and verified. The blockers are
scope (one positive endpoint), overlap with the companion manuscript, wording, and repository licensing.

## Supported claims (with evidence)

| Claim | Evidence | Label |
|---|---|---|
| Unlabeled topology adds MP information beyond element composition within known formulas: ΔR² +0.085 [0.065, 0.106], MAE −4.8 °C, 15/15 folds | E3a primary test, independent recomputation exact | exploratory amendment, primary test passed |
| …and across held-out formulas: +0.106 [0.087, 0.125] | E3b | exploratory |
| …and on a source-disjoint external set: +0.072 [0.019, 0.138] (known formula), +0.240 [0.203, 0.281] (unseen formula; composition extrapolates poorly) | E3-ext | exploratory, single fit, source-disjoint |
| No spurious gain when within-formula signal is destroyed: −0.028 [−0.052, −0.004] | negative control | control, not tested |
| Robust to the replicate filter: +0.088 [0.052, 0.123] | S-noflt | sensitivity |
| Within the cross-validated cells, 10 numbers (9 degree-pair counts + Δ) match the full index set: edgetypes_composition 0.663 vs TIa_composition 0.657 (MP E3a), 0.630 vs 0.629 (E3b); also on logP and external known formulas. On unseen external formulas the full set scored higher (0.471 vs 0.423) with different selected learners (extra trees vs SVR) and no interval on the difference, so whether the non-degree indices add information there is unresolved | edgetypes_composition vs TIa_composition, 6 contexts; ResearchOS T-a5cea8c5230c | exploratory; external gap learner-confounded |
| Better isomer discrimination does not mean better prediction: distance indices discriminate most and predict worst | E3a/E3b/ext + census | exploratory |
| No practically meaningful value beyond RDKit 2D (90% CI inside ±0.05). Zero or negative beyond ECFP (E3b −0.015 [−0.029, −0.002]) | E3 secondary | BH final (42 contrasts; E3b ECFP q 0.048) |
| logP (OPERA, E3): TI adds nothing practically meaningful beyond composition: E3a +0.016 [0.001, 0.032], q 0.056; E3b −0.001 [−0.023, +0.020]; 90% CIs inside ±0.05 | E3 logP | exploratory; equivalence shown |
| Drug cohort (E2): no detectable gain beyond composition for logP/logD/logS; MP underpowered (n = 711) and consistent with E3 | audit §7 | negative / inconclusive |
| All 104 indices verified (implementation 104/104; 44 against primary sources); registry rank 38 | 03, 03a, 03b | verification |

## Unsupported claims (must not appear)

"TIs improve prediction" (not beyond standard descriptors) · any new index or new theorem ·
degree-index linearity or rank as new (Rada 2026; companion manuscript) · "independent external
validation" · "preregistered" · "95% isomer discrimination → predictive superiority" · any positive
claim for logP/logD/logS · "the gain is carried by 10 numbers" (unqualified) · disease-specific predictive validity · "deployed v1.0 web tool".

## Answer to the owner's question

*What genuinely new understanding does TopoQSPR establish?* No new mathematics. One controlled
empirical result: on measured melting points, unlabeled topology carries real information beyond
composition. Within the cross-validated cells it is matched by the degree-pair edge counts alone (whether the other indices add anything on unseen external formulas is unresolved, see the table), and none of it survives once standard
labelled descriptors are present. A second, linked result: the index family that discriminates
isomers best is the worst predictor. A third, scoping result: the value is endpoint-dependent.
For logP, which is close to additive in atoms, it is practically zero. Everything else (linearity, rank, discrimination ordering,
collision floors, "2D descriptors win") is known and must be cited (08 §2).

## Conditions before submission (owner's decisions marked ★)

1. Done: E3 complete, BH family final. Frame the paper as endpoint-dependent (MP positive, logP null).
2. ★ Coordinate with `topological-index-orthogonality` (same authors): merge, or cite and narrow.
3. Rewrite the abstract and claims per the tables above.
4. ★ Repository: decide on the PDF history (11, Options A/B/C). Add DATA_LICENSES/ChEMBL
   attribution. Settle OPERA terms. Release no per-molecule OPERA data.
5. Recommended: one confirmatory replication of N1/N2 frozen before fitting (second measured
   endpoint, or untouched MP source).
6. ★ Nothing is merged into `main` without approval. All of this work is on `research-2026-10-13`.

## Venue

Journal of Cheminformatics (Research article). Fallback: Digital Discovery. Not SoftwareX or a
mathematics journal.
