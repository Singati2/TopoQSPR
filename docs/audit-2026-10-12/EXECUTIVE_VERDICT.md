# Executive verdict (audit 2026-10-12)

**Decision.** TopoQSPR is publishable as a *rigorous negative-leaning evaluation with
software*, not as evidence that topological indices are useful predictors.

**What survived.** On 351–711 measured approved drugs per endpoint, unlabeled-graph indices
show no detectable gain beyond an 11-count elemental composition. The audit's exploratory test
(TI+composition vs composition, same folds, fit-free asinh transform, random and scaffold
splits, frozen 102-index registry) gives ΔR² intervals that include zero for all 4 endpoints
× 2 splits. Point estimates are ≤ 0.02 except melting point (+0.05 / +0.09). The frozen
headline (TI beats size for logP, +0.15 [0.08, 0.23]) reproduces exactly but shrinks to
+0.10 [−0.05, +0.25] on scaffold folds with the *same* registry.

**Update 2026-10-13 (amendment E2, locked before running).** On a similarity-cluster split the
primary family (TI+composition − composition, 4 endpoints, BH) gives q = 0.45–0.84. **No
detectable gain beyond composition; equivalence within ±0.05 not shown** for any endpoint. ECFP
beats TI+composition for MP (+0.115 [+0.046, +0.183]). RDKit 2D beats it for MP, logP and logD
(+0.17 to +0.28). TI+composition still beats size by +0.14 to +0.43 everywhere, but that gain is
carried by composition. Supported (exploratory) message: unlabeled topological indices showed
no detectable gain beyond an 11-count element composition on these endpoints; equivalence is
not shown, so "no information" is not claimed. MP is
the only endpoint with a consistently positive (non-significant) estimate.

## Strengths (verified)

* Reproducibility: archived results regenerate to 1e-9 from the archived feature matrix. 1103
  tests pass at `fecb800` (1111 after the audit fixes). Run records and connectivity-grouped, leakage-safe nested CV are implemented as
  documented. The bootstrap keeps resampled copies in one fold.
* Honest frozen design: mandatory size and composition baselines, dated amendments, and the
  catastrophic failures are reported rather than hidden.
* A cohort-building feature (disease → MeSH → ChEMBL manifest) that the literature search did
  not find elsewhere (search scope in LITERATURE_2026.md).

## Fatal if left unfixed (all cheap to fix)

1. "Preregistered" is used for a plan never registered (F4).
2. The run records' commit IDs do not exist on GitHub after the 2026-10-12 rewrite (F0); the SHA map
   must be published.
3. Abstract counts: "1821 drugs with measured endpoints" (959 have any) (F3).

## Major

* The rdkit2d catastrophe is caused by RDKit `Ipc`, not Hosoya Z as the paper says (F1).
* The scaffold sentence misstates the logP change (F5). Experiment D non-significance is framed as
  equivalence (F6). The deployed site is still v0.x (F7). The ChemGraphX author overlap must be
  disclosed. ECFP baseline missing. The unlabeled-graph limitation needs to be stated (in-data
  collisions: ephedrine/pseudoephedrine MP 34 vs 117.5 °C with identical TI+composition).
* The frozen controls use another learner than the headline (F12).

## Publishability by venue (detail in PUBLICATION_ROADMAP.md)

| Venue | Now | After fixes |
|---|---|---|
| SoftwareX | no (length, private repo, licence file name) | yes, as a cut-down tool paper |
| J Cheminform (Research/Methodology) | no | **yes, recommended**: evaluation-led, with E1 and ECFP added |
| Digital Discovery | no | plausible fallback |

## Single next experiment (done 2026-10-13 as amendment E2; see RESULT_RECONCILIATION §7)

Planned as: E4 + E1 confirmation: declare in a dated amendment, *before running*, the
composition-conditional comparison (asinh, as run here) plus an ECFP baseline, on a
similarity-cluster split. Then run it once. The possible outcomes and their effect on the paper:

* MP gain over composition holds (interval > 0): a narrow, positive claim that topology carries
  melting-point information (packing and symmetry), consistent with chemical intuition.
* No endpoint gains: the paper's message becomes "unlabeled TIs are redundant with composition
  for these properties". That is a clear negative result and the stronger replication-style paper.
* ECFP ≫ TI+composition everywhere: TI-QSPR on drugs should be benchmarked against fingerprints
  by default, which is the recommendation to the genre.
