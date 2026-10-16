# Melting-point external gap: all indices vs edge-type counts

On the source-disjoint external test (Bradley PHYSPROP-free → OPERA PHYSPROP), unseen formulas,
TIa_composition reached R² ≈ 0.471 and edgetypes_composition ≈ 0.423. Earlier reports said the
ten degree-pair features carry the whole topology gain. This task asks what the existing archives
can and cannot say about that gap, without training anything.

Prose in this file is untrusted context. Only the block below is read by ResearchOS.

```researchos-task
{
 "title": "MP external unseen-formula gap: TIa_composition vs edgetypes_composition",
 "question": "Why does TIa_composition apparently outperform edgetypes_composition on the source-disjoint unseen-formula external test (R2 0.471 vs 0.423)? What can and cannot be inferred without new training?",
 "claims": ["C-EDGE10"],
 "steps": [
  {"id": "ingest", "adapter": "ingest_artifacts", "scope": "read-only",
   "description": "hash the archived external summaries against SHA256SUMS and their commit",
   "params": {"dir": "docs/research/out/E3/external", "files": ["summary.csv", "paired.csv"],
              "sums": "SHA256SUMS", "expect_commit": "4a8cb1d"}},
  {"id": "points", "adapter": "verify_external_points", "scope": "read-only", "deps": ["ingest"],
   "description": "re-derive every archived contrast and compare with the independent recomputation",
   "params": {"summary": "docs/research/out/E3/external/summary.csv",
              "paired": "docs/research/out/E3/external/paired.csv",
              "recomputed": "docs/research/out/E3/recomputed_external.csv"}},
  {"id": "learners", "adapter": "learner_family_audit", "scope": "read-only", "deps": ["ingest"],
   "description": "which learner each set used, in every cell where both sets exist",
   "params": {"set_a": "TIa_composition", "set_b": "edgetypes_composition",
              "external_summary": "docs/research/out/E3/external/summary.csv",
              "internal_summaries": ["docs/research/out/E3/mp_bradley_E3a/summary.csv",
                                     "docs/research/out/E3/mp_bradley_E3b/summary.csv",
                                     "docs/research/out/E3/logp_opera_E3a/summary.csv",
                                     "docs/research/out/E3/logp_opera_E3b/summary.csv"]}},
  {"id": "classify", "adapter": "classify_claim", "scope": "read-only", "deps": ["learners", "points"],
   "params": {"claim_id": "C-EDGE10", "from_step": "learners", "margin": 0.01,
              "narrowed_wording": "Within the cross-validated cells, edge-type counts plus maximum degree match the full index set (MP E3a 0.663 vs 0.657; E3b 0.630 vs 0.629). On unseen external formulas the full set scored higher (0.471 vs 0.423) with different selected learners and no interval on the difference; whether the non-degree indices add information there is unresolved."}},
  {"id": "review", "adapter": "adversarial_review", "scope": "read-only", "deps": ["learners"],
   "params": {"from_step": "learners"}},
  {"id": "fold_contrast", "adapter": "fold_level_contrast", "scope": "analysis-only", "deps": ["learners"],
   "description": "corrected-t interval for TIa_composition - edgetypes_composition from archived MP E3b fold metrics (outside Git; needs the owner to approve and to add the run archive to allow_roots)",
   "params": {"folds": "~/topoqspr_runs/research_2026-10-13_E3/cells/mp_bradley_E3b/folds.csv",
              "set_a": "TIa_composition", "set_b": "edgetypes_composition"}},
  {"id": "fixed_learner", "adapter": "fixed_learner_experiment", "scope": "compute", "deps": ["review"],
   "description": "proposed experiment; requires_owner_approval; never launched by v0.1",
   "params": {
    "protocol_inputs": ["docs/research/scripts/e3_isomer.py", "docs/research/06b_AMENDMENT_REVISION_BEFORE_RUN.md"],
    "protocol_draft": "DRAFT (not frozen, not approved). Fixed-learner comparison on the E3 external design.\n1. Data, splits and features exactly as E3-ext at protocol commit fc69883 (Bradley PHYSPROP-free training; OPERA PHYSPROP MP test; known/unseen formula subsets).\n2. Sets: composition, TIa_composition, edgetypes_composition.\n3. For EACH learner family in {ridge, elasticnet, pls, knn, svr, extratrees} separately: tune that family's default grid by the same grouped inner CV, refit, predict the test set. No selection across families.\n4. Estimand per family: dR2(TIa_composition - edgetypes_composition) on unseen formulas, formula-cluster bootstrap (2000 resamples, seed fixed before running). Report all six families.\n5. Decision (owner to confirm the margin before freezing): 'edge types carry the gain' only if the 90% CI lies inside +/-0.05 for the families that are competitive (R2 within 0.05 of the best); otherwise report the gap as unresolved or as information in the non-degree indices.\n6. Optional second arm (more costly): MP E3b repeated nested CV with the family fixed, seeds 0-2.\n7. Labelled post hoc: the question arose after seeing E3 results.",
    "resource_estimate": "Grounded in E3 archive timings: the external stage (10 sets, all six families, 14 workers) took 42 min wall on 2026-10-14, so arm 1 (3 sets) is roughly 15-25 min at 14 workers. Arm 2: E3b took 62,780 s for 10 sets (shared CPU), so 2-3 sets is roughly 3-5 h. During E3 the coordinating process used about 1.1 GB RSS and system memory use stayed below about 17 GB of 62 GB (observed 2026-10-14). Requires the machine to be free of the owner's jobs."
   }}
 ]
}
```
