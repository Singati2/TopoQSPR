# Review of `RECOMMENDATIONS.md` (2026-10-02)

Each recommendation was checked against the code, the expert audit and the draft experiment plan.

| § | Recommendation | Verdict | Where it lands |
|---|---|---|---|
| 1 | Finish audit, debug only on a quarantined pilot, freeze the plan with SHA/date, log amendments | adopt (already planned) | `EXPERIMENT_PLAN.md` §2, header, Amendments |
| 2 | Verify every index; fix Hosoya, DS1/DS2, reverse degree, missing values; provenance table; cross-check against ChemGraphX | adopt | v1.0 index module (being built and attacked now); `INDEX_VERIFICATION.md`. ChemGraphX code is MIT-licensed Python with callable functions, so it is used as a second oracle next to Mordred |
| 3 | Deterministic, controlled disease cohorts with manifest | adopt (already planned) | v1.0 cohort module (being built now) |
| 4 | Size is the critical confounder; ΔR²(TI − size) primary with CI and ΔMAE | **adopt, with one addition** | Primary endpoint changed to ΔR²(TI − size) + ΔMAE. **Kept** the element-composition baseline as a pre-specified secondary: the audit showed all indices are blind to element type and bond order, so without it a TI "gain" over size can simply be composition. This is a baseline, not a new method, so it does not conflict with "do not add new methods" |
| 4 | Quantify how much of each TI/family is explained by size | adopt (new) | `EXPERIMENT_PLAN.md` §7 |
| 5 | Everything data-dependent inside the fold; report ordinary → nested → scaffold CV | adopt; scaffold step new | Experiment C chain; evaluation module already supports grouped folds |
| 6 | Run the panel exactly per the frozen rule; keep negative results | adopt (already planned) | §3.3 |
| 7 | Replication of 3–5 studies with pre-specified eligibility; full chain incl. ordinary and scaffold CV | adopt; eligibility made concrete | Experiment D: fixed Europe PMC query, date-ordered screen, first eligible studies, chosen with co-authors |
| 8 | Paired fold-level stats, CIs, BH over the confirmatory family, confirmatory vs exploratory, n per endpoint, collinearity | adopt; explicit split new | §7, §10 |
| 9 | SHAP secondary, outer-fold only, stability reported | adopt (already planned) | §9 |
| 10 | Tagged release, provenance, CI, pinned env, CLI, one-command reproduction, Zenodo | adopt (already planned) | Roadmap M6 |
| 11 | Narrow, complementary novelty claim | adopt (already in paper) | Paper §1 positioning paragraph |
| 12 | Let the evidence decide; negative result can be central | adopt | §1 |
| — | GO / NO-GO checkpoint report before the main experiment | adopt (new) | `EXPERIMENT_PLAN.md` §11, roadmap M6b, `docs/CHECKPOINT.md` (to be written) |
| — | No GNNs or extra descriptor families unless the core result calls for it | adopt | The landscape audit's suggestion of learned baselines is deferred |

**One point where the plan is stricter than the recommendations.** The recommendations treat the
multi-disease panel as the main analysis. The audit found that approved-drug cohorts per disease are
small (roughly 40–100 parents) and that computed targets are close to circular. The confirmatory
family therefore uses experimental endpoints on the pooled approved-drug set, plus each main-panel
disease that has at least 50 measured values. Per-disease results on computed targets are exploratory.
