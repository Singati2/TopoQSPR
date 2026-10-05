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

## Literature cross-check (before implementation)

Each methodological recommendation was checked against primary sources (DOIs resolved on Crossref,
2026-10-02). All are legitimate. Six need a stated caveat, which is now part of the plan.

| Recommendation | Primary source | Legit? | Caveat added to the plan |
|---|---|---|---|
| Nested CV; tuning and model selection inside the fold | Varma & Simon 2006 (doi:10.1186/1471-2105-7-91); Cawley & Talbot 2010 | yes | — |
| Repeated CV, paired fold-level comparison against baselines | Ash et al. 2025, *JCIM* 65:9398 (doi:10.1021/acs.jcim.5c01609), the current small-molecule method-comparison protocol | yes; our 5×5 matches | Ash et al. stress *practical* significance: pre-specify a minimal meaningful ΔR² (proposed 0.05, to confirm with co-authors) instead of relying on p-values |
| CIs from the corrected repeated k-fold t-test | Nadeau & Bengio 2003 (doi:10.1023/A:1024068626366); Bouckaert & Frank 2004 | yes, approximate | Bates, Hastie & Tibshirani 2023, *JASA* (doi:10.1080/01621459.2023.2197686): naive CV intervals undercover, and CV estimates average error over training sets. Report the corrected t-test as primary and a bootstrap of the whole nested procedure as sensitivity |
| Benjamini–Hochberg over the confirmatory family | Benjamini & Hochberg 1995 | yes | BH assumes independence/PRDS; endpoints share drugs and folds. Report Benjamini–Yekutieli 2001 (doi:10.1214/aos/1013699998) as sensitivity |
| Scaffold-grouped outer CV | Bemis & Murcko 1996 (doi:10.1021/jm9602928) | yes | **Verified pitfall:** every acyclic molecule gets the empty scaffold (`''` for metformin, valproic acid, pregabalin), so they would form one giant group. Each acyclic molecule gets its own group; report the number of scaffold groups and the largest group's share |
| ΔR²(TI − size) vs a size baseline | standard confounder adjustment | yes | The size baseline must be fitted with the **same model families and tuning** as the TI model; a linear size model vs a tuned forest on TIs is not a fair comparison |
| RDKit 2D baseline minus target-defining descriptors | — | yes | Exclude whole derivation families, not single names (e.g. for Crippen logP: MolLogP, MolMR, SlogP_VSA*, SMR_VSA*) |
| Cross-check indices against ChemGraphX | Jacob et al. 2026 | yes | ChemGraphX is a *consistency* check (same senior author's group), not ground truth; ground truth stays closed forms, brute force and Mordred |
| Pilot, freeze, amendments log, GO/NO-GO | standard preregistration practice | yes | — |
| Replication of published studies | — | yes | Their properties are often ACD/Labs predictions and their index conventions may differ; with n ≈ 20 nested-CV estimates are very noisy, so CIs are mandatory and "not reproduced" must be distinguished from "too noisy to tell" |

**Verdict:** the recommendations are sound and are adopted, with the caveats above. They omit two
audit findings that the plan keeps: the composition baseline and the use of experimental rather than
computed targets for headline claims.
