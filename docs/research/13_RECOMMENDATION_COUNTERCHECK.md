# 13. Countercheck of doc 12 (independent publication assessment)

**Date:** 2026-10-16.

**Checked:** `docs/research/12_…_2026-10-16.md` at commit `411e252`, against the committed E3
outputs (`docs/research/out/E3/`), the local run archive (`~/topoqspr_runs/research_2026-10-13_E3`,
sha256-checked) and reports 01–11.

**Labels:**
- `VALID`: correct as written.
- `MODIFY`: right direction, but the wording or plan needs a change.
- `REJECT`: wrong.

No item was rejected.

## 1. Factual claims in doc 12

| doc 12 claim | check | verdict |
|---|---|---|
| Review basis `9dce4d3` (research) and `5fea2fa` (ResearchOS) | `git log` | VALID |
| MP E3a +0.085 [0.065, **0.106**], MAE 4.8 °C, 15/15 | `out/E3/mp_bradley_E3a/paired.csv`: +0.0852 [0.0650, **0.1055**], MAE 4.807, 15 | **MODIFY.** The upper bound rounds to **0.105**. The slip originates in 07 and `PUBLICATION_READINESS.md` (written by the assistant) and was copied into doc 12 |
| E3b +0.106 [0.087, 0.125] | +0.1062 [0.0870, 0.1254] | VALID |
| External +0.072 [0.019, 0.138] / +0.240 [0.203, 0.281] | `out/E3/external/paired.csv` | VALID |
| Unseen-formula R²: composition 0.23, TIa+comp 0.47, ECFP 0.57, RDKit 2D 0.67 | `external/summary.csv`: 0.231, 0.471, 0.565, 0.673 | VALID |
| logP E3a +0.016 [0.001, 0.032], q 0.056; E3b −0.001 [−0.023, +0.020] | `logp_opera_E3*/paired.csv`; q from 07 §3b | VALID |
| E3b ECFP_TIa − ECFP −0.015 [−0.029, −0.002], q 0.048 | paired.csv; 07 l.50 | VALID |
| Within-formula centred R² 0.18–0.19 (TIa+comp), 0.52–0.59 (labelled); pre-declared metric negative | `within_formula_posthoc_*.csv`: 0.189 / 0.184; ECFP 0.534 / 0.516, RDKit 2D 0.588 / 0.574; 07 §5 | VALID |
| Negative control −0.028 [−0.052, −0.004]; S-noflt +0.088 [0.052, 0.123], seed 0 | `negative_control_paired.csv` (0/5 folds); `sensitivity_noflt_paired.csv` (5/5, n = 6,916) | VALID |
| MP 6,725 / 758; logP 3,994 / 495; 3 seeds × 5 folds | 07 l.9–12; `run_info.json` | VALID |
| 104 indices, rank 38, 81 reducible | 03a l.12, 116; 08 K1–K2 | VALID |
| MP chosen after E1/E2 (outcome-informed) | 06 l.7–12 ("makes the choice **outcome-informed**") | VALID |
| The "independent recomputation" used the same predictions | `e3_recompute_and_failures.py` docstring | VALID |
| Copyrighted PDF gone from branch tips, still in history | `git ls-tree` on all 7 remote branches: only `paper/figures/fig_*.pdf` (own figures); `git log --all` still finds the related_papers files | VALID |
| N2: "missing paired uncertainty" on the full-vs-edge gap | true at `411e252`; resolved by §2 below | VALID then, superseded now |

## 2. P0 step 1, done (no training)

**Method.** `docs/research/scripts/p0_full_vs_edgetypes.py`, outputs in `docs/research/out/P0/`.
- Fold R² was recomputed from predictions, and equals the archived fold R² (max |diff| = 0).
- Internal cells use the E3 corrected resampled t.
- The external test uses the E3 formula-cluster bootstrap: 2,000 resamples, seed 2026, re-checked
  with seed 99. The two seeds agree within 0.001.

| context | R² full | R² edge | Δ full − edge [95% CI] | folds full better | same family in folds |
|---|---:|---:|---|---:|---:|
| MP E3a | 0.657 | 0.663 | −0.006 [−0.015, +0.003] | 3/15 | 12/15 |
| MP E3b | 0.629 | 0.630 | −0.001 [−0.014, +0.011] | 7/15 | 10/15 |
| logP E3a | 0.828 | 0.831 | −0.003 [−0.014, +0.008] | 6/15 | 11/15 |
| logP E3b | 0.766 | 0.786 | −0.020 [−0.043, +0.003] | 2/15 | 13/15 |
| External, known formula (1,014; 211 formulas) | 0.458 | 0.484 | −0.026 [−0.070, +0.007] | — | no (ET vs SVR) |
| **External, unseen formula (7,327; 4,382 formulas)** | 0.471 | 0.423 | **+0.048 [+0.025, +0.072]**; MAE 1.5 °C [0.2, 2.8] lower | — | **no (ET vs SVR)** |

**Reading:**
1. The external unseen-formula gap **is not test-sampling noise**: the formula-clustered interval
   excludes zero.
2. It is still **learner-confounded**. In the single context where it appears, extra trees
   (full set) is compared with SVR (edge types).
3. In every context where it can be compared, the edge-type set is **equal or better**. Within
   folds where both picked the same family, the mean Δ is −0.003 to −0.004 (MP) and −0.004 to
   −0.019 (logP). That result is post hoc and descriptive.
4. **The allowed claim is unchanged:** the narrowed `PUBLICATION_READINESS` wording at `9dce4d3`
   stays correct.

**Next step.** The fixed-learner experiment (doc 14) is now the decisive test. It is not a formality.

## 3. Recommendations

| doc 12 item | verdict | evidence / change |
|---|---|---|
| Executive: CONDITIONAL GO, narrow empirical article, not ready | VALID | agrees with 08 and PUBLICATION_READINESS |
| Venues: J Cheminform, then Digital Discovery; SoftwareX only for a separate software paper | VALID | 08 l.42 |
| N1 primary contribution: the design, not new descriptors | VALID | 08 |
| N3: discrimination ≠ prediction; check fixed learners first | VALID | the distance-index ranking also uses selected learners; add `distance_composition` as a secondary arm when doc 14 runs (MODIFY doc 14 if the owner wants it, at about +50% cost) |
| N4/N7: composition-relative, endpoint-dependent | VALID | 07 §3b |
| Known items (K1 linear in m_ij, rank 38, collisions) | VALID | 08 K1–K3, 02b |
| N2 unresolved | VALID | now quantified (§2) |
| Companion `topological-index-orthogonality` manuscript | VALID | owner decision; unchanged |
| P0.1 recover fold scores, external predictions, paired intervals | VALID | done (§2) |
| P0.2 selected learners and hyperparameters per fold | VALID | families in `out/P0/learner_families_full_vs_edgetypes.csv`; hyperparameters are in `data/E3/cells/*/folds.csv` (`params`) |
| P0.3 fixed learner, "at minimum SVR and extra trees" | **MODIFY** | Choosing those two because they won is itself post hoc. Run **all six** families, with ET and SVR as co-primary (pre-declared here, before fitting) and the other four reported. See doc 14 §3 |
| P0.3 the estimand is the same | **MODIFY** | A fixed-learner Δ is a **different estimand** (feature × learner) from E3's best-selected model per representation. Label it new and exploratory. It cannot revise E3 |
| P0.4 freeze before fitting, owner approval | VALID | doc 14 §9 |
| P0.5 three-way decision (equivalent / helps / unresolved) | VALID | doc 14 §6, margin δ = 0.02 proposed |
| P1 untouched measured source; target-blind census first | VALID | doc 15 |
| P1 "an independent MP collection is ideal" | **MODIFY** | The binding constraint is **no shared records with Bradley or PHYSPROP**. Bradley already aggregates PHYSPROP, Alfa Aesar, Karthikeyan, DrugBank and academic sites. The one large CC0 candidate found is the patent-mined set (Tetko 2016, figshare 2007426). Even that can only be labelled `source_disjoint_data` (doc 15 §2) |
| P1 keep RDKit 2D / ECFP at full strength | VALID | doc 15 §5 |
| P1 do not assume OPERA independence | VALID | 07 already uses "source-disjoint", never "independent" |
| P2 report the negative pre-declared and positive post hoc within-formula metrics together | VALID | 07 §5 already does; the manuscript must too (doc 16) |
| P2 absolute performance next to each ΔR²; no "non-significant = equivalent" | VALID | logP equivalence in 07 uses 90% CI inside ±0.05 (TOST-style), which is correct |
| P2 audit corrected t, cluster bootstrap, BH families | VALID | 04 audited the code. A coverage simulation for the corrected t under formula grouping has **not** been done; add it to the open list |
| P2 artifact risks (duplicates, isomorphic graphs, salts, tautomers) | VALID | 05 counts isomorphic unlabeled pairs (5.5% Bradley). Salts and parents are handled by the package pipeline. Tautomers are not separately audited (open) |
| P2 feature vs feature × learner distinction | VALID | doc 14 |
| P3 title and paper arc | VALID | doc 16 |
| P3 no parallel overlapping submissions | VALID | owner decision |
| P4 history, licences, private copy, no OPERA per-molecule values | VALID | the public dataset in `docs/research/data/E3/` follows this: Bradley (CC0) in full, OPERA values removed and restorable with hash checks |
| P4 website not presented as deployed v1.0 | VALID | branch `fix/website-claims` (`5c7d11d`, not merged) relabels the live site's claims |
| P5 ResearchOS is traceability, not evidence; do not expand it before P0/P1 | VALID | no ResearchOS work in this round |

## 4. Errors found in doc 12
1. The E3a upper CI is 0.105, not 0.106 (inherited, see §1).
2. "Missing paired uncertainty" for the N2 gap: correct when written, now superseded (§2).

No other factual errors were found.
