# Audit findings, 2026-10-12

Scope: branch `main` = `v1.0-dev` at `fecb800` (history rewritten 2026-10-12 to drop
co-author trailers; see F0), frozen tag `plan-freeze-v1` → commit `38eda80` (12:43:26,
2026-10-05), archived runs under `results/exp{A,B,C,D}` and the full run outputs in
`~/topoqspr_runs/` (predictions, features, molecules; not in git).
Audit work lives on the local branch `audit-2026-10-12`; nothing was pushed.

Statuses: **CONFIRMED BUG** (reproduced, with a counterexample), **METHODOLOGICAL RISK**
(design choice that a reviewer can attack; not a coding error), **UNVERIFIED CLAIM**
(stated somewhere, not checkable here or not yet checked), **VERIFIED CORRECT**.

Severity: S1 changes a headline claim, S2 changes a secondary claim or a reported
number, S3 hygiene.

---

## F0. Provenance: archived run records cite commit IDs that no longer exist on GitHub — CONFIRMED BUG, S1

* **What.** On 2026-10-12 the history was rewritten (git filter-repo, message-only change)
  and force-pushed. Every commit ID changed; file trees did not. The archived run records
  (`results/exp*/**/run_record.json`, `run_info.json`) and the paper macro
  `\runCommit{659fac38b35f}` (`paper/results/provenance.tex:1`) cite the *old* IDs. They
  no longer resolve on GitHub. The frozen tag `plan-freeze-v1` was re-created, now
  `f2ed9e5` → commit `38eda80`.
* **Counterexample.** `git cat-file -t 659fac38` works only in this local clone (the
  object has not been garbage-collected). On a fresh clone it fails.
* **Impact.** "Each number can be regenerated from the run record" (paper §5) is false for
  an outside reader until the mapping is published.
* **Fix (done, local).** `docs/audit-2026-10-12/sha_map_2026-10-12.tsv`: all 236
  old→new pairs from `.git/filter-repo/commit-map`. The nine IDs cited by archived
  records map as follows: `659fac3→61f63d8`, `2c5e1d1→325ecb2`, `2e56cfd→2c0294c`,
  `591dd18→611cbbb`, `5da2789→5ef35d7`, `963a362→70f3648`, `b59fa7e→ab2c30f`,
  `d353b46→aafb25a`, `fe92837→1a504a9`. Trees are byte-identical (`git diff` old↔new HEAD
  empty).
* **Still to do (author decision).** Cite the new ID plus the map in the paper. Do not
  edit the frozen run records. A Zenodo archive of the release is the durable fix.
* **Lesson.** Any later history rewrite (e.g. purging `related_papers/`) breaks the IDs
  again. Do it before the release tag and the Zenodo DOI, and only once.

## F1. Catastrophic out-of-fold predictions come from two unbounded descriptors, not one — CONFIRMED (mechanism); paper misattributes half of it — CONFIRMED BUG (manuscript), S2

* **Reproduction.** `docs/audit-2026-10-12/scripts/catastrophic_predictions.py`
  → `out/catastrophic_main_confirmatory.csv`. 8 selected-learner folds have R² < −1. All
  of them were fitted by **elastic net**. In each, one or two molecules carry 69–97% of the
  fold's squared error:

  | cell | molecules | y_true | y_pred | driving feature (robust z vs training fold) |
  |---|---|---|---|---|
  | logS, TI+composition (r2 f0) | digoxin, digitoxin | −4.08, −5.29 | −1406, −940 | Hosoya **Z** (3×10⁷) |
  | logS, rdkit2d (r2 f0; r1 f4) | digoxin, digitoxin | −4.08, −5.29 | −1594, −1070; +157, +104 | RDKit **Ipc** (6×10⁷) |
  | logP, TI+composition (4 folds) | lanatoside C | 0.07 | 39–58 | Hosoya **Z** (10¹⁰) |
  | MP, TI+composition | atracurium | 189.5 | −1215 | Hosoya **Z** (5×10⁹) |

* **Where the paper goes wrong.** `paper/main.tex:202` explains the two pooled R² < −1
  cells of the ladder (logS TI+composition −502 and logS rdkit2d −654) only through Hosoya's
  Z. The rdkit2d cell contains no TI. Its failure comes from RDKit's `Ipc` (Bonchev–Trinajstić
  information content, not log-scaled), which `features.rdkit2d_names`
  (`src/topoqspr/features.py:254`) keeps for every target.
* **Mechanism (verified).** Digoxin and digitoxin are near-analogues. When both fall in
  the same test fold, the training fold has nothing within orders of magnitude on Z or Ipc.
  The kNN applicability domain flags them (`in_domain = 0`), but flags are not used in the
  headline.
* **Fix.** Correct the sentence. Treat Ipc like Z: name it, and report the asinh
  sensitivity analysis (§B2 of the research plan, run below). Keep the frozen numbers.

## F2. Corrected t-test reports p = 0 and a zero-width CI when all fold differences are equal — CONFIRMED BUG, S3 (no archived result affected)

* `src/topoqspr/stats.py:66` (`corrected_ttest`). Before the fix, 25 identical differences
  of 0.1 returned `p = 0.0, CI = [0.1, 0.1]`. The existing test
  (`tests/test_stats.py::test_corrected_ttest_zero_variance`) asserted exactly this.
* **Why it is wrong.** With s² = 0 the corrected SE is not estimated, it is degenerate. A
  zero-width interval claims certainty that 25 correlated folds cannot give.
* **Archive check.** Of the paired groups in all four Experiment B runs (140 + 140 + 16 + 32),
  every zero-variance case is a reference set compared with itself (mean 0, p = 1). The
  p = 0 branch was **never reached**, so no frozen result changes.
* **Fix (local commit).** If s² = 0 and the mean is non-zero, t, p and the CI are NaN, with
  `reason = "zero variance ... not estimable"`. The all-zero case keeps p = 1 and CI [0, 0].
  A new regression test covers this, and the old assertion was changed.

## F3. Abstract: "On 1821 approved drugs with measured end points" — CONFIRMED BUG (manuscript), S2

`paper/main.tex:36` (macro `\nMolecules` = 1821). Only **959** of the 1821 standardised
parents have any measured endpoint. Per endpoint the count is 711 (MP), 505 (logP),
369 (logD) and 351 (logS) (`targets.csv` of the frozen run). Fix: "1821 approved drugs, 351–711
with a measured value per endpoint".

## F4. "preregistered" is used for a plan that was frozen in git but never registered — CONFIRMED (wording), S1 for reviewer trust

`paper/main.tex:38` ("leakage-free, preregistered statistics"), `:265` ("with
preregistered statistics"), `paper/supplement.tex:63` ("Preregistration and amendments"),
`docs/PAPER_VS_CHEMGRAPHX.md:15`, `docs/CHEMGRAPHX_COMPARISON.md:62`. The OSF
registration is still a `\todo` (`main.tex:305`), and that todo says registration happens
"before the main analyses". The main analyses ran on 2026-10-05.

A git tag carries the author's own clock and was re-created on 2026-10-12. It is evidence
of intent, not third-party evidence of timing. Fix: "pre-specified in a plan frozen as git
tag `plan-freeze-v1` before the main runs (not externally registered)". Registering on OSF
*now* is possible only as a registration after data collection, labelled so.

## F5. Scaffold-vs-random sentence: "log P's is the same size but wider" — CONFIRMED BUG (manuscript), S2

`paper/main.tex:235`. Random split (frozen, 102 indices): ΔR²(TI−size) = +0.153
[+0.076, +0.230]. Scaffold split (104 indices): +0.098 [−0.050, +0.246]. That is about a third
smaller, and the scaffold interval includes zero. The two runs also differ in registry. The
amended SVR run shows that the registry alone moves the random-split value by 0.003
(+0.150). So the drop is mostly the split, but the paper never compares the same registry
under both splits. The audit's sensitivity run does (RESULT_RECONCILIATION.md §3).

## F6. Experiment D: "statistically indistinguishable from the size baseline" is non-significance, not equivalence — METHODOLOGICAL RISK (framing), S2

`paper/main.tex:36` and the Experiment D text. TI_study − size, nested CV, 95% corrected-t
interval widths in R² units (`results/expD/*/nested.csv`):

| study | n | fits | median CI width | fits inside ±0.2 |
|---|---|---|---|---|
| Godlin 2026 | 21 | 27 | 0.24–0.27 | 13 |
| Kabeer 2026 | 25 | 24 | 0.69–0.91 | 1 |
| Almulla 2026 | 20 | 27 | 1.9–2.1 | 0 |

Only Godlin's data come close to supporting "no better than size". Almulla's and Kabeer's
are uninformative. Fix: "under nested CV none of the 78 single-study comparisons shows a
gain over size, and only the Godlin set is large enough to bound the gain (|ΔR²| ≤ 0.2 in
13 of 27 fits)".

## F7. The deployed site is still the v0.x prototype — CONFIRMED, S2

`https://moleculens.onrender.com/` serves `index.html` titled "MolecuLens — Drug QSPR
Pipeline" and the API "Drug QSPR Pipeline API" (`/openapi.json`). That is v0.x
`website.py`. Its `/health` returns `{"version":"1.0.0"}` because the string is
hard-coded in `website.py:28,56`, which is misleading. The v1.0 app
(`topoqspr.web.app`, title "TopoQSPR") is not deployed there. Checked 2026-10-12 ~16:40
EDT. The paper (`main.tex:156,168`) describes and screenshots the v1.0 interface; its own
`\todo` already says "once deployed". Fix: deploy `render.yaml` (needs the Render account),
and change v0.x `/health` to report `0.x`.

## F8. README is stale — CONFIRMED, S3

`README.md:11-22` says "`main` still holds the v0.x prototype" and lists `python main.py`
/ `uvicorn website:app` first. `main` and `v1.0-dev` have been identical since 2026-10-06.
Fix: v1.0 quick start first; v0.x as "legacy, do not use for results".

## F9. Representation collisions exist in the data and are large — VERIFIED (new evidence)

`scripts/representation_collisions.py` → `out/collision_*.csv`. The 102 frozen indices
(plus the 11 composition counts) cannot distinguish these approved-drug pairs. Each pair
has *identical* TI+composition vectors and different measured values:

| pair | relation | MP (°C, Bradley) |
|---|---|---|
| ephedrine / pseudoephedrine | diastereomers | 34 / 117.5 |
| codeine / hydrocodone | constitutional isomers (allylic alcohol vs ketone) | 155 / 198 |
| betamethasone / dexamethasone | C16 epimers | 232.5 / 262.5 |
| sulfamethazine / sulfisomidine | constitutional (ring N positions) | 198.5 / 243 |
| atenolol / practolol | constitutional (amide reversed) | 146.5 / 135 |
| epirubicin / doxorubicin | C4′ epimers | 344.53 / 204.5 (*see F10*) |

Ethanol (`CCO`) and dimethyl ether (`COC`): all 104 indices are identical
(`out/ethanol_dimethyl_ether.csv`). NIST WebBook gives normal boiling points of 351.5 ± 0.2 K
and 248.2–249.5 K. The implied empirical R² ceilings are **not binding** on these datasets,
however (TI alone: 0.94 MP, 0.98–0.997 elsewhere; TI+composition ≥ 0.995). The gap
between TI models (R² 0.2–0.4) and these ceilings is a learning and generalisation gap, not an
information wall. See RESULT_RECONCILIATION.md §4.

## F10. Possible source error: epirubicin melting point 344.53 °C — UNVERIFIED CLAIM (data)

The Bradley dataset row 16461 (source "DrugBank", DB00445) gives 344.53 °C, against
204.5 °C for its epimer doxorubicin. The two-decimal value suggests a predicted or
converted value. Not checked against a primary source. The row is in the frozen MP set. The
frozen data stay as they are; the row belongs in the data card's list of suspect values.

## F11. RDKit 2D baseline for experimental logP keeps Crippen `MolLogP` — METHODOLOGICAL RISK, S2

`src/topoqspr/target_providers.py:290-293` keeps alternative estimators of the target for
measured endpoints. The decision is documented, and a practitioner would make it too. But
Wildman–Crippen atom contributions were fitted to measured logP values, and those
overlap PHYSPROP, the source of OPERA logP. So rdkit2d's 0.84 on logP is partly a
recall of a fitted model. Report it as "conventional baseline incl. a fitted logP model",
or add rdkit2d minus Crippen as a sensitivity set.

## F12. Bootstrap and permutation controls use a different learner than the confirmatory estimand — METHODOLOGICAL RISK (documented), S2

Frozen: control learner ridge; amended: SVR. The confirmatory statistic uses in-fold
family selection (`selected`). The controls test "TI with ridge/SVR", not the selected
procedure, so they *support* but do not *validate* the confirmatory interval. The
amendment says so; the abstract should not imply more.

The resampling mechanics were verified correct (`pipeline._one_resample`): copies keep
their original CV group, connectivity groups resample as units, and the original-data row
reproduces the main analysis. Permutation resolution: p_min = 1/201 = 0.005 at K = 200,
enough for BH over 4 endpoints but not finer. The corrected t-test treats the 25 folds with
the Nadeau–Bengio n_test/n_train inflation; this is a heuristic correction for correlated
folds, not an exact interval (Bates et al. 2023), and the paper cites it as such.

## F13. Median imputation inside model pipelines vs "failures are never imputed" — METHODOLOGICAL RISK (documentation), S3

`features.py` docstring: "Failures are never imputed" (true at the feature level: NaN
plus a reason). `models.py:11` then median-imputes NaN inside every pipeline, fitted on the
training fold. This is leakage-safe but contradicts the shorthand. Fix: "never imputed in
the feature table; imputed with the training-fold median inside models".

## F14. Dataset licence attributions — UNVERIFIED CLAIM, S2

`src/topoqspr/expdata.py:219,231` record OPERA and Biogen data as "MIT (repository
LICENSE)". The code repository's MIT licence does not necessarily cover the data. OPERA's
curated PHYSPROP-derived sets have their own provenance (EPA). Bradley (CC0) and AqSolDB
(CC0, Harvard Dataverse) look right. ChEMBL (CC BY-SA 3.0) is right. Needs a check against
each dataset's own terms before redistributing processed values.

## F15. Archived CSVs re-read without exact float parsing — CONFIRMED, S3

pandas' default `read_csv` float parser is not round-trip exact. It changes 44,290 of the
archived `features.csv` cells in the last bit, which moved a re-fitted elastic-net R² by
4×10⁻⁶ (the first attempt of the E1 run caught this with its reproduction assertion). The
archived files are fine; `_write_csv` writes full `repr` precision. Readers must pass
`float_precision="round_trip"`. `tools/check_extrapolation.py:30` does not; it computes
z-scores, so the effect there is immaterial. Any future "re-analyse from archived CSV"
tool should use exact parsing.

ChemGraphX overlap (from LITERATURE_2026.md, verified in `paper/refs.bib`): Micheal
Arockiaraj is a co-author of both ChemGraphX (J Cheminform 18:39, 2026) and TopoQSPR.
This must be disclosed to the editor, and the cross-check of 41 shared indices is not
fully independent in authorship (it is independent in code).

---

## Verified correct (selected)

* Test suite at `fecb800`: **1103 passed, 3 skipped, 8 deselected** (92 s, venv Python 3.10).
  After the F2 fix: test_stats 23/23. The full re-run is in `TEST_LOG.md`.
* Archived Experiment B summaries reproduce exactly from the archived feature matrix
  (asserted to 1e-9 in `scripts/sensitivity_composition_asinh.py` for all raw sets).
* Connectivity grouping: 67 parents in 33 InChIKey-first-block groups. The random split
  uses grouped folds, so stereoisomers never straddle train and test. Scaffold groups give
  members their first member's scaffold. Bootstrap copies keep their group.
* Censoring and replicate rules (`expdata.py`): ChEMBL `standard_relation` ≠ "=" and assay
  limits become `censored_value`. Bradley `donotuse` rows are excluded. Parents whose
  replicate SD exceeds 0.5 log units or 10 °C are dropped.
* Manuscript build: `main.tex` 34 pp, `supplement.tex` 14 pp, 0 LaTeX errors, 0 undefined
  citations or references, 1 BibTeX warning (`lundberg2017`: empty pages). 26 + 1 `\todo`s,
  all author items.
* 102 vs 104 indices: the frozen confirmatory, scaffold-free ladder and Experiment C use
  102. `main_confirmatory_protocol_svr`, `main_scaffold` and `exploratory_families` use
  104. Experiment D uses 104. The paper states this (`main.tex:200`) but compares 102-random
  with 104-scaffold in prose (F5).
