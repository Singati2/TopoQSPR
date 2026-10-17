# E3 dataset (isomer study, run of 2026-10-16)

This folder holds the inputs and outputs of the E3 run that the committed results in
`docs/research/out/E3/` were computed from: inclusion lists, fold assignments with selected learners
and hyperparameters, and per-molecule predictions. It was built from the private run archive by
`docs/research/scripts/publish_e3_dataset.py`. Every file is listed in `SHA256SUMS`.

## What is included, and licences

| files | source | licence | measured values included? |
|---|---|---|---|
| `inclusion_mp_bradley.csv`, `cells/mp_bradley_E3a/`, `cells/mp_bradley_E3b/` | Jean-Claude Bradley Open Melting Point Dataset (figshare doi:10.6084/m9.figshare.1031637, file 1503990) | **CC0** | **yes** (`y_true`, °C) |
| `inclusion_logp_opera.csv`, `cells/logp_opera_E3a/`, `cells/logp_opera_E3b/` | OPERA curated PHYSPROP logP | OPERA repository MIT; upstream PHYSPROP terms unknown (`docs/research/01b_DATA_LICENCES.md`) | **no**: `y_true` removed |
| `external/` | OPERA curated PHYSPROP MP (test set); model trained on PHYSPROP-free Bradley | as above | **no**: `y` removed |

The OPERA tables keep identifiers (standardised parent InChIKeys), formulas, fold and split
assignments, and our model predictions. To get the measured values back, run
`docs/research/scripts/restore_opera_values.py`. It downloads the pinned public OPERA zip (sha256
checked), re-aggregates per parent exactly as E3 did, and verifies the restored column against the
sha256 recorded in `removed_columns.json`. `removed_columns.json` also records the sha256 of each
original archive file.

## Files

| file | columns |
|---|---|
| `run_info.json` | commit at start, versions, seeds, inclusion counts, feature-matrix hashes per representation |
| `inclusion_*.csv` | `id` (parent InChIKey), `formula`, `constitution` (InChIKey connectivity block) |
| `cells/*/folds.csv` | one row per representation × seed × outer fold: sizes, selected `family`, `params`, inner MSE, `r2`, `mae`, applicability-domain coverage |
| `cells/*/predictions*.csv.gz` | `repeat` (seed), `fold`, `row` (InChIKey), `pos`, [`y_true`], `y_pred`, `in_domain`, `set` |
| `external/predictions_no_y.csv.gz` | `id`, one prediction column per representation, `formula` |
| `external/test_split.csv` | `id`, `formula`, `test` (`known_formula` / `unseen_formula`: formula present in the training set or not) |
| `external/summary.csv`, `external/info.json` | archived per-representation R²/MAE, selected family, test sizes |

**Estimands:**
- E3a: groups = constitution, so the formula is seen in training.
- E3b: groups = formula, so whole formulas are held out.

The design is in `docs/research/06_PROSPECTIVE_STUDY_AMENDMENT.md` and `06b`; the results are in `07`.

## Attribution

- **Bradley:** Bradley J-C, Williams A, Lang A. *Jean-Claude Bradley Open Melting Point Dataset.*
  figshare (2014), doi:10.6084/m9.figshare.1031637.v2.
- **OPERA:** Mansouri K, Grulke CM, Judson RS, Williams AJ. *J Cheminform* 10:10 (2018),
  doi:10.1186/s13321-018-0263-1; Mansouri K et al. *SAR QSAR Environ Res* 27:911–937 (2016); data
  from github.com/kmansouri/OPERA (MIT), commit 22f4db60e385.

## Reproduce on Colab

`docs/research/colab/P0_full_vs_edgetypes.ipynb` covers the whole path:
1. clone the repository;
2. check `SHA256SUMS`;
3. restore the OPERA values;
4. re-run the P0 analysis;
5. compare with `docs/research/out/P0/`.
