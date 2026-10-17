# 17. Verification report for docs 13–16, the P0 analysis and the public E3 dataset

**Date:** 2026-10-16/17 (America/New_York).

**Machine:** Linux 6.8.0 x86_64, 16 CPUs, 62 GB RAM, Python 3.10.12, rdkit 2026.3.6,
chembl_structure_pipeline 1.2.4, numpy 2.2.6, pandas 2.3.3, scipy 1.15.3, scikit-learn 1.7.2.

**Git:**
- Work was done in a detached worktree (`~/moleculens-verify12`) at `origin/research-2026-10-13` = `411e252`.
- Commits were pushed to `research-2026-10-13` only, without force. The SHAs are in the commit log.
- `main` was not touched (`3ef0b50`).

**Running jobs:** none of the owner's protected jobs were running. Nothing was stopped or signalled.

## 1. What was actually done

| item | how | status |
|---|---|---|
| Doc 12 factual claims | compared with `out/E3/*` CSVs and reports 01–11 (doc 13 §1) | **read and recomputed**; one rounding error found (0.106 → 0.105) |
| Rounding fix | separate commit `b06f5b5` (07 ×2, 08, PUBLICATION_READINESS) | done |
| P0 step 1 (full vs edge types) | `scripts/p0_full_vs_edgetypes.py` on the local archive. Fold R² recomputed from predictions = archived fold R² (max \|diff\| 0). External split rebuilt with `e3_isomer.build` (data building only), and the split sizes asserted against the archived 1,014 / 7,327 | **rerun from predictions, no training** |
| Interval stability | external bootstrap with seed 2026 and 99 | interval ends agree within 0.001 |
| Determinism | P0 run twice on the archive, and once on the public dataset after restore | all three outputs identical (max \|diff\| 0.0) |
| Public dataset | `scripts/publish_e3_dataset.py` → `docs/research/data/E3/` (17 files + SHA256SUMS, 18 MB) | built |
| OPERA value removal and restore | `removed_columns.json` hashes; `scripts/restore_opera_values.py` re-derives the values from the pinned OPERA zip | **verified locally** (3/3 columns match, 0 unmatched). The local run used the package's download cache. Colab downloads fresh |
| Colab notebook | `docs/research/colab/P0_full_vs_edgetypes.ipynb` (pinned versions) | **structure validated, not executed on Colab**. The steps it runs were executed locally in the same order |
| Website Colab smoke test | `fix/website-claims` `221f958`: `colab/smoke_website_baseline.py` + notebook | script run locally: PASSED. Notebook not executed on Colab |
| Doc 14 (fixed-learner protocol) | draft | **not run**; awaits owner approval |
| Doc 15 (replication feasibility) | licence of the patents MP set checked via figshare API (CC0) | census **not run**; awaits approval |
| Doc 16 (claim ledger) | every number traced to a committed file | written |

## 2. Not attempted

- **No model training of any kind:** the fixed-learner experiment and the confirmation are pending.
- **No coverage simulation** for the corrected t under formula grouping (doc 13 §3, open).
- **No tautomer audit** (open).
- **Colab not run.** The main risk there is RDKit/InChIKey differences on Colab during the OPERA
  restore. Pinned versions reduce that risk, and the hash check catches it: the notebook reports
  FAIL rather than producing wrong numbers.

## 3. File identities (sha256)

| file | sha256 |
|---|---|
| `out/P0/internal_full_vs_edgetypes.csv` | 44ba4dca30d10a0d84759ddac44867fc0141c749ac397ad52a808a7a62dd4648 |
| `out/P0/external_full_vs_edgetypes.csv` | 8e6db4c97690b6db1c8999b686cef0407d82b1888176f6c0cc92c9dbe12ec078 |
| `out/P0/learner_families_full_vs_edgetypes.csv` | c5b199f8e084a874759f14ef9772c8aed965e55fc8c692b0e99897c105ab5657 |
| `out/P0/external_split.csv` | 90ead03f38fec4b1663d97f3c25889eed0f9a1110bb946c9c6a33ae8833b2ded |
| `scripts/p0_full_vs_edgetypes.py` | bb77ce1c61aea218b1538239aa22bb152d76e4a87200ed1a9043564538530bce |
| `scripts/restore_opera_values.py` | 6fa81a25671b99f753e18f3f1aebc03be05c95e34864760b2f5786676334b305 |
| `scripts/publish_e3_dataset.py` | cfa835b687d6af23699b751a9bf5ffce52dfac152e6e81d769667565c123146d |
| dataset | see `data/E3/SHA256SUMS`; original archive files in `data/E3/removed_columns.json` |

## 4. Independence of this verification

- The same assistant wrote the E3 code, the recompute and this check. The P0 script reuses the
  E3 interval formulas on purpose, so its intervals are comparable to E3's.
- That makes this a consistency check plus a new paired contrast over the **same predictions**. It
  is labelled `consistency_check` / `independent_recomputation`, never independent replication.
- The Colab notebooks exist so the owner can re-run it on separate hardware.
