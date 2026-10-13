# 10 — Copyright and licensing risk review of the public repository (2026-10-13)

The repository `Singati2/TopoQSPR` is **public** (verified with `gh repo view`, 2026-10-13).
This review is read-only. **No change has been made.** Every action below needs the owner's
approval. I am not a lawyer; this is a technical risk assessment, not legal advice.

## Findings

| # | Item on public `main` (`a0fd2a6`) | Licence status | Risk | Evidence |
|---|---|---|---|---|
| L1 | `related_papers/Arockiaraj_2023_SARQSAR_hybrid_GHZ_descriptors_photos.pdf` (3.0 MB; 29 photos of pp. 570–585: text, Tables 1–6, Figures 1–5) | Article: Arockiaraj et al., *SAR QSAR Environ. Res.* 34(7):569–589 (2023), doi:10.1080/1062936X.2023.2239149, publisher Informa (Taylor & Francis). Crossref lists **no open licence**. The folder's own README says "Publisher PDFs must not be committed unless their license (e.g., CC BY) permits public redistribution" | **HIGH**: substantial reproduction of a copyrighted, apparently non-OA article, publicly downloadable. A co-author being an author of the article does not by itself grant redistribution rights (rights are usually assigned to the publisher; T&F author rights typically allow sharing the accepted manuscript after an embargo, not publisher pages) | `git show origin/main:related_papers/README.md`; added in `ec19cb1` (2026-10-03) |
| L2 | The same file in **git history** (every clone, every fork, GitHub's cached views) | as L1 | HIGH even after deletion from `HEAD`, because history keeps it | `git log -- related_papers` |
| L3 | ChEMBL-derived content: ChEMBL IDs in `results/expC/*/target_failures.csv`, run records and a test fixture. No ChEMBL structures or activity values in `results/` on `main` | ChEMBL: CC BY-SA 3.0. Redistribution allowed with attribution and share-alike | **LOW**, but **attribution is missing**: the repo has only an MIT `LICENSE` for code | `git ls-tree`, grep for CHEMBL IDs |
| L4 | `data/replication/*/` (10 studies: drugs, extracted properties, index values) | 7 sources CC BY 4.0. **3 sources CC BY-NC-ND 4.0** (ahmed2026, supriya2026, tang2026), quoted from each paper's JATS `<permissions>` in `SOURCE.md` | **LOW–MEDIUM**. Each `SOURCE.md` provides attribution. Extracted factual values are generally not protected expression, but for ND sources a reproduced table could be argued to be an adaptation, and NC restricts commercial reuse of anything that is protected | `data/replication/*/SOURCE.md` |
| L5 | Experimental endpoints (Bradley CC0, AqSolDB CC0, OPERA/PHYSPROP, AZ via ChEMBL) | Not redistributed on `main` as per-molecule values; downloaded at run time with pinned hashes | LOW. OPERA's data terms are unresolved (01b), so keep them un-redistributed | `src/topoqspr/expdata.py` |
| L6 | New branches `audit-2026-10-12`, `research-2026-10-13` | Contain scripts, docs, small summary CSVs (fold-level metrics, census counts, formula-group statistics). No per-molecule measured values, no structures of external datasets, no PDFs | LOW | pre-push scan, 2026-10-13 |

## Recommended correction plan (ordered; each step needs approval)

1. **Immediate, non-destructive:** make the repository **private** until L1/L2 are resolved,
   *or* at least delete `related_papers/*.pdf` from `main` in a normal commit and keep the
   bibliographic README. A normal commit stops new downloads from `HEAD` but leaves the file
   in history, so step 2 is still needed.
2. **History purge (destructive; explicit approval required):** `git filter-repo --path
   related_papers/Arockiaraj_2023_SARQSAR_hybrid_GHZ_descriptors_photos.pdf --invert-paths` on
   a fresh mirror clone. Then force-push **all** branches and tags that contain it (`main`,
   `v1.0-dev`, `audit-tests-20261012`, `audit-2026-10-12`, `research-2026-10-13`,
   `plan-freeze-v1`). Ask GitHub Support to purge cached views and the PR refs (PR #1, PR #2).
   Consequences: every commit ID after 2026-10-03 changes *again*. The run records' commit
   IDs need a second SHA map, just as after the 2026-10-12 trailer rewrite (audit F0). Do this
   **once**, before any release tag or Zenodo DOI, and publish the map.
3. **Attribution files (normal commit):** add `DATA_LICENSES.md` (or `NOTICE`) stating that
   the code is MIT, ChEMBL-derived identifiers and results are CC BY-SA 3.0 with the required
   sentence ("ChEMBL data is from http://www.ebi.ac.uk/chembl - the version of ChEMBL is
   chembl_37") plus the Mendez et al. 2019 citation, the replication extracts are under each
   source's licence as recorded in `SOURCE.md` (3 are CC BY-NC-ND 4.0), and Bradley/AqSolDB are CC0.
4. **Replication extracts (judgement call):** for the three CC BY-NC-ND sources, either keep the
   extracted values with attribution plus an explicit NC-ND notice, or replace them with a
   fetch-on-demand script plus DOIs. The latter is the conservative choice.
5. **Process:** add a pre-commit or CI check that rejects `*.pdf` under `related_papers/`, and a
   short CONTRIBUTING note: "bibliographic references only".

## What I did not do

No deletion, no history rewrite, no force-push, no visibility change, no contact with GitHub
Support. The T&F author-rights policy for this specific article was not checked against the
publication agreement. The authors should confirm what rights they retained.
