# Colab notebooks: reproduce and check the TopoQSPR research results

Each notebook:
- runs top to bottom in a fresh Google Colab CPU runtime;
- checks out an exact commit, installs a pinned dependency set, and rebuilds everything from public,
  licensed sources;
- calls the same scripts as a Linux/HPC run, via `python -I docs/research/scripts/...`.

Any `FAIL` stops the notebook, and its numbers must not be used. The science lives in tested modules
(`docs/research/scripts/*.py`, `tests/test_repro_colab.py`). The notebooks are thin wrappers, generated
by `make_notebooks.py`.

**Code commit checked out by every notebook:** `104899afec01c366b9c31d11ba311aadbdd432c0` (set as `GIT_SHA` in the configuration cell). Branch:
`repro/colab-first`.

## Map: question → notebook

| question / claim | notebook | runtime (2 vCPU) | inputs (provenance) | outputs (`OUT=…/topoqspr_repro/out`) |
|---|---|---|---|---|
| Are the environment, sources, published files and OPERA restore intact? | [00](https://colab.research.google.com/github/Singati2/TopoQSPR/blob/repro/colab-first/notebooks/colab/00_environment_and_data_audit.ipynb) | ≈ 5–10 min, < 2 GB RAM | Bradley xlsx (figshare 1503990, CC0, sha256 pinned); OPERA_Data.zip (GitHub kmansouri/OPERA @22f4db6, MIT, sha256 pinned); `docs/research/data/E3` (SHA256SUMS) | `00_audit/checks_*.csv/json` |
| Do the committed E3 / P0 numbers follow from the published predictions? | [01](https://colab.research.google.com/github/Singati2/TopoQSPR/blob/repro/colab-first/notebooks/colab/01_E3_P0_reproduction.ipynb) | ≈ 5–10 min | as 00 | `01_reproduction/checks_*.csv`, recomputed contrasts |
| Fixed-learner, four-arm comparison (doc 14 r2) | [02](https://colab.research.google.com/github/Singati2/TopoQSPR/blob/repro/colab-first/notebooks/colab/02_fixed_learner_four_arm.ipynb) | smoke ≈ 30–45 min (incl. feature build); full: doc 18 §5 | as 00; features rebuilt and hash-checked | `fixed_learner/{data,units_*,results_*}`, `feature_manifest.json` |
| Is the patent MP set usable for a confirmation? (doc 15) | [03](https://colab.research.google.com/github/Singati2/TopoQSPR/blob/repro/colab-first/notebooks/colab/03_external_replication.ipynb) | metadata < 1 min; census ≈ 1–2 h (approval needed) | figshare 2007426 file 3575841 (CC0, md5 pinned) | `03_patent/figshare_metadata.json`, `census.json` |
| Statistical robustness: bootstrap unit, corrected-t coverage, training seeds, both within-formula metrics | [04](https://colab.research.google.com/github/Singati2/TopoQSPR/blob/repro/colab-first/notebooks/colab/04_statistical_robustness.ipynb) | smoke ≈ 15–30 min | as 02 | `04_robustness/{a,b,c,d}_*.csv` |
| Website branch: smoke test and public-claims check | [05](https://colab.research.google.com/github/Singati2/TopoQSPR/blob/repro/colab-first/notebooks/colab/05_website_smoke_test.ipynb) | ≈ 5 min | branch `fix/website-claims` at `SITE_SHA` | `05_website/checks_*.csv` |

## Smoke vs full

- **`smoke` / `bench`** (notebooks 02 and 04) train on **permuted labels**, or use too few replicates
  to mean anything. They test the code path and timing. **Their numbers are never results.**
- **`full`** is the frozen protocol. Notebook 02 and the training seeds in notebook 04 refuse to start
  unless the configuration cell holds `APPROVAL_SHA`, the commit that froze doc 14 (owner approval).
- The **patent census** (notebook 03) needs `APPROVE_DOWNLOAD=True`, which is owner approval for a
  638 MB download.

## Long runs, interruptions, resume

- Set `USE_DRIVE=True`. Then the cache, feature data and every finished unit live on your Drive, under
  `MyDrive/topoqspr_repro`.
- Each unit (context, learner, arm, seed, fold) is written atomically, with a sha256 record. Running the
  notebook again skips only units whose record verifies. A half-written unit is recomputed.
- **To split the work** across sessions, accounts or machines, set `SHARD="i/N"`. Each value of i
  covers a fixed subset of units. `MAX_UNITS` stops a session cleanly after that many units.
- **To finish elsewhere,** copy the `fixed_learner` folder (from Drive, or downloaded) to the other
  machine and run the same command there.
- `fixed_learner.py status OUT --mode full` lists what is done.

## Local / HPC equivalents (same code path)

```bash
git checkout 104899afec01c366b9c31d11ba311aadbdd432c0 && pip install <pins below> && pip install --no-deps -e .
python -I docs/research/scripts/restore_opera_values.py docs/research/data/E3 R
python -I docs/research/scripts/audit_environment.py docs/research/data/E3 R OUT/00_audit
python -I docs/research/scripts/reproduce_e3.py R OUT/01_reproduction
python -I docs/research/scripts/fixed_learner.py build OUT/fixed_learner
python -I docs/research/scripts/fixed_learner.py run OUT/fixed_learner --mode smoke --n-jobs 8
python -I docs/research/scripts/fixed_learner.py analyse OUT/fixed_learner --mode smoke --restored R
python -I docs/research/scripts/patent_census.py OUT/03_patent                 # metadata only
python -I docs/research/scripts/robustness.py R OUT/fixed_learner OUT/04_robustness --mode smoke --n-jobs 8
python -I docs/research/scripts/website_claims_check.py SITE OUT/05_website
```

**Pins:**
- `rdkit==2026.3.6 chembl_structure_pipeline==1.2.4 numpy==2.2.6 pandas==2.3.3 scipy==1.15.3
  scikit-learn==1.7.2 openpyxl==3.1.5 joblib==1.6.0 xgboost==3.2.0`.
- The E3 results were produced with Python 3.10.12. The notebooks were tested with Python 3.12.15,
  Colab's current major version.
- If Colab moves to a Python version for which a pinned wheel does not exist, `pip` fails in the
  setup cell and the notebook stops. Do not loosen the pins silently. Report it; a validated
  alternative set has to be checked against notebook 01 first.

## Verification status

| notebook | local (clean Python 3.12.15 + pinned set, nbclient) | Google Colab |
|---|---|---|
| 00–05 | **all six OK** at `104899a`; also from an empty work folder (00), and by resuming a persisted one (02) | **unverified**: not executed on a real Colab runtime |

Details, timings and the acceptance checklist: `docs/research/18_COLAB_FIRST_COUNTERCHECK_AND_VERIFICATION.md` §6.
