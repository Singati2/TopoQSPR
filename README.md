# TopoQSPR

Disease-anchored QSPR modelling with graph-theoretic topological indices.
TopoQSPR turns a disease concept into a reproducible QSPR study: a controlled,
versioned drug manifest from ChEMBL (confirmed MeSH IDs, clinical-phase rule,
parent/salt normalisation), 104 topological indices each traced to its primary source
and unit-tested (`src/topoqspr/`), and leakage-free repeated nested cross-validation
against null, size, composition and RDKit-descriptor baselines.

> Formerly developed as *MolecuLens* (renamed to avoid a clash with an existing product of that name).

**Status: v1.0 in development on branch `v1.0-dev`.** `main` still holds the v0.x
prototype (repository root: `main.py`, `ml.py`, `fetcher.py`, `website.py`, ...), which
has the known issues listed below and must not be used for results. The v1.0 package in
`src/topoqspr/` fixes them; see `docs/AUDIT_TOPOQSPR.md` and `docs/PAPER_ROADMAP.md`.

## Run locally

```bash
pip install -r requirements.txt
python main.py --disease "asthma"          # command line, writes outputs/
uvicorn website:app --port 8000            # web app at http://localhost:8000
```

## Quick start (v1.0, in development)

```bash
pip install -r requirements-lock.txt && pip install --no-deps -e .
topoqspr headings "migraine"                          # find the MeSH descriptor ID
topoqspr cohort --mesh D008881 --out cohort/          # controlled cohort manifest
topoqspr indices --smiles-csv my.csv --out idx.csv    # verified indices (id,smiles columns)
topoqspr run --config examples/pilot_migraine.json --out results/pilot
topoqspr rerun results/pilot/run_record.json --out results/pilot-rerun   # byte-identical check
```

```bash
docker build -t topoqspr \
    --build-arg GIT_COMMIT="$(git rev-parse HEAD)" \
    --build-arg GIT_DIRTY="$(git status --porcelain --untracked-files=no | grep -q . && echo 1 || echo 0)" .
```

gives the same CLI from the lock file (the commit and the dirty flag of the working tree are
baked in so that run records from the image name the code; do not hard-code `GIT_DIRTY=0`;
docs/PIPELINE.md §7). The pipeline, the
study-config schema, the outputs and the run record are described in
[docs/PIPELINE.md](docs/PIPELINE.md).

## Repository layout

| Path | Purpose |
|---|---|
| `config.py` | All settings: drug limits, CV folds, index list, targets |
| `disease_validator.py` | Input sanitising and fuzzy disease-name matching |
| `fetcher.py` | ChEMBL / PubChem / UniChem retrieval with caching and retries |
| `ml.py` | Topological indices, correlation, model evaluation, SHAP, drug-likeness filters |
| `main.py` | Pipeline orchestrator and Excel export |
| `website.py`, `index.html` | FastAPI service and web client |
| `paper/` | SoftwareX manuscript draft (Overleaf-ready) and the ChEMBL indication audit script |

## Known issues of the v0.x prototype (fixed in the v1.0 package)

- **Drug-set definition.** Disease matching is a substring match on MeSH headings with no clinical-phase filter, so most retrieved molecules are investigational (e.g. *diabetes*: 580 molecules, of which 103 reached phase 4 for a matched indication, and only 71 of those are small molecules, including non-antidiabetic drugs attached via label text mining; ChEMBL 37) and unrelated headings are included (*Diabetes Insipidus*). Plain-language names such as *breast cancer* do not match the MeSH term *Breast Neoplasms*. Truncation to `MAX_DRUGS` takes an unordered set, so the selected drugs can differ between runs. Salts and parents are not merged. See `paper/audit_chembl_indications.py`.
- **Evaluation.** Hyperparameters are tuned on the full data before cross-validation (not truly nested). No null or molecular-size baseline. Duplicate and size-determined targets.
- **Indices.** Hosoya Z is incorrect for molecules with more than 18 bonds. Index values are not yet unit-tested against reference values.
- **Annotations.** The `ADMET_*` columns are RDKit descriptors and rule-of-thumb flags (GI absorption, BBB), not validated predictions; the former `ADMET_Caco2_approx` heuristic was removed. Bioactivity columns take the first reported value regardless of target or unit.
- The interface previously advertised 11 models; 10 are implemented, and it now says 10.
- Also found by the audit: four indices with wrong numerators, duplicate indices (DS1 = DS2), Balaban J wrong for salts, failed indices filled with 0, a shared model instance reused across targets. Details and fixes: `docs/AUDIT_TOPOQSPR.md`, `docs/INDEX_VERIFICATION.md`.

## Citation

Manuscript in preparation for SoftwareX (see `paper/` and the [paper roadmap](docs/PAPER_ROADMAP.md)).

## License

Code: MIT; see [LICENSE](LICENSE). The original author of the v0.x code (Advik Natarajan) has granted full flexibility on licensing (`docs/ADVIK_RECOMMENDATIONS.md`).

Data: drug sets and properties retrieved from ChEMBL are licensed CC BY-SA 3.0 by EMBL-EBI. Redistributed manifests or result files that contain ChEMBL-derived data must carry attribution and the same licence.

## Required Reading

See [`related_papers/README.md`](related_papers/README.md) for the papers every contributor should read, starting with Arockiaraj, Greeni & Abul Kalaam (2023), *Int. J. Quantum Chem.* 123:e27136 ([doi:10.1002/qua.27136](https://doi.org/10.1002/qua.27136)).
