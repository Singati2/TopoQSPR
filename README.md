# TopoQSPR

Disease-anchored QSPR modelling with graph-theoretic topological indices.
Type a disease name; TopoQSPR retrieves the associated drugs from ChEMBL, computes
135 degree-, reverse-degree-, neighbourhood-degree-sum- and distance-based
topological indices per molecule, and evaluates regression models predicting
physicochemical properties from those indices.

> Formerly developed as *MolecuLens* (renamed to avoid a clash with an existing product of that name).

**Status: pre-release (v0.x).** The v1.0 work below is in progress and the
current outputs should not be used for publication yet.

## Run locally

```bash
pip install -r requirements.txt
python main.py --disease "asthma"          # command line, writes outputs/
uvicorn website:app --port 8000            # web app at http://localhost:8000
```

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

## Known issues being fixed for v1.0

- **Drug-set definition.** Disease matching is a substring match on MeSH headings with no clinical-phase filter, so most retrieved molecules are investigational (e.g. *diabetes*: 580 molecules, of which 103 reached phase 4 for a matched indication, and only 71 of those are small molecules, including non-antidiabetic drugs attached via label text mining; ChEMBL 37) and unrelated headings are included (*Diabetes Insipidus*). Plain-language names such as *breast cancer* do not match the MeSH term *Breast Neoplasms*. Truncation to `MAX_DRUGS` takes an unordered set, so the selected drugs can differ between runs. Salts and parents are not merged. See `paper/audit_chembl_indications.py`.
- **Evaluation.** Hyperparameters are tuned on the full data before cross-validation (not truly nested). No null or molecular-size baseline. Duplicate and size-determined targets.
- **Indices.** Hosoya Z is incorrect for molecules with more than 18 bonds. Index values are not yet unit-tested against reference values.
- **Annotations.** `ADMET_Caco2_approx` is a heuristic formula, not a model. Bioactivity columns take the first reported value regardless of target or unit.
- The interface advertises 11 models; 10 are implemented.

## Citation

Manuscript in preparation for SoftwareX (see `paper/` and the [paper roadmap](docs/PAPER_ROADMAP.md)).

## License

Code: MIT; see [LICENSE](LICENSE) (pending written agreement of the original author of the v0.x code).

Data: drug sets and properties retrieved from ChEMBL are licensed CC BY-SA 3.0 by EMBL-EBI. Redistributed manifests or result files that contain ChEMBL-derived data must carry attribution and the same licence.
