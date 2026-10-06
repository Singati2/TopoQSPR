# The TopoQSPR website against chemgraphx.vercel.app: direct comparison

Measured 2026-10-06. TopoQSPR: `topoqspr.web.app:app` at commit of this date, served locally
(`uvicorn ... --port 8000`; the Render deployment waits for an account, see Status). ChemGraphX:
the live site https://chemgraphx.vercel.app/ and the API its client calls
(`https://vpqebjeug2dtqoi3exx6nlfdta0jnxyl.lambda-url.us-east-1.on.aws`, hard-coded in
`client/src/App.jsx` of github.com/aditya-jindal/chemgraphx; request format read from
`client/src/tabs/HomeTab.jsx` and `server/graph.py`: a `textareaData` form field holding a JSON
adjacency list, 1-based).

## 1. Same graph, both services (`tools/compare_sites.py`, reproducible)

Graph: aspirin, hydrogen-suppressed, 13 vertices, 13 edges, sent as an adjacency matrix to
TopoQSPR and as ChemGraphX's JSON adjacency list to ChemGraphX.

| | TopoQSPR | ChemGraphX |
|---|---|---|
| Requests needed for all values | 1 (`/api/indices`) | 3 (`/graph_information`, `/degree_based`, `/distance_based`) |
| Values returned | 191 (104 indices + 91 entropies, with a reason for every undefined value) | 149 (52 + 52 degree, 24 + 21 distance) |
| Wall time per request (3 repeats) | 0.10–0.11 s (local) | 0.12–0.18 s per endpoint (live Lambda, warm) |
| Shared indices under a fixed name map | 21 of 21 agree to 1e-9 (edge Wiener after ChemGraphX's documented convention, which counts adjacent edges at distance 0) | |
| Shared degree entropies | 13 of 13 agree to 1e-9 | |
| Graph information | 13 vertices, 13 edges on both; degree edge partition identical | |

So on what both compute, the two services return the same numbers; the earlier code-level
cross-check (`docs/CHEMGRAPHX_CROSSCHECK.md`: 41 shared indices, 1,763 values; 1,414 entropies)
holds on the live service too.

## 2. What each site offers (verified against the live site and the repository)

| Capability | ChemGraphX | TopoQSPR |
|---|---|---|
| Inputs | .mol, .pdb files; JSON adjacency list pasted in a textarea (the placeholder says "adjacency matrix", the server parses only the JSON list) | SMILES lists, SDF, MOL, PDB, adjacency matrix or list (validated: symmetric, 0/1, no loops, connected); batch of 50 molecules in one request |
| Disease → drug set | no | MeSH heading search, ChEMBL cohort with the full manifest and exclusion reasons shown and downloadable before anything is computed |
| Indices / entropies | 76 indices, 73 entropies | 104 indices, 91 bond-wise entropies; distance entropies not implemented |
| Formula and primary source per index | no | yes: searchable catalogue of all 195 entries (`/api/catalog`) |
| Undefined values | "indeterminate" / inf | NaN with a reason |
| QSPR modelling | no | exploratory baseline ladder (null, size, composition, indices, RDKit 2D) with ΔR², ΔMAE and 95 % CIs; the confirmatory protocol runs from the command line |
| Provenance | time taken per request | `run_record.json` per job (commit, versions, ChEMBL release, config hash, seeds), results zip, data-licence notes |
| Downloads | no | CSV/zip of indices, manifest, results |
| Structure view | graph plot | RDKit depiction and degree edge partition per structure |
| Input validation and limits | none in the server (malformed JSON gives HTTP 500; my first attempt with a plain matrix did) | request-size and molecule caps, per-IP rate limits, timeouts, strict parsing with error messages |
| Rendering of server output | via `innerHTML` (chart) | no `innerHTML`; `textContent` only (CSP-safe) |
| Persistence and health | stateless Lambda | SQLite job queue, one worker, unguessable job ids, `/health` |
| Tests / CI | none found | 1,101 offline tests incl. 36 web tests; CI |
| Accessibility | not assessed | every input labelled, live-region status, keyboard operable (tested) |
| Look and feel | blue header, three tabs, white cards, serif text | the same language on purpose (co-author request, 2026-10-05): side-by-side in `docs/WEB.md` |
| Cold start | AWS Lambda, 3.4 s observed 2026-10-04 | Render free tier sleeps after 15 min idle; the cold start will be slower than Lambda's — stated, not hidden |

Screenshots of every TopoQSPR view are in `docs/screenshots/` (01–10); `10_home_chemgraphx_style.png`
next to `reference_chemgraphx_home.png` in `docs/WEB.md` shows the matched layout.

## 3. Status and what is still open

* TopoQSPR is measured on a local deployment of the same code that `render.yaml` deploys; the
  public URL does not exist yet because deployment needs a Render account (the lead author has
  none; Advik's hosts the v0.x prototype). Once deployed, re-run `tools/compare_sites.py` with
  `--topoqspr https://<host>` and paste the output here.
* An adversarial review of the website by an independent agent is scheduled after the agent usage
  limit resets (2026-10-07); until then the evidence is this document, the test suite and the
  screenshots.

Verdict (lead author, 2026-10-06): the TopoQSPR site does everything the ChemGraphX site does on the
inputs both accept, returns identical numbers on the shared indices, and adds the inputs, the
disease workflow, the catalogue, the modelling, the downloads, the provenance, the validation and
the tests that ChemGraphX lacks; its one functional gap is the missing distance entropies, and its
one unverified property is cold-start latency on the chosen host.
