# Cohort definition (v1.0)

How TopoQSPR turns a disease into a list of drug structures, and why each default is what it
is. Code: `src/topoqspr/chembl.py` (REST client and cache) and `src/topoqspr/cohort.py`
(cohort builders). Tests: `tests/test_cohort.py`, offline from fixtures recorded on ChEMBL_37.

All counts in this document were measured on **ChEMBL_37** (release date 2026-05-01, read
from `/status`) on 2026-10-02. The reproduction snippet is in section 8.

## 1. Procedure (controlled cohort, `build_cohort`)

1. **Release.** Read `chembl_db_version` from `/status` once per client. Every cache entry
   is keyed by (release, endpoint, parameters) and stored under `<cache>/<release>/`, so
   data from two releases never mix. New entries are first written to a private staging
   directory (`<cache>/.staging-*/`). At the end of a build the client re-reads `/status`.
   If the release is unchanged, the staged entries are moved into `<cache>/<release>/`. If
   it changed, they are deleted and the build raises, so data served by a newer release is
   never filed under the old one. A client that never calls `verify_release()` (or exits
   its `with` block with an exception) promotes nothing.
2. **Disease to MeSH IDs.** `suggest_headings(free_text)` ranks the distinct
   (`mesh_id`, `mesh_heading`) pairs that occur in `drug_indication` for this release.
   The score is `0.6 * token overlap + 0.4 * difflib ratio` on lower-cased tokens.
   Optionally, it also scores MeSH entry terms (synonyms) from NLM's lookup API. Terms are
   looked up with `match=contains` on the text and with `match=exact` on its MeSH-style
   inversions ("high blood pressure" → "blood pressure, high", which is an entry term of
   Hypertension). A match through an entry term is multiplied by 0.95, so a heading that
   matches as well ranks first ("diabetes" → *Diabetes Mellitus* before *Hemochromatosis*,
   whose entry terms include "Bronze Diabetes"). The output is only a suggestion. **The
   cohort is defined only by the descriptor IDs the user confirms** (`D######`). Free
   text, IDs that are not descriptor IDs, and IDs with no rows in the release are all
   rejected. A single ID may be passed as a string; IDs are stripped and upper-cased.
   There is no substring retrieval and no keyword fallback.
3. **Optional tree expansion** (`include_descendants=True`, default off). Descendants of
   each ID are read from the NLM MeSH SPARQL endpoint
   (`meshv:treeNumber/meshv:parentTreeNumber+`). Only descendants that have rows in this
   release are added. All of them are recorded in the sidecar
   (`mesh_descendants_added`, `mesh_descendants_without_rows`).
4. **Indication rows.** Fetch every `drug_indication` row with `mesh_id` equal to each
   confirmed ID, fully paginated and ordered by `drugind_id`. The client checks that the
   number of records equals `total_count`.
5. **Evidence filter (optional, default off).** A row *qualifies* unless it is removed by
   * `evidence_sources=S`: at least one `indication_refs.ref_type` must be in `S`, or
   * `exclude_sole_evidence=E`: the row is dropped if *all* its ref types are in `E`.
     For example, `{"DailyMed"}` drops indications that are supported only by
     text-mined label sections. A row with no refs counts as sole-evidence.

   A single string is one ref type. Every name must be one of the seven ref types that
   occur in ChEMBL_37 (`ATC, ClinicalTrials, DailyMed, EMA, FDA, INN, USAN`; together they
   cover all 60,055 indication rows). Unknown or misspelled names raise `ValueError`, so a
   typo cannot turn the filter into a no-op or an exclude-all.
6. **Parent mapping.** Each molecule is mapped to `molecule_hierarchy.parent_chembl_id`.
   If a record has no hierarchy, the mapping falls back to
   `drug_indication.parent_molecule_chembl_id`, then to the molecule itself, and is
   followed to a fixed point. In ChEMBL_37 the two parent sources agree for every row of
   the four diseases below (0 mismatches; the count is stored in the sidecar).
7. **Phase of a parent.**
   `phase(P) = max{ max_phase_for_ind(r) : r a qualifying row of a molecule in P's family }`,
   where the family is the parent plus its salt and hydrate forms. ChEMBL phases are
   decimals (0.5 = early phase 1); `-1`/null mean unknown and are treated as missing, never
   as low. A parent with no known phase is excluded as `phase_unknown`. Keep the parent
   if `phase >= min_phase` (default 4).
8. **Molecule type.** The parent record must have `molecule_type == 'Small molecule'`.
9. **Structure (RDKit, on the parent's ChEMBL `canonical_smiles`).** The first failing
   check is reported:
   `no_smiles` → `unparseable_smiles` → `no_carbon` → `metal_containing` →
   `unsupported_element` → `mixture` (more than one fragment left in the *parent*) →
   `inorganic_carbon` (at most 3 heavy atoms and no carbon bonded to H or to another
   carbon: elemental carbon, CO, CO2, cyanide, CS2) → `too_small` (one heavy atom, e.g.
   methane, a graph with no edges). Allowed elements: H, B, C, N, O, F, Si, P, S, Cl, Se,
   Br, I. Urea and CCl4 pass (the 3-atom cap keeps them). A parent that passes these
   checks but has ChEMBL `inorganic_flag == 1` is excluded as `inorganic_flag`.
   Optional sensitivity rule `exclude_diagnostics=True` (default off): exclude as
   `diagnostic_agent` a parent whose family's ATC codes are *all* V04 (diagnostic agents)
   or V09 (radiopharmaceuticals), or whose structure has a non-deuterium isotope label
   (e.g. [18F]). "All" keeps tolbutamide (A10BB03 and V04CA01) and folic acid.
10. **De-duplication.** Among the remaining parents, keep one per full standard InChIKey,
    the one with the lowest ChEMBL ID. The others get `duplicate_inchikey:<kept ID>`.
11. **Order.** Rows are sorted by the numeric part of the parent ChEMBL ID. The parent's
    own row comes first, then its family members by ChEMBL ID.

Each candidate molecule gets exactly one manifest row:

* each parent gets one decision row (`molecule_chembl_id == parent_chembl_id`), and
* each salt or hydrate form gets an `exclude / collapsed_to_parent` row.

A parent that has no indication row of its own (its phase comes only from a salt) still
gets a decision row.

## 2. Defaults and why

| Setting | Default | Justification |
|---|---|---|
| Disease identity | explicit MeSH descriptor IDs | `mesh_heading__icontains` mixes diseases: "diabetes" matches 8 headings, including Diabetes Insipidus (section 6). `drug_indication/search` returns HTTP 404, so v0.x's keyword fallbacks never returned anything. |
| `include_descendants` | off | MeSH is a polyhierarchy. Below *Breast Neoplasms* (D001943), SPARQL returns *Carcinoma, Lobular* and *Hereditary Breast and Ovarian Cancer Syndrome* as well as the subtypes. Expansion must be an explicit, recorded choice. |
| `min_phase` | 4 | The plan (EXPERIMENT_PLAN §3.1) studies approved drugs for the indication. With all phases, the T2D cohort grows from 33 to 256 parents, mostly phase 1–3 investigational compounds. |
| Phase aggregation | max over the parent's salt family | ChEMBL often records the phase on a salt. Without this rule, 8 of 33 T2D parents would be lost (metformin, sitagliptin, pioglitazone, alogliptin, …), as would 29 of 85 for hypertension, 15 of 36 for asthma and 16 of 35 for breast neoplasms. |
| Evidence filter | none | Dropping DailyMed-only indications removes some real noise (biotin, tyrosine, pyridoxine for T2D). It also removes genuine drugs: tolbutamide, chlorpropamide and tolazamide for T2D; propranolol, verapamil, labetalol, nadolol and fosinopril for hypertension. Requiring FDA/EMA/ATC refs is worse: it drops albuterol, budesonide and fluticasone for asthma. Section 5 has the measurements. Both filters are available for sensitivity analyses. |
| Molecule type | `Small molecule` | Topological indices are defined on small-molecule graphs. Proteins, antibodies and oligosaccharides are excluded. |
| Structure source | ChEMBL parent `canonical_smiles` | The ChEMBL parent is already salt-stripped by the ChEMBL structure pipeline. Using it keeps the manifest traceable to one ChEMBL record. |
| Metal salts | excluded (`metal_containing`) | Some ChEMBL parents keep a metal counter-ion, for example sodium butyrate, potassium citrate and zinc gluconate. We exclude these rather than strip the ion silently. **Deviation from EXPERIMENT_PLAN §3.1**, which says "parent compound → largest organic fragment": see the note below the table. |
| Mixtures | excluded | A multi-fragment *parent* is a co-crystal or multi-component product, such as omega-3-acid ethyl esters or ipragliflozin L-proline. Its graph is not one molecule. Same deviation as above. |
| De-duplication key | full standard InChIKey | Stereoisomers are distinct compounds with distinct properties, so we do not merge on the first InChIKey block. This choice does not affect any of the four cohorts below: there were 0 InChIKey duplicates. Pairs that share the first block (identical 2D graph) are kept, and the pipeline keeps each such pair in one CV fold (PIPELINE.md §1, *Connectivity groups*). |
| Tie-break / order | lowest ChEMBL ID; numeric ChEMBL ID sort | Deterministic and independent of API order and of `PYTHONHASHSEED`. |

**Deviation from the experiment plan.** EXPERIMENT_PLAN §3.1 specifies "parent compound →
largest organic fragment". The builder instead *excludes* parents that still have several
fragments (`mixture`) or a metal (`metal_containing`), as the cohort task specification
asked. The two documents conflict; we followed the task specification. The reason: the
ChEMBL parent is already salt-stripped, so a fragment that remains is part of the
registered substance, and reducing it to its largest fragment would assign the record's
properties to a different chemical entity without trace. At the default phase 4 this
choice changes no cohort (0 mixture or metal exclusions in all four diseases). With all
phases it removes 3 (T2D), 5 (hypertension), 1 (asthma) and 3 (breast neoplasms) parents
(section 5). The plan should be amended to match before it is frozen.

## 3. Manifest

`Cohort.write(csv_path, timestamp_utc)` writes two files:

* **CSV**, UTF-8 with `\n` line endings. Columns: `disease_label, mesh_ids, chembl_release,
  molecule_chembl_id, parent_chembl_id, pref_name, canonical_smiles, standard_inchi_key,
  molecule_type, phase_for_indication, indication_ref_types, efo_ids, decision, reason`.
  * `phase_for_indication` is formatted with one decimal and is empty for unknown.
  * `indication_ref_types` is the `;`-joined set of ref types over the family's rows. It
    lists all rows, not only those that qualified.
  * `efo_ids` is the `;`-joined set of `drug_indication.efo_id` values over the same rows
    (e.g. T2D rows carry MONDO:0005148, and two carry Orphanet:552, MODY). The sidecar
    maps each EFO ID to its terms (`efo_terms`).
  * `reason` is one of `included, collapsed_to_parent, no_qualifying_evidence,
    phase_unknown, phase_below_min, molecule_record_missing, not_small_molecule,
    no_smiles, unparseable_smiles, no_carbon, metal_containing, unsupported_element,
    mixture, inorganic_carbon, too_small, inorganic_flag, diagnostic_agent,
    duplicate_inchikey:<ID>`.
* **JSON sidecar.** It holds:
  * the kind and the release with its date;
  * every parameter;
  * the MeSH headings used;
  * counts per reason;
  * the number of indication rows and the parent-mismatch count;
  * `topoqspr_version`;
  * `created_utc`, which is **passed in by the caller**;
  * the SHA-256 of the CSV.

Same inputs and the same release give a byte-identical CSV. The tests check this
in-process and across two interpreter runs with different `PYTHONHASHSEED`.

## 4. Alternative definitions

* **`build_atc_cohort(prefixes)`** selects molecules with a WHO ATC level-5 code that
  starts with any of the prefixes (ChEMBL `molecule.atc_classifications`). The parent,
  structure and de-duplication rules are the same; there is no phase rule.
  `indication_ref_types` lists the matched codes as `ATC:<code>`. It is used for
  cohort-sensitivity experiments. ATC classifies drugs by therapeutic use, not by
  indication: R03 also covers COPD, and L02 is endocrine therapy for any cancer. It is a
  proxy, not a gold standard.
* **`naive_cohort(free_text)`** is the v0.x retrieval frozen for Experiment A:
  * it uses `mesh_heading__icontains`;
  * every distinct `molecule_chembl_id` whose ChEMBL SMILES parses in RDKit is included,
    across all phases and all molecule types, with salts kept as separate molecules.

  It differs from v0.x only where v0.x was nondeterministic or dead code:
  * there is no `max_drugs` cut of an unordered `set`;
  * there is no early stop at `2 * max_drugs` IDs;
  * the output is sorted by ChEMBL ID;
  * the keyword fallbacks are omitted, because `drug_indication` is not searchable
    (HTTP 404 on `/drug_indication/search`, checked 2026-10-02).

## 5. Measured counts, ChEMBL_37

Naive search terms are the ones used in the v0.x paper draft (Table 1). The ATC proxies are:
T2D `A10`; hypertension `C02+C03+C07+C08+C09`; asthma `R03`; breast neoplasms `L02`.

**Cohort sizes**

| Disease (MeSH ID) | indication rows | naive: rows / molecules / included (distinct parents) | controlled | controlled, no DailyMed-only | controlled, FDA/EMA/ATC refs only | controlled, all phases | ATC proxy |
|---|---|---|---|---|---|---|---|
| Diabetes Mellitus, Type 2 (D003924) | 376 | "diabetes": 821 / 580 / 423 (390) | **33** | 26 | 19 | 256 | 48 |
| Hypertension (D006973) | 277 | "hypertension": 550 / 407 / 370 (318) | **85** | 69 | 38 | 194 | 167 |
| Asthma (D001249) | 269 | "asthma": 288 / 275 / 209 (190) | **36** | 28 | 8 | 177 | 57 |
| Breast Neoplasms (D001943) | 1171 | "breast neoplasms": 1407 / 868 / 624 (570) | **35** | 34 | 22 | 524 | 24 |

The naive queries match these headings:

* **"diabetes"** matches 8 headings: T2D, T1D, Diabetes Mellitus, Diabetes Insipidus,
  Diabetes Insipidus Neurogenic, Gestational, LADA, and Diabetes Complications.
* **"hypertension"** matches 11 headings, including Ocular, Portal, Intracranial and
  Pulmonary hypertension.
* **"asthma"** matches 4 headings.
* **"breast neoplasms"** matches 4 headings: the main heading plus Male, Inflammatory and
  Triple Negative.

**Controlled exclusions (default settings)**

| Disease | candidates | parents | collapsed salts | phase < 4 | not small molecule | structural |
|---|---|---|---|---|---|---|
| T2D | 389 | 353 | 36 | 305 | 13 (11 protein, 1 oligosaccharide, 1 unknown) | 2 (cyanocobalamin: no SMILES; omega-3-acid ethyl esters: mixture) |
| Hypertension | 296 | 239 | 57 | 150 | 3 | 1 (nitroprusside: no SMILES) |
| Asthma | 289 | 251 | 38 | 208 | 7 | 0 |
| Breast Neoplasms | 881 | 785 | 96 | 736 | 14 | 0 |

At phase 4 none of the four cohorts has a metal-containing, no-carbon or duplicate
exclusion. Those rules matter when the phase rule is off, because the phase check runs
first and is the reason reported. With `min_phase=None`, the structural exclusions are:

| Disease | metal_containing | no_carbon | mixture | no_smiles | inorganic_carbon |
|---|---|---|---|---|---|
| T2D | 1 (zinc gluconate) | 3 (e.g. sodium chloride, magnesium oxide) | 2 (e.g. ipragliflozin L-proline) | 5 | 1 (activated charcoal) |
| Hypertension | 4 (e.g. sodium butyrate, potassium citrate) | 6 (e.g. calcium, potassium chloride) | 1 | 1 | 0 |
| Asthma | 1 | 4 | 0 | 0 | 1 (activated charcoal) |
| Breast Neoplasms | 2 | 9 | 1 | 18 | 0 |

Activated charcoal (CHEMBL2106049, SMILES `[C]`, ChEMBL `inorganic_flag` 0) was included
before the `inorganic_carbon` rule existed; its ATC code A07BA01 would also have put it in an
ATC A07 cohort. No
parent in these cohorts is caught by `too_small` or by `inorganic_flag` alone.

The rules keep boron (dutogliptin) and selenium (ebselen) compounds.

**Effect of the DailyMed-only filter (`exclude_sole_evidence={"DailyMed"}`)**

| Disease | drugs dropped |
|---|---|
| T2D (−7) | biotin, pyridoxine, tyrosine, fenofibric acid, **chlorpropamide, tolazamide, tolbutamide** |
| Hypertension (−16) | atorvastatin, brimonidine, latanoprost, **acebutolol, betaxolol, eprosartan, esmolol, fenoldopam, fosinopril, labetalol, metolazone, moexipril, nadolol, nisoldipine, propranolol, verapamil** |
| Asthma (−8) | menthol, betamethasone acetate, cortisone acetate, methylprednisolone acetate, **ephedrine, metaproterenol, terbutaline, zafirlukast** |
| Breast Neoplasms (−1) | methyltestosterone |

Bold entries are standard therapies for the disease, so the filter is not a clean noise
filter. It stays off by default and is reported as a sensitivity analysis.

**Overlap (Jaccard)**

| Disease | naive vs controlled, InChIKey of the structures used | naive vs controlled, parent ChEMBL ID | controlled ∩ naive (InChIKey) | ATC vs controlled, InChIKey | controlled ∩ ATC |
|---|---|---|---|---|---|
| T2D | **0.073** | 0.085 | 31 of 33 | 0.397 | 23 |
| Hypertension | **0.213** | 0.267 | 80 of 85 | 0.457 | 79 |
| Asthma | **0.129** | 0.189 | 28 of 36 | 0.240 | 18 |
| Breast Neoplasms | **0.053** | 0.061 | 33 of 35 | 0.135 | 7 |

The two naive columns use different keys:

* The **InChIKey** column compares the structures each pipeline would actually hand to
  feature computation. Naive keeps salt structures; controlled uses parent structures.
* The **parent ChEMBL ID** column maps the naive molecules to their ChEMBL parents first.

The naive cohort has 4.4–17.8× as many included molecules as the controlled cohort
(370/85 for hypertension up to 624/35 for breast neoplasms). Under either key the Jaccard
overlap is at most 0.27. Most controlled parent structures also appear in the naive set:
31 of 33 (94%, T2D), 80 of 85 (94%, hypertension), 28 of 36 (78%, asthma) and 33 of 35
(94%, breast neoplasms), i.e. 78–94%. The rest appear there only as a salt structure. The
difference comes from what the naive query adds, not from what it misses.

**Independent cross-check.** A separate script, which does not use `topoqspr`, rebuilt the
four default cohorts. It used the REST filter `max_phase_for_ind=4`, the
`drug_indication.parent_molecule_chembl_id` field and the same structure rules. It gave
identical parent sets: 33, 85, 36 and 35.

## 6. Known limitations

* **Prodrug and ester parents** are separate ChEMBL parents, so they are not merged with
  their active moiety. `molecule_hierarchy.active_chembl_id` differs from the parent for
  1 of 33 (T2D), 14 of 85 (hypertension), 8 of 36 (asthma) and 3 of 35 (breast) included
  parents. Both members of the pair are in the cohort in five cases:
  * fenofibrate and fenofibric acid;
  * dexamethasone phosphate and dexamethasone;
  * hydrocortisone hemisuccinate and hydrocortisone;
  * methylprednisolone hemisuccinate and methylprednisolone;
  * prednisone and prednisolone.

  We keep the administered chemical entity because physicochemical QSPR targets belong to
  it. A scaffold-grouped CV sensitivity (EXPERIMENT_PLAN §6) covers the analog leakage.
* **A phase-4 indication record is not "a therapy for disease X".** The default cohorts
  contain agents whose ChEMBL indication mapping is diagnostic, supportive or about a
  related condition:
  * asthma: methacholine (bronchoprovocation diagnostic), menthol;
  * breast neoplasms: fluoroestradiol (the [18F] PET tracer), dexrazoxane
    (cardioprotectant);
  * hypertension: ambrisentan, bosentan, macitentan, riociguat and sitaxentan (pulmonary
    arterial hypertension), and atorvastatin, brimonidine and latanoprost (comorbidity or
    ocular hypertension labels);
  * T2D: ticagrelor, irbesartan, finerenone and folic acid (cardiovascular or renal
    complications, supplementation).

  We keep them in the default cohort because removing them needs a judgement per drug
  that ChEMBL does not record. `exclude_diagnostics=True` removes the clear cases
  (diagnostic or radiopharmaceutical ATC only, or an isotope label): methacholine
  (asthma 36 → 35) and fluoroestradiol (breast 35 → 34); T2D and hypertension are
  unchanged. With all phases it also removes fludeoxyglucose and piflufolastat (breast
  524 → 521). Comorbidity drugs are not caught by any rule and are listed here as known
  noise.
* **The evidence filter acts on rows, not on references.** ChEMBL does not attach a phase
  to each reference, so the filter cannot tell which reference supports phase 4.
* **Synonym suggestions** use NLM `lookup/term` with `match=contains`, capped at 50 terms.
  They help only when the free text is a contiguous part of an entry term. "breast
  cancer" finds D001943 and, through the inversion lookup, "high blood pressure" finds
  Hypertension. "adult-onset diabetes" still does not find D003924, whose entry term is
  "Diabetes Mellitus, Adult-Onset" (no contiguous or inverted match), and abbreviations
  such as "T2D" return nothing. For "diabetes", Hemochromatosis ("Bronze Diabetes")
  still appears in the top 4, below Diabetes Mellitus and Diabetes Insipidus.
  Suggestions never define a cohort.
* **The heading index** needs a full scan of `drug_indication`: 61 pages, about 45 s live,
  then cached per release.
* **MeSH version.** NLM serves the current MeSH year. Descendant lists are cached by
  descriptor but not by MeSH year, so delete `<cache>/mesh` when MeSH is updated.

## 7. Cache and offline use

* `ChemblClient(cache_dir, offline=False)`:
  * retries connection errors, timeouts, HTTP 429/5xx and non-JSON 200 responses with
    backoff `backoff·2^attempt`, honouring `Retry-After`;
  * raises immediately on other HTTP errors;
  * never caches a failure.
* Cache files are gzip JSON with `mtime=0` and sorted keys. Each file stores its own key,
  which is verified when the file is read.
* An online client stages new entries and promotes them only after `verify_release()`
  succeeds (every builder calls it; `with ChemblClient(...) as c:` calls it on exit).
  Leftover `.staging-*` directories from a crashed run are never read and can be deleted.
* `offline=True` serves from the cache only and raises `CacheMiss` on a miss. The test
  fixtures in `tests/fixtures/cohort/` are such a cache. Re-record them with
  `PYTHONPATH=src python tests/fixtures/cohort/record.py`.

## 8. Reproduction

```python
from topoqspr.chembl import ChemblClient
from topoqspr.cohort import build_cohort, build_atc_cohort, naive_cohort, jaccard

c = ChemblClient()                      # ~/.cache/topoqspr/chembl
ctl = build_cohort(["D003924"], c)      # Diabetes Mellitus, Type 2
nodm = build_cohort(["D003924"], c, exclude_sole_evidence={"DailyMed"})
reg = build_cohort(["D003924"], c, evidence_sources={"FDA", "EMA", "ATC"})
nav = naive_cohort("diabetes", c)
atc = build_atc_cohort("A10", c)
print(len(ctl.included), len(nodm.included), len(reg.included), len(atc.included))
print(jaccard(nav.included.standard_inchi_key, ctl.included.standard_inchi_key))
ctl.write("t2d_manifest.csv", timestamp_utc="2026-10-02T00:00:00Z")
```

Repeat with D006973 / "hypertension" / C02, C03, C07, C08, C09; D001249 / "asthma" /
R03; and D001943 / "breast neoplasms" / L02.
