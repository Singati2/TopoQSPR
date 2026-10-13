# 01b — Data licences and measurement provenance (verified 2026-10-13)

Scope: every external dataset that TopoQSPR uses or cites. Each entry was checked against primary pages: repository APIs, licence files, dataset metadata APIs, and the data files themselves in `~/.cache/topoqspr/expdata/`. Nothing in the package was edited.

Confidence scale:
- **high**: licence and measurement both read from the primary source.
- **medium**: one of them is inferred or depends on an upstream source whose terms were not found.
- **low**: secondhand.

## Summary

| Dataset | Licence (as published) | Redistribution | Measured or predicted | `expdata.py` record correct? | Confidence |
|---|---|---|---|---|---|
| OPERA `OPERA_Data.zip` (LogP/MP/WS_QR.sdf) | Repo `LICENSE`: **MIT** (a software licence; no separate data licence in the repo or the zip). The 2018 article's data supplement falls under BMC's **CC0** data waiver. Upstream PHYSPROP (SRC, distributed by EPA in EPI Suite): no licence terms found | Probably yes, under MIT with the notice kept. Upstream terms are unknown | **Measured** literature values (PHYSPROP experimental), curated | **Partly.** "MIT (repository LICENSE)" is literally what the repo states, but it misdescribes the data's status (see §1) | medium |
| Bradley Open Melting Point (figshare 1031637 v2, file 1503990) | **CC0** | yes | **Measured**, aggregated (vendors, literature, DrugBank, PhysProp, …) | yes ("CC0") | high |
| AqSolDB (Dataverse doi:10.7910/DVN/OVHAW8 V1) | **CC0 1.0** | yes | **Measured**, aggregated from 9 public sets. One value per compound; for 2-value compounds the chosen value is the one closest to an ALOGPS prediction | yes ("CC0 1.0") | high |
| AstraZeneca deposition CHEMBL3301361 (via ChEMBL 37) | **CC BY-SA 3.0 Unported** (ChEMBL) | yes, with attribution and share-alike | **Measured** (AZ in-house shake-flask logD7.4; solid-start solubility, pH 7.4) | yes ("CC BY-SA 3.0 (ChEMBL)") | high |
| Biogen ADME (`ADME_public_set_3521.csv`) | Repo `LICENSE`: **MIT** (no separate data licence) | yes, under MIT with the notice kept | **Measured** (Biogen in-house in vitro assays) | **Yes**, with the same caveat that MIT is a software licence | high (licence text) / high (measured) |
| ChEMBL 37 (general) | **CC BY-SA 3.0 Unported** | yes, with attribution and share-alike | mixed (deposited and extracted measurements, plus computed properties) | n/a | high |
| QM9 (figshare collection 978904) | **CC0** on the data files. Article *Sci Data* 1:140022 is CC BY-NC-SA | yes | **Computed** (quantum chemistry), not measured | not used by the package (no `qm9` in code) | high (licence) |

---

## 1. OPERA data files

- **URL (pinned in `expdata.py`):** `https://github.com/kmansouri/OPERA/raw/22f4db60e38548f03921c44aecc84eb1b792e955/OPERA_Data.zip`
- **Repository licence:** `gh api repos/kmansouri/OPERA` → `license.spdx_id = "MIT"`. `LICENSE` text: "MIT License / Copyright (c) 2016 kmansouri / Permission is hereby granted … to any person obtaining a copy of this **software and associated documentation files** (the "Software") …" (read 2026-10-13; repo pushed 2026-10-09).
- **README:** "OPERA is a free and open-source/open-data suite of QSAR models …". It states no data-specific licence.
- **Zip contents:** 36 members (AOH…WS_QR.sdf, plus `*_2.9_update*.csv/xlsx`). **No licence or readme file inside the zip.** The QR SDFs used by TopoQSPR are dated 2019-12 to 2020-01.

**Article licence.** Mansouri K, Grulke CM, Judson RS, Williams AJ. *J Cheminform* 10:10 (2018), doi:10.1186/s13321-018-0263-1 (PMC5843579, read in full).
- Article: CC BY 4.0.
- "The Creative Commons Public Domain Dedication waiver (http://creativecommons.org/publicdomain/zero/1.0/) applies to the data made available in this article, unless otherwise stated." The curated training/test sets are "available in the Additional file 1: S1, a GitHub repository, ResearchGate, and the JRC's QMRF Inventory".
- Code: "OPERA executables … and the associated MATLAB source code are available for free on Github under the … (MIT) license … together with the data and QMRFs".
- **Caveat.** The CC0 waiver covers the 2018 article's supplementary data, which is the v1.5-era sets. The zip pinned by TopoQSPR is a later revision (files dated 2019-12/2020-01, repo tag v2.9.5). It is distributed only under the repository's MIT licence.

**Upstream provenance.**
- The values come from PHYSPROP. The article's ref. 53 is "EPA (2014) EPI Suite Data. http://esc.syrres.com/interkow/EpiSuiteData_ISIS_SDF.htm"; ref. 55 is "SRC (2016) Scientific databases"; curation is Mansouri et al. *SAR QSAR Environ Res* 27:911–937 (2016), doi:10.1080/1062936X.2016.1253611 (not OA, not read).
- EPA describes PhysProp (maintained by Syracuse Research Corporation) as the "experimental database" behind EPI Suite (EPA EPI Suite pages, via search summary [secondhand]).
- **No licence or redistribution terms for PHYSPROP data were found.**

**Measured or predicted.**
- **Measured.** The article calls them "experimental data" curated from PHYSPROP (Table 2: logP 15,806 → 14,041 QSAR-ready; MP 10,051 → 8,656; WS 5,764 → 4,224).
- Duplicates were averaged, or removed if they disagreed.
- Evidence in the files themselves:
  - `LogP_QR.sdf` (13,963 records) has a `Kow Reference` field on every record. Top sources: HANSCH,C ET AL. (1995) 3,936; BIOBYTE STARLIST (2009) 1,767; HANSCH & LEO (1985) 1,287; BIOBYTE (1995) 793; SANGSTER (1993/1994) 1,244. None contain "estim", "calc", "pred" or "kowwin" (grep count 0).
  - `WS_QR.sdf` (5,269) has `WS Reference` and `WS Temperature`.
  - `MP_QR.sdf` (8,698) has **no reference field**, so MP's "measured" status rests on the article (medium confidence).
- **Notes for the methods section:**
  - WS temperatures are not all 25 °C: 2,234 at 25, 1,028 at 20, 1,549 blank, the rest 21–30. The OPERA README says "Water solubility at 25C".
  - 160 MP, 90 WS and 91 logP records carry a `Salt_Solvent` component. TopoQSPR already drops multi-component solid-state records (`salt_or_solvate_form`, `expdata.py` ~l.34, ~l.724).
  - OPERA logP training data overlap with KOWWIN training data (article, logP comparison section). This is irrelevant to TopoQSPR, which uses the measured values.

**Is `licence="MIT (repository LICENSE)"` at `expdata.py` ~l.219 correct?**
- It is accurate as a description of what the repository declares.
- It is incomplete as a data licence statement: MIT is worded for software, and the upstream PHYSPROP terms are unknown.
- Suggested record: `"MIT (repository LICENSE; no data-specific licence). 2018 article data supplement: BMC CC0 waiver. Upstream PHYSPROP (SRC/EPA EPI Suite): terms not stated"`.
- `redistribution=True`: keep it, but because the 78 MB zip is never copied into the repo, this only matters for derived tables. If processed OPERA tables are ever published (e.g. on Zenodo), keep the MIT notice and cite Mansouri 2018 and 2016 plus PHYSPROP/EPI Suite.

**Attribution text:** "OPERA curated PHYSPROP data: Mansouri K, Grulke CM, Judson RS, Williams AJ. J Cheminform 10:10 (2018), doi:10.1186/s13321-018-0263-1; Mansouri K et al. SAR QSAR Environ Res 27:911–937 (2016), doi:10.1080/1062936X.2016.1253611; data from github.com/kmansouri/OPERA (MIT), commit 22f4db60e385."

## 2. Bradley Open Melting Point dataset

- **figshare API** `articles/1031637`:
  - title "Jean-Claude Bradley Open Melting Point Dataset";
  - **license {name: "CC0", url: creativecommons.org/publicdomain/zero/1.0/}**;
  - DOI 10.6084/m9.figshare.1031637.v2; published 2014-05-20; version 2;
  - single file `BradleyMeltingPointDataset.xlsx`, **file id 1503990** (so "figshare 1503990" is the file, not the article);
  - authors Bradley, Williams, Lang.
- Description: "28,645 measurements including those found to be incorrect (marked as 'do not use'). csid corresponds to Chemspider ID."
- **Measured:** yes. Columns: key, name, smiles, csid, link, source, donotuse, donotusebecause.
  - Source links by domain: alfa.com (Alfa Aesar) 8,782; drugbank.ca 875; scripts.iucr.org 262; chemspider.com 150; pubs.acs.org 52; NIST WebBook 37; and others.
  - Source labels include "Alfa Aesar", "PhysProp", "Oxford", and Karthikeyan, Glen & Bender *JCIM* 45:581 (2005) records.
  - These are aggregated literature and vendor measurements of variable quality, hence the donotuse flag.
- **Redistribution:** yes (CC0). The repo's test excerpts are fine.
- **Attribution (requested, not legally required):** "Bradley J-C, Williams A, Lang A. Jean-Claude Bradley Open Melting Point Dataset. figshare (2014). doi:10.6084/m9.figshare.1031637.v2".
- `expdata.py` (licence "CC0", doi …1031637.v2, url files/1503990): **correct.** Confidence high.

## 3. AqSolDB

- **Dataverse API** `datasets/:persistentId/?persistentId=doi:10.7910/DVN/OVHAW8`:
  - latest version 1.0, released 2019-07-18;
  - **license {name: "CC0 1.0", rightsIdentifier: "CC0-1.0"}**, termsOfUse empty;
  - file `curated-solubility-dataset.tab` id 3407241 (the `?format=original` download is the CSV).
- Description: "aqueous solubility values of 9,982 unique compounds curated from 9 different publicly available aqueous solubility datasets … reference source of measured and calculated solubility data". "Calculated" refers to the RDKit 2D descriptor columns, not the LogS values.

**What it aggregates** (Sorkun, Khetan, Er, *Sci Data* 6:143 (2019), doi:10.1038/s41597-019-0151-1, PMC6687799, read):
- A: eChemPortal (OECD), "experimental studies" filter, 6,110;
- B: EPI Suite Data website, 4,651;
- C: Raevsky et al., 2,603;
- D: EPI Suite Data website, 2,115;
- E: Huuskonen (AQUASOL + PHYSPROP), 1,291;
- F: Wang et al. (Beilstein), 1,210;
- G: Delaney (ESOL), 1,144;
- H: Wang/Jain–Yalkowsky, 578;
- I: Goodman/Llinas solubility challenge, 94.

Temperature filter 25 ± 5 °C where reported.

**Measured or predicted:** measured ("experimental aqueous solubility value (LogS)", Table 2).
- Caveat from Methods: for compounds with exactly two differing values (1,417), "we selected the closest value to an external reference value, which is obtained using the solubility prediction tool ALOGPS".
- The retained value is experimental, but a model chose which one. Groups G1–G5 encode reliability.
- AqSolDB inherits EPI Suite/PHYSPROP solubility data (sets B, D, E), so **it overlaps OPERA WS.** Do not treat the two as independent replications.

**Redistribution:** yes (CC0). Article CC BY 4.0. Curation code MIT (Code Ocean).

**Attribution:** "Sorkun MC, Khetan A, Er S. AqSolDB, a curated reference set of aqueous solubility and 2D descriptors for a diverse set of compounds. Sci Data 6:143 (2019), doi:10.1038/s41597-019-0151-1; data doi:10.7910/DVN/OVHAW8 (V1)."

`expdata.py`: **correct.** Confidence high.

## 4. AstraZeneca ChEMBL deposition CHEMBL3301361

- **ChEMBL API** `document/CHEMBL3301361`:
  - title "Experimental in vitro DMPK and physicochemical data on a set of publicly disclosed compounds"; authors Mark Wenlock and Nicholas Tomkinson; doc_type DATASET; src_id 27; DOI 10.6019/CHEMBL3301361; first in CHEMBL_20 (2015-01-14).
  - Abstract: "Experimental data determined at AstraZeneca on a set of compounds in the following assays: pKa, lipophilicity (LogD7.4), aqueous solubility, plasma protein binding …, intrinsic clearance …"
- **Assays:**
  - `CHEMBL3301363`: "ASTRAZENECA: Octan-1-ol/water (pH7.4) distribution coefficent measured by a shake flask method described in J. Biomol. Screen. 2011, 16, 348-355. Experimental range -1.5 to 4.5".
  - `CHEMBL3301364`: "Solubility in pH7.4 buffer using solid starting material using the method described in J. Assoc. Lab. Autom. 2011, 16, 276-284. Experimental range 0.10 to 1500 uM".
  - `expdata.py` `assay_range` values match these: logD (−1.5, 4.5); solubility 100 nM to 1.5136 mM.
- **Measured:** yes (high).
- This deposition is also the source of MoleculeNet "Lipophilicity" (4,200 compounds, random split recommended; Wu et al. 2018, Table 1).
- **Licence:** ChEMBL's CC BY-SA 3.0 Unported (§6).
- **Redistribution:** allowed with attribution. Derived tables must stay CC BY-SA 3.0 (share-alike).

**Attribution:** "AstraZeneca. Experimental in vitro DMPK and physicochemical data on a set of publicly disclosed compounds. ChEMBL deposited dataset CHEMBL3301361, doi:10.6019/CHEMBL3301361. ChEMBL data is from http://www.ebi.ac.uk/chembl - the version of ChEMBL is chembl_37."

`expdata.py`: **correct.**

## 5. Biogen ADME dataset

- **Repository:** `gh api repos/molecularinformatics/Computational-ADME` → `license.spdx_id = "MIT"`; `LICENSE`: "MIT License / Copyright (c) 2020 molecularinformatics / … software and associated documentation files …". Pushed 2024-02-15. Top level: `ADME_public_set_3521.csv`, FCNN, ML, MPNN, LICENSE, README.md.
- **README:** "we have disclosed a collection of 3521 diverse compounds selected from commercially available compound libraries (i.e. Enamine, eMolecules, WuXi LabNetwork, Mcule) and tested them against our internal six ADME in vitro assays … The file 'ADME_public_set_3521.csv' contains … the experimental log(properties) for six endpoints: HLM, RLM, Solubility, MDR1-MDCK ER, hPPB, and rPPB".
- **Measured:** yes (high).
- **Paper:** Fang C et al. *JCIM* 63:3263–3274 (2023), doi:10.1021/acs.jcim.3c00160. Crossref licence links are ACS standard policies (10.15223/policy-029/037/045); Europe PMC isOpenAccess = N. The paper itself was not read.
- **Redistribution:** yes, under MIT with the copyright and permission notice kept.

**Attribution:** "Fang C, Wang Y, Grater R, et al. Prospective Validation of Machine Learning Algorithms for ADME Prediction: An Industrial Perspective. J Chem Inf Model 63:3263–3274 (2023), doi:10.1021/acs.jcim.3c00160; data github.com/molecularinformatics/Computational-ADME (MIT), commit 685c9c828b23."

**Is `licence="MIT (repository LICENSE)"` at `expdata.py` ~l.231 correct?** **Yes**, as a description of the only licence the data carry. Optionally add "(no data-specific licence)". The OPERA upstream ambiguity does not apply: Biogen generated the data itself.

## 6. ChEMBL 37 licence and attribution

- **Release:** `chembl_37`, dated 01/05/2026, i.e. 1 May 2026 (FTP `README` "Date: 01/05/2026"; API `status.json`: chembl_db_version "ChEMBL_37", chembl_release_date "2026-05-01").
- **LICENSE file** (`ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/LICENSE`): full text of "Creative Commons Attribution-ShareAlike 3.0 Unported".
- **README:** "The data content in ChEMBL is licensed under … 'CC Attribution-ShareAlike 3.0 Unported license' … The required attribution should contain the url of the ChEMBL resource, and also the release version, e.g.: *ChEMBL data is from http://www.ebi.ac.uk/chembl - the version of ChEMBL is chembl_37.* This should be visible on the entry portal for a web resource in which ChEMBL is integrated, or contained with the documentation for any further distribution."
- **REQUIRED.ATTRIBUTION** (verbatim excerpts):
  - "Under the -BY clause, we request attribution for subsequent use of ChEMBL."
  - Primary citation: "Mendez D, Gaulton A, Bento AP, … Leach AR. ChEMBL: towards direct deposition of bioassay data. Nucleic Acids Res. 2019 47(D1):D930-D940. DOI: 10.1093/nar/gky1075".
  - "If ChEMBL is incorporated into other works, we ask that the ChEMBL IDs are preserved, and that the release number of ChEMBL is clearly displayed."

**Implications for TopoQSPR:**
1. The web UI and any published cohort manifest must show "ChEMBL data is from http://www.ebi.ac.uk/chembl - the version of ChEMBL is chembl_37", or whichever release is used.
2. Keep ChEMBL IDs in manifests.
3. Any redistributed derived ChEMBL table (cohort CSVs, AZ logD records) must be released under CC BY-SA 3.0, not CC0 or MIT. Keep these separate from MIT-licensed code in a Zenodo deposit.
4. Cite Mendez 2019, the file's stated primary citation. Optionally also cite Zdrazil et al. *NAR* 52:D1180 (2024), doi:10.1093/nar/gkad1004.

## 7. QM9

- **figshare collection 978904:** "Quantum chemistry structures and properties of 134 kilo molecules", DOI 10.6084/m9.figshare.c.978904.v5 (2019-06-18). The collection object has no licence field.
- Member articles carry the licence:
  - 1057646 "Data for 133885 GDB-9 molecules" (`dsgdb9nsd.xyz.tar.bz2`): **CC0**;
  - 1057645 "Data for 6095 constitutional isomers of C7H10O2" (`dsC7O2H10nsd.xyz.tar.bz2`): **CC0**.
- **Article:** Ramakrishnan R, Dral PO, Rupp M, von Lilienfeld OA, *Sci Data* 1:140022 (2014), doi:10.1038/sdata.2014.22. Europe PMC licence "cc by-nc-sa".
- **Measured:** no. These are computed quantum-chemical properties. Do not present them as experimental.
- **Redistribution:** yes (CC0).
- **Attribution:** cite Ramakrishnan 2014 and the figshare DOI.
- **Use in package:** none found (grep for `qm9`/`978904` in `*.py`/`*.md` outside the audit folder returned nothing). It is relevant only as a formula-controlled computed comparator for the isomer study (see `02_NOVELTY_PRIOR_ART_MATRIX.md` §2.11).

---

## Corrections to make (suggestions only; nothing edited)

1. `expdata.py` ~l.219 (OPERA): extend the licence string as in §1. Do not change `redistribution`, but note the unresolved upstream PHYSPROP terms in `docs/EXPERIMENTAL_DATA.md`. Its "Licence notes" already mention PHYSPROP/SRC; add the BMC CC0 statement and its version caveat.
2. `expdata.py` ~l.231 (Biogen): correct as is. Optionally append "(no data-specific licence)".
3. `docs/EXPERIMENTAL_DATA.md` and the paper's data-availability statement:
   - add the exact ChEMBL attribution sentence and the Mendez 2019 citation;
   - state that ChEMBL-derived tables are CC BY-SA 3.0.
4. Methods text:
   - AqSolDB and OPERA WS share EPI Suite/PHYSPROP provenance, so they are not independent;
   - AqSolDB's two-value tie-break uses ALOGPS;
   - OPERA WS temperatures range 20–30 °C, with 1,549 unrecorded.

## Access log (2026-10-13)

- `gh api repos/kmansouri/OPERA` (+ `/license`, `/readme`, `/contents`);
- `unzip -l` and field greps on the cached `OPERA_Data.zip` (sha256 as pinned);
- Europe PMC fullTextXML PMC5843579 (OPERA 2018) and PMC6687799 (AqSolDB);
- figshare API `articles/1031637`, `collections/978904`, `collections/978904/articles`, `articles/1057645`, `articles/1057646`;
- XML of the cached Bradley xlsx (`xl/sharedStrings.xml`);
- Harvard Dataverse API for doi:10.7910/DVN/OVHAW8;
- ChEMBL REST `document/CHEMBL3301361`, `assay/CHEMBL3301363`, `assay/CHEMBL3301364`, `status.json`;
- ChEMBL FTP `chembl_37/{LICENSE, README, REQUIRED.ATTRIBUTION, chembl_37_release_notes.txt}`;
- `gh api repos/molecularinformatics/Computational-ADME` (+ `/license`, `/readme`);
- Crossref licence fields for 10.1021/acs.jcim.3c00160;
- WebSearch for PHYSPROP licence terms (none found; EPA pages via search summary).
