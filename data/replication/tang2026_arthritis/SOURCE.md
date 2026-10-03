# Source and licence

**Article.** Tang, Jiang-Hua; Noureen, Sadia; Ashraf, Nazma; Azam, Saood; Aslam, Adnan. Degree-Based Topological Indices and Machine Learning for QSPR Modeling of Arthritis Drugs. *ACS Omega* 11 (2026) 33216. doi:[10.1021/acsomega.6c03590](https://doi.org/10.1021/acsomega.6c03590). PMC13261421, PMID 42294207.

**Shortlist row** 7 of `docs/REPLICATION_SHORTLIST.md`.

**Retrieved** 2026-10-03 from Europe PMC (`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13261421/fullTextXML`); the SHA-256 of the XML is in `study.json` (`source.xml_sha256`).

**Article licence (quoted from the JATS `<permissions>` element):**

> © 2026 The Authors. Published by American Chemical SocietyThis article is licensed under CC-BY-NC-ND 4.0

**What is stored here.** The article is licensed CC BY-NC-ND 4.0, which does not permit sharing adapted material. To stay on the safe side this folder stores only the drug names, the structure resolution, the table layout (table and column identifiers), the index formulas and the reported fit statistics needed to cite the published result. The per-drug property values are **not** committed. `python tools/extract_replication.py fetch-xml <study>` followed by `extract <study>` rebuilds `properties.csv` in the local cache (`~/.cache/topoqspr/replication/`) from the Europe PMC XML.

**Our additions** (not from the article): PubChem CIDs, SMILES and InChIKeys (PubChem PUG-REST, public domain), resolution rules and notes, consistency checks.

**Attribution.** Data extracted from Tang, Jiang-Hua; Noureen, Sadia; Ashraf, Nazma; Azam, Saood; Aslam, Adnan (2026), doi:10.1021/acsomega.6c03590.
