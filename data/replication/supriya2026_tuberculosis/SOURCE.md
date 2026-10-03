# Source and licence

**Article.** Supriya; Iyer, Radha Rajamani. Predictive performance and QSPR analysis of SARS-CoV-2 and tuberculosis drugs using distance-based topological descriptors. *Scientific Reports* 16 (2026) 20043. doi:[10.1038/s41598-026-49132-8](https://doi.org/10.1038/s41598-026-49132-8). PMC13319808, PMID 42062328.

**Shortlist row** 10 of `docs/REPLICATION_SHORTLIST.md`.

**Retrieved** 2026-10-03 from Europe PMC (`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13319808/fullTextXML`); the SHA-256 of the XML is in `study.json` (`source.xml_sha256`).

**Article licence (quoted from the JATS `<permissions>` element):**

> © The Author(s) 20262026https://creativecommons.org/licenses/by-nc-nd/4.0/Open Access This article is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License, which permits any non-commercial use, sharing, distribution and reproduction in any medium or format, as long as you give appropriate credit to the original author(s) and the source, provide a link to the Creative Commons licence, and indicate if you modified the licensed material. You do not have permission under this licence to share adapted material derived from this article or parts of it. The images or other third party material in this article are included in the article’s Creative Commons licence, unless indicated otherwise in a credit line to the material. If material is not included in the article’s Creative Commons licence and your intended use is not permitted by statutory regulation or exceeds the permitted use, you will need to obtain permission directly from the copyright holder. To view a copy of this licence, visit http://creativecommons.org/licenses/by-nc-nd/4.0/.

**What is stored here.** The article is licensed CC BY-NC-ND 4.0, which does not permit sharing adapted material. To stay on the safe side this folder stores only the drug names, the structure resolution, the table layout (table and column identifiers), the index formulas and the reported fit statistics needed to cite the published result. The per-drug property values are **not** committed. `python tools/extract_replication.py fetch-xml <study>` followed by `extract <study>` rebuilds `properties.csv` in the local cache (`~/.cache/topoqspr/replication/`) from the Europe PMC XML.

**Our additions** (not from the article): PubChem CIDs, SMILES and InChIKeys (PubChem PUG-REST, public domain), resolution rules and notes, consistency checks.

**Attribution.** Data extracted from Supriya; Iyer, Radha Rajamani (2026), doi:10.1038/s41598-026-49132-8.
