"""Regenerate the verbatim CC0 excerpts used by tests/test_expdata.py.

Both sources are CC0 (AqSolDB: Harvard Dataverse doi:10.7910/DVN/OVHAW8; Bradley:
figshare doi:10.6084/m9.figshare.1031637.v2), so storing an excerpt is permitted.
Run once with the pinned files in the user cache:

    PYTHONPATH=src python tests/fixtures/expdata/make_excerpts.py
"""

from pathlib import Path

import pandas as pd

from topoqspr.expdata import fetch

HERE = Path(__file__).resolve().parent

aq = pd.read_csv(fetch("aqsoldb", offline=True), dtype={"ID": str})
parts = [
    aq.head(20),
    aq[aq["Group"].isin(["G2", "G4"])].head(4),  # source SD > 0.5
    aq[aq["Group"].isin(["G3", "G5"])].head(3),
    aq[aq["SMILES"].str.contains(".", regex=False)].head(5),  # salts / hydrates
]
pd.concat(parts).drop_duplicates("ID").sort_index().to_csv(HERE / "aqsoldb_excerpt.csv", index=False)

br = pd.read_excel(fetch("bradley", offline=True), sheet_name="Sheet1", engine="openpyxl")
good = br[br["donotuse"].isna() & (br["csid"] > 0)]
multi = good[good.duplicated("csid", keep=False)].sort_values(["csid", "key"])
csids = multi["csid"].drop_duplicates().head(6)
parts = [br.head(10), multi[multi["csid"].isin(csids)], br[br["donotuse"].notna()].head(4)]
pd.concat(parts).drop_duplicates("key").sort_index().to_csv(HERE / "bradley_excerpt.csv", index=False)
