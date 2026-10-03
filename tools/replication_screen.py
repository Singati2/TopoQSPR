"""Candidate list for Experiment D (replication of published disease-specific TI-QSPR studies).

The query is fixed here before any screening (EXPERIMENT_PLAN.md, Experiment D). Candidates are
listed newest first; eligibility (compound list resolvable to structures, property values
tabulated, index formulas stated, n >= 10) is screened in that order, without regard to the
R^2 a study reports.

Usage: python tools/replication_screen.py OUT.csv
"""

import csv
import sys

import requests

QUERY = (
    '("topological indices" OR "topological index") AND (QSPR OR "structure-property") '
    'AND (drug OR drugs) AND (regression) AND PUB_YEAR:[2020 TO 2026] AND OPEN_ACCESS:y'
)
URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def fetch(query: str = QUERY, page_size: int = 1000) -> list[dict]:
    rows, cursor = [], "*"
    while True:
        r = requests.get(URL, params={"query": query, "format": "json", "resultType": "lite",
                                      "pageSize": page_size, "cursorMark": cursor}, timeout=60)
        r.raise_for_status()
        d = r.json()
        rows += d["resultList"]["result"]
        nxt = d.get("nextCursorMark")
        if not nxt or nxt == cursor or not d["resultList"]["result"]:
            return rows
        cursor = nxt


def main(out: str) -> None:
    rows = fetch()
    rows.sort(key=lambda x: (x.get("firstPublicationDate") or "", x.get("id") or ""), reverse=True)
    cols = ["firstPublicationDate", "title", "authorString", "journalTitle", "doi", "pmcid"]
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols + ["eligible", "reason"])
        for x in rows:
            w.writerow([x.get(c, "") for c in cols] + ["", ""])
    print(f"{len(rows)} candidates written to {out}")


if __name__ == "__main__":
    main(sys.argv[1])
