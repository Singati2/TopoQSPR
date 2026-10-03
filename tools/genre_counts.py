"""Europe PMC hit counts for drug topological-index QSPR papers per year (paper, Section 1).

Usage: python tools/genre_counts.py [first_year] [last_year]
"""
import sys

import requests

BASE = '("topological indices" OR "topological index") AND QSPR AND (drug OR drugs)'
ML = ' AND ("machine learning" OR "random forest" OR XGBoost)'
URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def hits(query: str) -> int:
    r = requests.get(URL, params={"query": query, "format": "json", "pageSize": 1}, timeout=60)
    r.raise_for_status()
    return int(r.json()["hitCount"])


def main(first: int = 2020, last: int = 2026) -> None:
    print("year\tall\twith_ML")
    for year in range(first, last + 1):
        q = f"{BASE} AND PUB_YEAR:{year}"
        print(f"{year}\t{hits(q)}\t{hits(q + ML)}")


if __name__ == "__main__":
    main(*(int(a) for a in sys.argv[1:3]))
