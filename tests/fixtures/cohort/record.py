"""Re-record the offline fixtures used by tests/test_cohort.py (needs network).

    PYTHONPATH=src python tests/fixtures/cohort/record.py

The fixture directories are ordinary topoqspr caches: ``chembl/`` (ChEMBL REST, keyed by
release) and ``mesh/`` (NLM MeSH).  Re-recording against a newer ChEMBL release creates a
new ``chembl/<release>/`` directory; delete the old one and update the expected manifest
(``expected_gout.csv``) in the same commit.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from topoqspr.chembl import ChemblClient, MeshClient
from topoqspr.cohort import build_atc_cohort, build_cohort, naive_cohort, suggest_headings

HERE = Path(__file__).resolve().parent
GOUT, ARTHRITIS_GOUTY, T2D = "D006073", "D015210", "D003924"


def record(fresh: bool = False) -> None:
    if fresh:
        shutil.rmtree(HERE / "chembl", ignore_errors=True)
        shutil.rmtree(HERE / "mesh", ignore_errors=True)
    client = ChemblClient(HERE / "chembl")
    mesh = MeshClient(HERE / "mesh")
    client.mesh_headings()
    suggest_headings("gout", client, mesh_client=mesh)
    suggest_headings("breast cancer", client, mesh_client=mesh)
    gout = build_cohort([GOUT], client)
    build_cohort([GOUT], client, include_descendants=True, mesh_client=mesh)
    build_cohort([T2D], client)
    naive_cohort("gout", client)
    build_atc_cohort("M04", client)
    (HERE / "expected_gout.csv").write_bytes(gout.to_csv_bytes())


if __name__ == "__main__":
    record(fresh=True)
