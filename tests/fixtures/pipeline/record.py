"""Re-record the offline ChEMBL cache used by tests/test_pipeline.py (needs network).

    PYTHONPATH=src python tests/fixtures/pipeline/record.py

The cache holds the quarantined pilot cohort only (EXPERIMENT_PLAN §2: Migraine
Disorders, MeSH D008881) plus the release's MeSH heading index. It is an
ordinary topoqspr ChEMBL cache keyed by release; re-recording against a newer
release creates a new ``chembl/<release>/`` directory (delete the old one in
the same commit). ChEMBL data: CC BY-SA 3.0, EMBL-EBI.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from topoqspr.chembl import ChemblClient
from topoqspr.cohort import build_cohort

HERE = Path(__file__).resolve().parent
MIGRAINE = "D008881"


def record(fresh: bool = False) -> None:
    if fresh:
        shutil.rmtree(HERE / "chembl", ignore_errors=True)
    client = ChemblClient(HERE / "chembl")
    client.mesh_headings()
    build_cohort([MIGRAINE], client)  # promotes the staged entries (verify_release)


if __name__ == "__main__":
    record(fresh=True)
