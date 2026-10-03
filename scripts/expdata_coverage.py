"""Drug coverage of the measured-property targets (counts only, no modelling).

Usage::

    PYTHONPATH=src python scripts/expdata_coverage.py [--offline] [--out docs/expdata_coverage.json]

Downloads (or reads from the user cache) every source in ``topoqspr.expdata.TARGETS``,
builds the pooled approved-drug cohort (``approved_drug_pool``) and the controlled
phase-4 cohort of every disease on the EXPERIMENT_PLAN candidate list plus the
quarantined pilot, standardises their parents and counts, per target, how many have a
measured value.  Writes a JSON record with the ChEMBL release, source hashes and counts.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from topoqspr import __version__
from topoqspr.chembl import ChemblClient
from topoqspr.cohort import build_cohort
from topoqspr.expdata import (
    SOURCES,
    TARGETS,
    aggregate,
    approved_drug_pool,
    coverage,
    get_target,
    load_records,
    standardised_keys,
)

#: EXPERIMENT_PLAN.md section 3.3 candidate list (exact MeSH headings), then the pilot.
CANDIDATES = (
    "Diabetes Mellitus, Type 2",
    "Hypertension",
    "Asthma",
    "Breast Neoplasms",
    "Epilepsy",
    "Depressive Disorder, Major",
    "Schizophrenia",
    "Arthritis, Rheumatoid",
    "Hypercholesterolemia",
    "Heart Failure",
    "Pulmonary Disease, Chronic Obstructive",
    "HIV Infections",
    "Parkinson Disease",
    "Alzheimer Disease",
    "Prostatic Neoplasms",
    "Psoriasis",
)
PILOT = "Migraine Disorders"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--out", default="docs/expdata_coverage.json")
    args = ap.parse_args(argv)

    client = ChemblClient(offline=args.offline)
    result: dict = {"topoqspr_version": __version__, "chembl_release": client.release, "targets": {}}

    tgt_default, tgt_sens = {}, {}
    for name, spec in TARGETS.items():
        recs = load_records(name, offline=args.offline, client=client)
        agg = aggregate(recs, spec)
        tgt_default[name] = agg[~agg["replicates_disagree"]]
        salts = (
            get_target(name, offline=args.offline, client=client, keep_salt_forms=True, drop_disagreeing=False)
            if spec.solid_state
            else agg
        )
        tgt_sens[name] = salts
        result["targets"][name] = {
            "source": spec.source,
            "unit": spec.unit,
            "raw_hash": recs.attrs["raw_sha256"],
            "records": int(len(recs)),
            "status_counts": {k: int(v) for k, v in sorted(recs["status"].value_counts().items())},
            "parents": int(len(agg)),
            "parents_disagreeing": int(agg["replicates_disagree"].sum()),
            "parents_used": int(len(tgt_default[name])),
            "parents_incl_salt_forms_and_disagreeing": int(len(salts)),
        }
        print(name, result["targets"][name], file=sys.stderr)

    pool = approved_drug_pool(client)
    keys = standardised_keys(pool.included)
    ok = keys[keys["status"] == "ok"]
    result["pooled"] = {
        "definition": "drug_indication max_phase_for_ind=4, any MeSH heading; build_cohort curation",
        "indication_rows": pool.extra["indication_rows"],
        "mesh_ids": pool.extra["n_mesh_ids"],
        "curated_parents": int(len(pool.included)),
        "standardised_ok": int(len(ok)),
        "standardised_unique_inchikeys": int(ok["inchikey"].nunique()),
        "standardisation_rejections": {k: int(v) for k, v in sorted(keys["status"].value_counts().items()) if k != "ok"},
        "coverage": coverage(tgt_default, ok["inchikey"]),
        "coverage_incl_salt_forms_and_disagreeing": coverage(tgt_sens, ok["inchikey"]),
    }
    print("pooled", result["pooled"], file=sys.stderr)

    heading_to_id = {h["mesh_heading"]: h["mesh_id"] for h in client.mesh_headings() if h["mesh_id"]}
    result["diseases"] = {}
    for heading in (*CANDIDATES, PILOT):
        mid = heading_to_id.get(heading)
        if mid is None:
            result["diseases"][heading] = {"mesh_id": None, "note": "heading not in drug_indication"}
            continue
        coh = build_cohort([mid], client)
        k = standardised_keys(coh.included)
        k = k[k["status"] == "ok"]
        result["diseases"][heading] = {
            "mesh_id": mid,
            "pilot": heading == PILOT,
            "curated_parents": int(len(coh.included)),
            "standardised_unique_inchikeys": int(k["inchikey"].nunique()),
            "coverage": {t: c["n"] for t, c in coverage(tgt_default, k["inchikey"]).items()},
        }
        print(heading, result["diseases"][heading], file=sys.stderr)
    client.verify_release()
    result["sources"] = {
        k: {"doi": s.doi, "version": s.version, "licence": s.licence, "sha256": s.sha256 or dict(s.chembl_sha256 or {})}
        for k, s in SOURCES.items()
    }
    Path(args.out).write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
