"""Regenerate the fixtures used by tests/test_indices.py.

    PYTHONPATH=src python tests/fixtures/indices/make_fixtures.py [--fetch]

``drugs.json``  - drug panel: name, PubChem CID and PubChem isomeric SMILES.
                  Only rewritten with ``--fetch`` (needs network: PubChem
                  PUG REST, ``/compound/name/<name>/property/IsomericSMILES``).
``mordred_reference.json`` - Mordred-community values of the descriptors
                  that overlap TopoQSPR indices, computed on the *standardised
                  parent* molecule (so both tools see the same graph).

Mordred descriptor -> TopoQSPR index:
    ABC      -> ABC   (atom-bond connectivity)
    BalabanJ -> J
    WPath    -> W     (Wiener path number)
    Zagreb1  -> M1    (ZagrebIndex(lambda=1, variable=1))
    Zagreb2  -> M2    (ZagrebIndex(lambda=2, variable=1))
    Xp-1d    -> R     (path chi index of order 1 on simple degrees = Randic)
"""

from __future__ import annotations

import argparse
import importlib.metadata as md
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent

DRUG_NAMES = [
    "aspirin", "ibuprofen", "acetaminophen", "caffeine", "metformin", "imatinib",
    "atorvastatin", "donepezil", "memantine", "sitagliptin", "hydrochlorothiazide",
    "gemcitabine", "diazepam", "warfarin", "lisinopril", "amlodipine",
    "ciprofloxacin", "morphine", "nicotine", "sildenafil", "vancomycin",
    "cyclosporine", "rifampicin", "paclitaxel",
]

MORDRED_TO_TOPOQSPR = {
    "ABC": "ABC",
    "BalabanJ": "J",
    "WPath": "W",
    "Zagreb1": "M1",
    "Zagreb2": "M2",
    "Xp-1d": "R",
}


def fetch_drugs() -> list[dict]:
    out = []
    for name in DRUG_NAMES:
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            + urllib.parse.quote(name)
            + "/property/IsomericSMILES/JSON"
        )
        with urllib.request.urlopen(url, timeout=30) as fh:
            prop = json.load(fh)["PropertyTable"]["Properties"][0]
        smiles = prop.get("IsomericSMILES") or prop.get("SMILES")
        out.append({"name": name, "pubchem_cid": prop["CID"], "smiles": smiles})
        time.sleep(0.25)
    return out


def mordred_reference(drugs: list[dict]) -> dict:
    from mordred import ABCIndex, BalabanJ, Calculator, Chi, WienerIndex, ZagrebIndex

    from topoqspr.graph import standardize

    calc = Calculator(
        [
            ABCIndex.ABCIndex(),
            BalabanJ.BalabanJ(),
            WienerIndex.WienerIndex(False),
            ZagrebIndex.ZagrebIndex(1, 1),
            ZagrebIndex.ZagrebIndex(2, 1),
            Chi.Chi("path", 1, "d", False),
        ]
    )
    names = [str(d) for d in calc.descriptors]
    assert names == list(MORDRED_TO_TOPOQSPR), names
    values = {}
    for d in drugs:
        res = standardize(d["smiles"])
        assert res.ok, (d["name"], res.rejection)
        row = calc(res.mol)
        values[d["name"]] = {
            "canonical_smiles": res.canonical_smiles,
            "inchikey": res.inchikey,
            "mordred": {k: float(v) for k, v in zip(names, row)},
        }
    return {
        "generated_with": {
            "mordredcommunity": md.version("mordredcommunity"),
            "rdkit": md.version("rdkit"),
            "chembl_structure_pipeline": md.version("chembl_structure_pipeline"),
        },
        "mapping": MORDRED_TO_TOPOQSPR,
        "values": values,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true", help="re-download drugs.json from PubChem")
    args = ap.parse_args()
    drugs_path = HERE / "drugs.json"
    if args.fetch:
        drugs = fetch_drugs()
        payload = {
            "source": "PubChem PUG REST, property IsomericSMILES",
            "retrieved": time.strftime("%Y-%m-%d"),
            "drugs": drugs,
        }
        drugs_path.write_text(json.dumps(payload, indent=1) + "\n")
    drugs = json.loads(drugs_path.read_text())["drugs"]
    ref = mordred_reference(drugs)
    (HERE / "mordred_reference.json").write_text(json.dumps(ref, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
