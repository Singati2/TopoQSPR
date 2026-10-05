"""Convert an extracted replication study folder into the Experiment D study spec.

Usage: python tools/replication_to_spec.py data/replication/<study_key> OUT.json

Reads study.json (verbatim extraction), drugs.csv (resolved structures) and properties.csv,
and writes the spec that ``topoqspr.experiments.experiment_d`` accepts. Every index label is
normalised (``ND1(G)`` -> ``ND_1``, ``ABC^d`` -> ``ABC``, ``ReZ_{1}`` -> ``ReZ1``) and mapped
through ``map_index_name``; a label that does not map is reported, and the study is refused
if any *modelled* predictor does not map, so coverage gaps are explicit rather than silent.
Model forms the harness cannot reproduce (ridge, logarithmic, ...) are kept as reported fits
flagged ``reproducible: false``; they are listed, never compared.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from topoqspr.experiments import _FORMS, map_index_name  # noqa: E402

SUPPORTED_FORMS = set(_FORMS)


def normalise_label(label: str) -> str:
    """Strip typographic decoration from a paper's index symbol."""
    s = str(label).strip()
    s = re.sub(r"\(G\)$", "", s)            # ND1(G) -> ND1
    s = re.sub(r"\^\{?(d|nd)\}?$", "", s)   # ABC^d, ABC^{nd} -> ABC   (neighbourhood handled below)
    s = s.replace("_{", "").replace("}", "").replace("*", "")
    s = s.replace("^{", "").replace("^", "")
    m = re.match(r"^(ND|NM|NF)(\d)$", s)    # ND1 -> ND_1 (the registry's literature names)
    if m:
        s = f"{m.group(1)}_{m.group(2)}"
    return s


def is_neighbourhood(label: str) -> bool:
    return bool(re.search(r"\^\{?nd\}?$", str(label).strip()))


def map_label(label: str) -> tuple[str | None, str]:
    norm = normalise_label(label)
    if is_neighbourhood(label):
        # Godlin-style X^nd: the neighbourhood-sum version of a degree index
        name, how = map_index_name("S_" + norm)
        if name:
            return name, how
    return map_index_name(norm)


def convert(folder: Path) -> tuple[dict, dict]:
    study = json.loads((folder / "study.json").read_text(encoding="utf-8"))
    drugs = list(csv.DictReader((folder / "drugs.csv").open(encoding="utf-8")))
    props = list(csv.DictReader((folder / "properties.csv").open(encoding="utf-8")))

    values: dict[str, dict[str, float]] = {}
    for r in props:
        if r.get("value", "") in ("", None):
            continue
        values.setdefault(r["drug_id"], {})[r["property"]] = float(r["value"])

    compounds, skipped = [], []
    for d in drugs:
        if d.get("status") != "resolved" or not d.get("isomeric_smiles"):
            skipped.append(dict(name=d.get("name_verbatim"), reason=d.get("status") or "unresolved"))
            continue
        compounds.append(dict(name=d["name_verbatim"].strip(), smiles=d["isomeric_smiles"],
                              pubchem_cid=d.get("pubchem_cid"), properties=values.get(d["drug_id"], {})))

    index_map, unmapped = {}, []
    for i in study["indices"]:
        lab = i["symbol_verbatim"]
        name, how = map_label(lab)
        if name:
            index_map[lab] = name
        else:
            unmapped.append(dict(label=lab, name=i.get("name_verbatim", ""), reason=how))

    reported_fits, unsupported = [], []
    for m in study["models"]:
        if len(m["predictors"]) != 1:
            unsupported.append(dict(id=m["id"], reason=f"{len(m['predictors'])} predictors"))
            continue
        lab = m["predictors"][0]
        stats = m.get("statistics") or {}
        r2 = stats.get("R2") if stats.get("R2") not in (None, "") else None
        rr = stats.get("R") if stats.get("R") not in (None, "") else None
        fit = dict(id=m["id"], property=m["property"], index=index_map.get(lab, normalise_label(lab)),
                   study_label=lab, form=m["form"], source=m.get("source", ""))
        if r2 is not None:
            fit["r2"] = float(r2)
        elif rr is not None:
            fit["r"] = float(rr)
        if m["form"] not in SUPPORTED_FORMS or lab not in index_map:
            fit["reproducible"] = False
            unsupported.append(dict(id=m["id"], reason=("form " + m["form"]) if m["form"] not in SUPPORTED_FORMS
                                    else f"predictor {lab} unmapped"))
        reported_fits.append(fit)

    modelled = {m["predictors"][0] for m in study["models"] if len(m["predictors"]) == 1}
    blocking = sorted(lab for lab in modelled if lab not in index_map)
    forms = sorted({f["form"] for f in reported_fits if f["form"] in SUPPORTED_FORMS}) or ["linear"]

    spec = dict(
        name=study["study_key"],
        reference=(study.get("citation") or {}).get("text") if isinstance(study.get("citation"), dict) else str(study.get("citation", "")),
        notes=f"converted by tools/replication_to_spec.py from data/replication/{folder.name}; "
              f"licence {study.get('licence', '')}",
        compounds=compounds,
        properties=sorted({p for c in compounds for p in c["properties"]}),
        indices=sorted(set(index_map.values())),
        model_forms=forms,
        reported_fits=reported_fits,
    )
    report = dict(study=study["study_key"], n_compounds=len(compounds), skipped_compounds=skipped,
                  index_map=index_map, unmapped_indices=unmapped, blocking_unmapped_predictors=blocking,
                  unsupported_fits=unsupported, n_reported_fits=len(reported_fits),
                  n_comparable_fits=sum(1 for f in reported_fits if f.get("reproducible", True)))
    return spec, report


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("out")
    ap.add_argument("--allow-blocking", action="store_true",
                    help="write the spec even if a modelled predictor does not map (its fits are not compared)")
    a = ap.parse_args()
    spec, report = convert(Path(a.folder))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).with_suffix(".report.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    if report["blocking_unmapped_predictors"] and not a.allow_blocking:
        print(f"refused: modelled predictors without a registry index: {report['blocking_unmapped_predictors']}",
              file=sys.stderr)
        return 2
    Path(a.out).write_text(json.dumps(spec, indent=1) + "\n", encoding="utf-8")
    print(f"{report['study']}: {report['n_compounds']} compounds, {len(report['index_map'])} indices mapped, "
          f"{report['n_comparable_fits']}/{report['n_reported_fits']} fits comparable; wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
