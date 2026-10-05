"""Convert an extracted replication study folder into the Experiment D study spec.

Usage: python tools/replication_to_spec.py data/replication/<study_key> OUT.json

Reads study.json (verbatim extraction), drugs.csv (resolved structures), properties.csv and,
when present, index_values.csv (the study's printed index values, passed on as
``reported_index_values`` so the harness checks them against its own), and writes the spec that ``topoqspr.experiments.experiment_d`` accepts. Every index label is
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
        # Godlin-style X^nd is the neighbourhood-degree-sum version of a degree index; it maps
        # to the registry's S_X and to nothing else (falling back to the degree index X would
        # silently compare the paper's neighbourhood index with a different one).
        name, how = map_index_name("S_" + norm)
        return (name, how) if name else (None, f"no neighbourhood-sum index S_{norm} in the registry")
    return map_index_name(norm)


def _stat(stats: dict, *keys: str) -> float | None:
    """A printed statistic by case-insensitive key (papers write R2, r2, R^2 ...); None if absent."""
    low = {str(k).strip().lower(): v for k, v in stats.items()}
    for k in keys:
        v = low.get(k)
        if v not in (None, ""):
            try:
                return float(str(v).strip())
            except ValueError:
                return None
    return None


def convert(folder: Path) -> tuple[dict, dict]:
    study = json.loads((folder / "study.json").read_text(encoding="utf-8"))
    drugs = list(csv.DictReader((folder / "drugs.csv").open(encoding="utf-8")))
    props = list(csv.DictReader((folder / "properties.csv").open(encoding="utf-8")))

    values: dict[str, dict[str, float]] = {}
    for r in props:
        if r.get("value", "") in ("", None):
            continue
        values.setdefault(r["drug_id"], {})[r["property"]] = float(r["value"])

    printed: dict[str, dict[str, float]] = {}
    if (folder / "index_values.csv").exists():
        for r in csv.DictReader((folder / "index_values.csv").open(encoding="utf-8")):
            if r.get("value", "") not in ("", None):
                printed.setdefault(r["drug_id"], {})[r["index_label"]] = float(r["value"])

    index_map, unmapped = {}, []
    for i in study["indices"]:
        lab = i["symbol_verbatim"]
        name, how = map_label(lab)
        if name:
            index_map[lab] = name
        else:
            unmapped.append(dict(label=lab, name=i.get("name_verbatim", ""), reason=how))

    compounds, skipped = [], []
    reported_index_values: dict[str, dict[str, float]] = {}
    for d in drugs:
        if d.get("status") != "resolved" or not d.get("isomeric_smiles"):
            skipped.append(dict(name=d.get("name_verbatim"), reason=d.get("status") or "unresolved"))
            continue
        compounds.append(dict(name=d["name_verbatim"].strip(), smiles=d["isomeric_smiles"],
                              pubchem_cid=d.get("pubchem_cid"), properties=values.get(d["drug_id"], {})))
        if printed.get(d["drug_id"]):
            # keyed by registry name when mapped, else by the paper's own label (never by the
            # normalised label: M1^nd would then overwrite the mapped M1^d -> M1 entry)
            reported_index_values[d["name_verbatim"].strip()] = {
                index_map.get(lab, lab): v for lab, v in printed[d["drug_id"]].items()}

    reported_fits, unsupported, reported_loocv, reported_other = [], [], [], []
    for m in study["models"]:
        if len(m["predictors"]) != 1:
            unsupported.append(dict(id=m["id"], reason=f"{len(m['predictors'])} predictors"))
            continue
        lab = m["predictors"][0]
        stats = m.get("statistics") or {}
        fit_type = str(m.get("fit", "in_sample"))
        if fit_type == "loocv":
            # A printed LOO statistic is not an in-sample fit: it is kept in the report and set
            # beside the harness's own LOO Q^2 of the same (property, index, form) by make_expD.
            reported_loocv.append(dict(id=m["id"], property=m["property"],
                                       index=index_map.get(lab, normalise_label(lab)), study_label=lab,
                                       form=m["form"], q2_loo=_stat(stats, "q2_loo", "q2", "q2loo"),
                                       r2=_stat(stats, "r2", "r^2", "r²"), source=m.get("source", "")))
            continue
        if fit_type != "in_sample":
            reported_other.append(dict(id=m["id"], fit=fit_type, property=m["property"], study_label=lab))
            continue
        r2, rr = _stat(stats, "r2", "r^2", "r²"), _stat(stats, "r")
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
        reported_index_values=reported_index_values,
    )
    report = dict(study=study["study_key"], n_compounds=len(compounds), skipped_compounds=skipped,
                  n_reported_index_values=sum(len(v) for v in reported_index_values.values()),
                  index_map=index_map, unmapped_indices=unmapped, blocking_unmapped_predictors=blocking,
                  unsupported_fits=unsupported, n_reported_fits=len(reported_fits),
                  reported_loocv=reported_loocv, reported_other_validation=reported_other,
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
