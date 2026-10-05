"""Experiment D replication data: extraction, structure resolution and validation.

Each study lives in ``data/replication/<study_key>/``:

* ``study.json``      citation, licence, what the paper reports (property table layout, index
                      formulas, model forms, fit statistics, verbatim), consistency notes and a
                      reproducibility verdict.  The block ``automated_checks`` is written by
                      this tool and must match a recomputation.
* ``drugs.csv``       the drug list exactly as printed plus the resolution *inputs*
                      (``query_name``, ``rule``) and *outputs* (PubChem CID, SMILES, InChIKey,
                      status) written by this tool.
* ``properties.csv``  long format: one row per (drug, property) with the printed string and its
                      parsed value (NaN + reason when missing).  Absent when the article licence
                      does not allow redistribution (``licence.data_redistribution ==
                      "withheld"``); ``extract`` then writes it to a local cache instead.
* ``SOURCE.md``       attribution and licence.

Structure resolution is deterministic: every PubChem PUG-REST response is cached in
``data/replication/_pubchem_cache.json`` (keyed by request path), so ``resolve --offline``
reproduces ``drugs.csv`` byte for byte.  Resolution rules (docs/REPLICATION_DATA.md §2):

* ``as_printed``        query the printed name; one CID -> ``resolved``.
* ``typo_corrected``    the printed name is misspelt; ``query_name`` holds the correction and
                        ``note`` the evidence (e.g. the paper's own formula or MW).
* ``spelling_variant``  INN/USAN or British/American variant (e.g. azacytidine/azacitidine).
* ``code_mapped``       the table prints a code (BD1, A1, G1); the name comes from the text.
* ``parent_form``       the printed name denotes an ion or salt while the paper's own formula or
                        values describe the neutral parent; ``query_name`` is the parent.
* ``component_by_printed_index``  a combination product whose printed index values identify
                        one component; ``query_name`` is that component and ``note`` the evidence
                        (a printed index value recomputed exactly for it and not for the other).
* ``unresolvable``      combination product, polymer or mixture: no single structure, status
                        ``unresolved`` (never guessed).
If PubChem returns several CIDs for a name the first (PubChem's ranking) is taken and the row is
marked ``ambiguous_first_cid``; a multi-component record is marked ``resolved_multicomponent``.

Usage::

    python tools/extract_replication.py validate            # offline schema + consistency
    python tools/extract_replication.py resolve [--offline] [study ...]
    python tools/extract_replication.py fetch-xml [study ...]          # network, cached
    python tools/extract_replication.py extract [--record] [study ...]  # needs cached XML
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
from decimal import ROUND_DOWN, Decimal
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Callable

REPO = Path(__file__).resolve().parents[1]
DATA_DIR = REPO / "data" / "replication"
CACHE_PATH = DATA_DIR / "_pubchem_cache.json"
XML_CACHE = Path.home() / ".cache" / "topoqspr" / "replication"
SCHEMA_VERSION = 1

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_FIELDS = ("SMILES", "ConnectivitySMILES", "InChIKey", "MolecularFormula",
                  "MolecularWeight", "ExactMass", "HeavyAtomCount", "Complexity", "TPSA",
                  "XLogP", "Title")
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"

DRUG_INPUT_COLUMNS = ["drug_id", "row_label_verbatim", "name_verbatim", "name_source",
                      "query_name", "rule", "formula_verbatim", "cid_verbatim", "note"]
DRUG_OUTPUT_COLUMNS = ["status", "pubchem_cid", "n_cids", "pubchem_title", "isomeric_smiles",
                       "connectivity_smiles", "inchikey", "molecular_formula",
                       "molecular_weight", "exact_mass", "heavy_atom_count", "complexity", "tpsa",
                       "xlogp", "n_components"]
DRUG_COLUMNS = DRUG_INPUT_COLUMNS + DRUG_OUTPUT_COLUMNS
PROPERTY_COLUMNS = ["drug_id", "row_label_verbatim", "property", "value_verbatim", "value",
                    "missing_reason"]

RULES = {"as_printed", "typo_corrected", "spelling_variant", "code_mapped", "parent_form",
         "component_by_printed_index", "unresolvable"}
STATUSES = {"resolved", "ambiguous_first_cid", "resolved_multicomponent", "unresolved"}
SOURCE_LABELS = {"predicted_acd_chemspider", "computed_pubchem", "predicted_chemdraw",
                 "measured", "measured_or_predicted_unclear", "mixed", "unknown"}
# 'computed_pubchem' must be confirmed by automated_checks: at least this fraction of the
# compared values agree with PubChem (otherwise the column is labelled 'mixed' or 'unknown').
COMPUTED_PUBCHEM_MIN_MATCH = 0.8
MODEL_FORMS = {"linear", "quadratic", "cubic", "logarithmic", "exponential", "power",
               "multiple_linear", "ridge", "elastic_net", "random_forest", "xgboost", "ann",
               "fuzzy_ann"}
FIT_TYPES = {"in_sample", "loocv", "kfold_cv", "test_split", "external", "bootstrap"}
VERDICTS = {"exact", "approximate", "not_reproducible"}
STUDY_STATUS = {"proposed", "reserve"}
REQUIRED_STUDY_KEYS = {
    "schema_version": int, "study_key": str, "shortlist_row": int, "status": str,
    "citation": dict, "licence": dict, "source": dict, "n_reported": int,
    "n_reported_source": str, "drug_list_source": str, "property_table": dict,
    "properties": list, "indices": list, "index_values_table": (str, type(None)),
    "models": list, "validation_reported": list, "consistency": list,
    "resolution_issues": list, "reproducibility": dict, "automated_checks": dict,
}
# PubChem computed field used to label a property column as 'computed_pubchem'.
PUBCHEM_COMPARABLE = {"MolecularWeight": "molecular_weight", "HeavyAtomCount": "heavy_atom_count",
                      "Complexity": "complexity", "TPSA": "tpsa", "XLogP": "xlogp"}
# Native decimals of PubChem's computed values.  PUG-REST JSON drops trailing zeros (an XLogP
# of 3.0 arrives as 3, a TPSA of 56.0 as 56), so the reference is read at this precision.
# MolecularWeight is a string carrying its own precision; HeavyAtomCount is an integer.
# Complexity is special: PubChem computes it to 1 dp (PUG-View: valproic acid 93.4,
# phenelzine 77.3, ethanol 2.8) but PUG-REST returns it truncated to an integer (93, 77, 2).
PUBCHEM_NATIVE_DECIMALS = {"XLogP": 1, "TPSA": 1, "HeavyAtomCount": 0}
# An equation_verbatim must be a substring of the rendered XML; location or interpretation
# remarks belong in 'source' or 'note', never in the verbatim field.
_ANNOTATION = re.compile(r"\((?:Table|Figure|Fig\.|Eq\.|equation)[^)]*\)\s*$", re.I)


# --------------------------------------------------------------------------- JATS rendering

MML = "{http://www.w3.org/1998/Math/MathML}"
XLINK = "{http://www.w3.org/1999/xlink}href"


def _mathml(e: ET.Element) -> str:
    """Linearise a MathML element (fractions as (a)/(b), scripts as _{} and ^{})."""
    tag = e.tag.replace(MML, "")
    kids = list(e)
    t = (e.text or "").strip()
    if tag in ("mi", "mn", "mo", "mtext", "ms"):
        return t
    if tag == "mfrac" and len(kids) == 2:
        return f"({_mathml(kids[0])})/({_mathml(kids[1])})"
    if tag == "msqrt":
        return "sqrt(" + "".join(_mathml(k) for k in kids) + ")"
    if tag == "mroot" and len(kids) == 2:
        return f"root[{_mathml(kids[1])}]({_mathml(kids[0])})"
    if tag in ("msup", "mover") and len(kids) == 2:
        return f"{_mathml(kids[0])}^{{{_mathml(kids[1])}}}"
    if tag in ("msub", "munder") and len(kids) == 2:
        return f"{_mathml(kids[0])}_{{{_mathml(kids[1])}}}"
    if tag in ("msubsup", "munderover") and len(kids) == 3:
        return f"{_mathml(kids[0])}_{{{_mathml(kids[1])}}}^{{{_mathml(kids[2])}}}"
    if tag == "mspace":
        return " "
    if tag == "mfenced":
        sep = e.get("separators", ",")
        inner = (sep[:1] or "").join(_mathml(k) for k in kids)
        return e.get("open", "(") + inner + e.get("close", ")")
    return t + "".join(_mathml(k) for k in kids)


def _tex(s: str) -> str:
    m = re.search(r"\\begin\{document\}(.*)\\end\{document\}", s, re.S)
    s = m.group(1) if m else s
    return re.sub(r"\s+", " ", s.strip().strip("$").strip())


def render(e: ET.Element) -> str:
    """Text of a JATS element with formulas rendered inline as $...$ (TeX or linear MathML)."""
    out = [e.text or ""]
    for k in e:
        if k.tag == "alternatives":
            tex = k.find("tex-math")
            mm = [c for c in k if c.tag.startswith(MML)]
            out.append(f" ${_tex(tex.text or '')}$ " if tex is not None
                       else (f" ${_mathml(mm[0])}$ " if mm else render(k)))
        elif k.tag == "tex-math":
            out.append(f" ${_tex(k.text or '')}$ ")
        elif k.tag.startswith(MML):
            out.append(f" ${_mathml(k)}$ ")
        elif k.tag in ("inline-graphic", "graphic"):
            out.append(f" [IMG {k.get(XLINK)}] ")
        else:
            out.append(render(k))
        out.append(k.tail or "")
    s = "".join(out).replace("\u00ad", "").replace("\u200b", "")
    return re.sub(r"\s+", " ", s).strip()


def read_table(xml: bytes, table_id: str) -> dict:
    """Return {'label', 'caption', 'rows'} for the <table-wrap> with this id (rows = list of
    rendered cell strings, header rows included, in document order)."""
    root = ET.fromstring(xml)
    for tw in root.iter("table-wrap"):
        if tw.get("id") == table_id:
            lab, cap = tw.find("label"), tw.find("caption")
            rows = [[render(c) for c in tr if c.tag in ("td", "th")] for tr in tw.iter("tr")]
            return {"label": render(lab) if lab is not None else "",
                    "caption": render(cap) if cap is not None else "", "rows": rows}
    raise KeyError(f"table {table_id!r} not found")


# --------------------------------------------------------------------------- values

_MINUS = str.maketrans({"\u2212": "-", "\u2013": "-", "\u2014": "-", "\u2009": "", "\u00a0": "",
                        "\u202f": "", " ": ""})
_DASHES = {"-", "\u2013", "\u2014", "\u2212", "N/A", "NA", "n/a"}


def parse_value(s: str) -> tuple[float, str]:
    """(value, missing_reason).  Empty or dash cells give NaN with a reason, never 0."""
    raw = (s or "").strip()
    if raw == "":
        return math.nan, "empty_cell"
    if raw in _DASHES:
        return math.nan, "dash_in_table"
    t = raw.translate(_MINUS)
    try:
        return float(t), ""
    except ValueError:
        return math.nan, "unparseable"


def fmt_float(x: Any) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return repr(float(x))


def extract_properties(xml: bytes, spec: dict, drugs: list[dict]) -> list[dict]:
    """Build properties.csv rows from the JATS property table described by ``spec``.

    spec = {"xml_table_id": .., "skip_rows": k, "label_column": i,
            "columns": [{"property": "BP", "column": j}, ...]}
    Rows are matched to drugs.csv by ``row_label_verbatim`` (in table order)."""
    tab = read_table(xml, spec["xml_table_id"])
    body = tab["rows"][spec.get("skip_rows", 1):]
    body = [r for r in body if any(c.strip() for c in r)]
    if len(body) != len(drugs):
        raise ValueError(f"{spec['xml_table_id']}: {len(body)} body rows, {len(drugs)} drugs")
    out = []
    for row, drug in zip(body, drugs):
        label = row[spec["label_column"]]
        if label != drug["row_label_verbatim"]:
            raise ValueError(f"row label {label!r} != drugs.csv {drug['row_label_verbatim']!r}")
        for col in spec["columns"]:
            cell = row[col["column"]] if col["column"] < len(row) else ""
            v, why = parse_value(cell)
            out.append({"drug_id": drug["drug_id"], "row_label_verbatim": label,
                        "property": col["property"], "value_verbatim": cell,
                        "value": fmt_float(v), "missing_reason": why})
    return out


# --------------------------------------------------------------------------- PubChem cache

class CacheMissError(RuntimeError):
    pass


class PubChemCache:
    """Request-path keyed cache of PUG-REST responses: {path: {"status": int, "json": ...}}."""

    def __init__(self, path: Path = CACHE_PATH, offline: bool = True,
                 getter: Callable[[str], tuple[int, Any]] | None = None):
        self.path = Path(path)
        self.offline = offline
        self.getter = getter or _http_get
        self.data: dict = json.loads(self.path.read_text()) if self.path.exists() else {}
        self.dirty = False

    def get(self, rel: str) -> tuple[int, Any]:
        if rel in self.data:
            e = self.data[rel]
            return e["status"], e["json"]
        if self.offline:
            raise CacheMissError(rel)
        status, js = self.getter(PUG + rel)
        if status not in (200, 404):
            raise RuntimeError(f"PubChem {status} for {rel}")
        self.data[rel] = {"status": status, "json": js}
        self.dirty = True
        return status, js

    def save(self) -> None:
        if self.dirty:
            self.path.write_text(json.dumps(self.data, indent=1, sort_keys=True,
                                            ensure_ascii=False) + "\n")
            self.dirty = False


def _http_get(url: str) -> tuple[int, Any]:  # pragma: no cover - network
    import requests
    for attempt in range(5):
        r = requests.get(url, timeout=60)
        if r.status_code == 503:
            time.sleep(2 + 3 * attempt)
            continue
        time.sleep(0.25)  # PubChem asks for <= 5 requests / s
        try:
            js = r.json()
        except ValueError:
            js = None
        return r.status_code, js
    return 503, None


def name_path(name: str) -> str:
    return "/compound/name/" + urllib.parse.quote(name, safe="") + "/cids/JSON"


def props_path(cid: int) -> str:
    return f"/compound/cid/{cid}/property/{','.join(PUBCHEM_FIELDS)}/JSON"


def resolve_drug(row: dict, cache: PubChemCache) -> dict:
    """Return the output columns for one drugs.csv row (pure function of row inputs + cache)."""
    out = {c: "" for c in DRUG_OUTPUT_COLUMNS}
    if row["rule"] == "unresolvable":
        out["status"] = "unresolved"
        return out
    status, js = cache.get(name_path(row["query_name"]))
    cids = (js or {}).get("IdentifierList", {}).get("CID", []) if status == 200 else []
    if not cids:
        out["status"] = "unresolved"
        out["n_cids"] = "0"
        return out
    cid = int(cids[0])
    _, pj = cache.get(props_path(cid))
    p = pj["PropertyTable"]["Properties"][0]
    smi = p.get("SMILES", "")
    ncomp = len(smi.split(".")) if smi else 0
    out.update({
        "pubchem_cid": str(cid), "n_cids": str(len(cids)), "pubchem_title": p.get("Title", ""),
        "isomeric_smiles": smi, "connectivity_smiles": p.get("ConnectivitySMILES", ""),
        "inchikey": p.get("InChIKey", ""), "molecular_formula": p.get("MolecularFormula", ""),
        "molecular_weight": str(p.get("MolecularWeight", "")),
        "exact_mass": str(p.get("ExactMass", "")),
        "heavy_atom_count": str(p.get("HeavyAtomCount", "")),
        "complexity": str(p.get("Complexity", "")), "tpsa": str(p.get("TPSA", "")),
        "xlogp": str(p.get("XLogP", "")), "n_components": str(ncomp)})
    out["status"] = ("ambiguous_first_cid" if len(cids) > 1 else
                     "resolved_multicomponent" if ncomp > 1 else "resolved")
    return out


# --------------------------------------------------------------------------- CSV io

def read_csv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def csv_text(rows: list[dict], columns: list[str]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=columns, lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in columns})
    return buf.getvalue()


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.write_text(csv_text(rows, columns), encoding="utf-8")


def properties_sha256(rows: list[dict]) -> str:
    """SHA-256 of the canonical properties.csv text (as ``write_csv`` writes it).  ``extract
    --record`` stores it as ``source.properties_sha256``; ``validate`` then detects any later
    edit of a committed value offline, and ``extract`` detects drift for withheld studies."""
    return hashlib.sha256(csv_text(rows, PROPERTY_COLUMNS).encode("utf-8")).hexdigest()


def equation_problems(xml: bytes, study: dict) -> list[str]:
    """Model ids whose ``equation_verbatim`` is not a substring of any rendered paragraph,
    caption, title, formula or table cell of the JATS XML."""
    root = ET.fromstring(xml)
    texts = [render(e) for e in root.iter()
             if e.tag in ("p", "td", "th", "caption", "title", "label", "disp-formula")]
    return [m["id"] for m in study["models"] if m.get("equation_verbatim")
            and not any(m["equation_verbatim"] in t for t in texts)]


def study_dirs(keys: list[str] | None = None) -> list[Path]:
    """All study folders, or the requested ones; an unknown key raises ``KeyError``."""
    ds = sorted(p for p in DATA_DIR.iterdir() if p.is_dir() and (p / "study.json").exists())
    if keys:
        unknown = sorted(set(keys) - {d.name for d in ds})
        if unknown:
            raise KeyError(f"unknown study key(s): {', '.join(unknown)}")
        ds = [d for d in ds if d.name in keys]
    return ds


def load_study(d: Path) -> dict:
    return json.loads((d / "study.json").read_text(encoding="utf-8"))


def dump_study(d: Path, study: dict) -> None:
    (d / "study.json").write_text(json.dumps(study, indent=2, ensure_ascii=False) + "\n",
                                  encoding="utf-8")


def load_properties(d: Path, study: dict) -> list[dict] | None:
    """Committed properties only (None for licence-withheld studies), so that
    ``automated_checks`` never depends on a local cache."""
    if study["licence"].get("data_redistribution") == "withheld":
        return None
    p = d / "properties.csv"
    return read_csv(p) if p.exists() else None


# --------------------------------------------------------------------------- automated checks

def _decimals(s: str) -> int:
    t = (s or "").translate(_MINUS)
    return len(t.split(".")[1]) if "." in t else 0


def matches_printed(printed: str, ref: str, field: str | None = None) -> bool:
    """True if the printed number and the reference agree at the coarser of their precisions:
    the more precise one, rounded to the coarser number of decimals, gives the other (a tie at
    exactly half a unit is accepted under either rounding convention).  Numbers given at the
    same precision must be equal, so integers are compared exactly.

    ``field`` names the PubChem field of ``ref``; its native precision is used where PUG-REST
    dropped trailing zeros (``PUBCHEM_NATIVE_DECIMALS``).  For ``Complexity`` the reference is
    PUG-REST's truncated integer: a printed integer must equal it and a printed decimal must
    truncate to it.  (PubChem itself prints, e.g., sertraline MW as '306.2'.)"""
    _, why = parse_value(printed)
    _, why2 = parse_value(ref)
    if why or why2:
        return False
    p = Decimal(printed.strip().translate(_MINUS))
    r = Decimal(ref.strip().translate(_MINUS))
    dp = _decimals(printed)
    dr = max(_decimals(ref), PUBCHEM_NATIVE_DECIMALS.get(field or "", 0))
    if field == "Complexity" and dp > 0 and _decimals(ref) == 0:
        return p.to_integral_value(rounding=ROUND_DOWN) == r
    if dp == dr:
        return p == r
    fine, coarse = (p, r) if dp > dr else (r, p)
    return abs(fine - coarse) <= Decimal(5).scaleb(-(min(dp, dr) + 1))


def _hill_norm(f: str) -> str:
    f = re.sub(r"[\s$_{}]", "", f or "")
    return f


def automated_checks(study: dict, drugs: list[dict], props: list[dict] | None) -> dict:
    """Deterministic consistency checks recorded in study.json['automated_checks']."""
    res: dict[str, Any] = {}
    res["n_drugs_listed"] = len(drugs)
    res["n_matches_reported"] = len(drugs) == study["n_reported"]
    st: dict[str, int] = {}
    for r in drugs:
        st[r["status"]] = st.get(r["status"], 0) + 1
    res["resolution_status_counts"] = dict(sorted(st.items()))
    # duplicate structures and stereo/connectivity duplicates
    by_key: dict[str, list[str]] = {}
    by_block: dict[str, list[str]] = {}
    for r in drugs:
        if r["inchikey"]:
            by_key.setdefault(r["inchikey"], []).append(r["drug_id"])
            by_block.setdefault(r["inchikey"][:14], []).append(r["drug_id"])
    res["duplicate_structures"] = sorted(v for v in by_key.values() if len(v) > 1)
    res["shared_connectivity_groups"] = sorted(v for v in by_block.values() if len(v) > 1)
    # printed formula / CID vs PubChem
    fm = [r for r in drugs if r["formula_verbatim"] and r["molecular_formula"]]
    res["formula_checked"] = len(fm)
    res["formula_mismatch"] = sorted(r["drug_id"] for r in fm
                                     if _hill_norm(r["formula_verbatim"]) != r["molecular_formula"])
    cm = [r for r in drugs if r["cid_verbatim"] and r["pubchem_cid"]]
    res["cid_checked"] = len(cm)
    res["cid_mismatch"] = sorted(r["drug_id"] for r in cm
                                 if r["cid_verbatim"].replace(" ", "") != r["pubchem_cid"])
    if props is None:
        res["properties"] = "withheld (licence); run `extract` to rebuild locally"
        return res
    # per property counts and PubChem comparisons
    by_drug = {r["drug_id"]: r for r in drugs}
    pinfo = {p["code"]: p for p in study["properties"]}
    per: dict[str, dict] = {}
    for code in [p["code"] for p in study["properties"]]:
        rows = [r for r in props if r["property"] == code]
        vals = [r for r in rows if r["value"] != ""]
        miss = sorted(r["drug_id"] for r in rows if r["value"] == "")
        d = {"n_values": len(vals), "n_missing": len(miss), "missing": miss}
        fld = pinfo[code].get("pubchem_field")
        if fld:
            col = PUBCHEM_COMPARABLE[fld]
            cmp_ = [(r, by_drug[r["drug_id"]][col]) for r in vals
                    if by_drug[r["drug_id"]][col] not in ("", None)]
            ok = [r["drug_id"] for r, ref in cmp_
                  if matches_printed(r["value_verbatim"], ref, fld)]
            bad = sorted(r["drug_id"] for r, _ in cmp_ if r["drug_id"] not in ok)
            d["pubchem_field"] = fld
            d["pubchem_compared"] = len(cmp_)
            d["pubchem_match"] = len(ok)
            d["pubchem_mismatch"] = bad
            if fld == "MolecularWeight":
                # diagnostic: a printed 'molecular weight' that is the monoisotopic mass
                d["mismatch_equal_to_exact_mass"] = sorted(
                    r["drug_id"] for r, _ in cmp_
                    if r["drug_id"] in bad and by_drug[r["drug_id"]]["exact_mass"]
                    and matches_printed(r["value_verbatim"], f"{float(by_drug[r['drug_id']]['exact_mass']):.6f}"))
        per[code] = d
    res["properties"] = per
    # identical property vectors for different drugs (copy errors or stereo duplicates)
    vec: dict[tuple, list[str]] = {}
    for did in [r["drug_id"] for r in drugs]:
        v = tuple(r["value_verbatim"] for r in props if r["drug_id"] == did)
        vec.setdefault(v, []).append(did)
    res["identical_property_rows"] = sorted(v for v in vec.values() if len(v) > 1)
    return res


# --------------------------------------------------------------------------- validation

def validate_study(d: Path, cache: PubChemCache | None = None) -> list[str]:
    """Return a list of problems (empty = valid).  Offline; uses the committed cache."""
    errs: list[str] = []
    try:
        s = load_study(d)
    except Exception as exc:  # noqa: BLE001
        return [f"{d.name}: study.json unreadable: {exc}"]
    for k, typ in REQUIRED_STUDY_KEYS.items():
        if k not in s:
            errs.append(f"missing key {k}")
        elif not isinstance(s[k], typ):
            errs.append(f"{k} has type {type(s[k]).__name__}")
    if errs:
        return [f"{d.name}: {e}" for e in errs]
    if s["schema_version"] != SCHEMA_VERSION:
        errs.append("schema_version")
    if s["study_key"] != d.name:
        errs.append("study_key != folder name")
    if s["status"] not in STUDY_STATUS:
        errs.append(f"status {s['status']!r}")
    for k in ("authors", "title", "journal", "year", "doi", "pmcid"):
        if not s["citation"].get(k):
            errs.append(f"citation.{k} missing")
    lic = s["licence"]
    if lic.get("data_redistribution") not in ("stored", "withheld"):
        errs.append("licence.data_redistribution")
    if not lic.get("statement_verbatim") or not lic.get("spdx"):
        errs.append("licence statement/spdx missing")
    if lic.get("data_redistribution") == "stored" and not lic["spdx"].startswith("CC-BY-4.0"):
        errs.append("values stored for a non CC-BY article")
    for k in ("fulltext_url", "retrieved", "xml_sha256", "properties_sha256"):
        if not s["source"].get(k):
            errs.append(f"source.{k} missing")
    pt = s["property_table"]
    for k in ("xml_table_id", "label", "label_column", "columns"):
        if k not in pt:
            errs.append(f"property_table.{k} missing")
    codes = [p.get("code") for p in s["properties"]]
    if len(set(codes)) != len(codes):
        errs.append("duplicate property codes")
    if sorted(c["property"] for c in pt.get("columns", [])) != sorted(codes):
        errs.append("property_table.columns do not match properties")
    for p in s["properties"]:
        for k in ("code", "name_verbatim", "unit_verbatim", "source_stated_verbatim",
                  "source_label", "source_evidence"):
            if k not in p:
                errs.append(f"property {p.get('code')}: {k} missing")
        if p.get("source_label") not in SOURCE_LABELS:
            errs.append(f"property {p.get('code')}: source_label {p.get('source_label')!r}")
        if p.get("pubchem_field") and p["pubchem_field"] not in PUBCHEM_COMPARABLE:
            errs.append(f"property {p.get('code')}: pubchem_field")
    if not s["indices"]:
        errs.append("no indices")
    for i in s["indices"]:
        for k in ("symbol_verbatim", "name_verbatim", "formula_verbatim", "source"):
            if k not in i:
                errs.append(f"index {i.get('symbol_verbatim')}: {k} missing")
    if not s["models"]:
        errs.append("no models")
    for m in s["models"]:
        for k in ("id", "property", "form", "predictors", "equation_verbatim", "statistics",
                  "fit", "source"):
            if k not in m:
                errs.append(f"model {m.get('id')}: {k} missing")
        if m.get("form") not in MODEL_FORMS:
            errs.append(f"model {m.get('id')}: form {m.get('form')!r}")
        if m.get("fit") not in FIT_TYPES:
            errs.append(f"model {m.get('id')}: fit {m.get('fit')!r}")
        if not isinstance(m.get("statistics"), dict) or not all(
                isinstance(v, str) for v in m.get("statistics", {}).values()):
            errs.append(f"model {m.get('id')}: statistics must map to verbatim strings")
        if m.get("property") not in codes:
            errs.append(f"model {m.get('id')}: property {m.get('property')!r} not in properties")
        eq = m.get("equation_verbatim")
        if eq is not None and (not isinstance(eq, str) or _ANNOTATION.search(eq)):
            errs.append(f"model {m.get('id')}: equation_verbatim carries a non-verbatim "
                        "annotation (move it to source or note)")
        if "equation_interpreted" in m and not m.get("interpretation_source"):
            errs.append(f"model {m.get('id')}: equation_interpreted needs interpretation_source")
        symbols = {i.get("symbol_verbatim") for i in s["indices"]}
        unknown = [p for p in m.get("predictors", []) if p not in symbols]
        if unknown:
            errs.append(f"model {m.get('id')}: predictors {unknown} not among index symbols")
    ids = [m.get("id") for m in s["models"]]
    if len(set(ids)) != len(ids):
        errs.append("duplicate model ids")
    rv = s["reproducibility"]
    if rv.get("verdict") not in VERDICTS or not rv.get("reasons"):
        errs.append("reproducibility verdict/reasons")
    # drugs.csv
    dp = d / "drugs.csv"
    if not dp.exists():
        errs.append("drugs.csv missing")
        return [f"{d.name}: {e}" for e in errs]
    with open(dp, newline="", encoding="utf-8") as fh:
        header = next(csv.reader(fh))
    if header != DRUG_COLUMNS:
        errs.append("drugs.csv columns")
    drugs = read_csv(dp)
    if len({r["drug_id"] for r in drugs}) != len(drugs):
        errs.append("duplicate drug_id")
    for r in drugs:
        if r["rule"] not in RULES:
            errs.append(f"drug {r['drug_id']}: rule {r['rule']!r}")
        if r["status"] not in STATUSES:
            errs.append(f"drug {r['drug_id']}: status {r['status']!r}")
        if r["rule"] not in ("as_printed", "code_mapped") and not r["note"]:
            errs.append(f"drug {r['drug_id']}: rule {r['rule']} needs a note")
        if r["rule"] == "as_printed" and r["query_name"] != r["name_verbatim"]:
            errs.append(f"drug {r['drug_id']}: as_printed but query_name differs")
        if r["rule"] == "as_printed" and _ion_name_for_acid(r["name_verbatim"],
                                                            r["pubchem_title"]):
            errs.append(f"drug {r['drug_id']}: ion/salt name {r['name_verbatim']!r} resolved to "
                        f"{r['pubchem_title']!r}; use rule parent_form with a note")
        if r["status"] != "unresolved" and not (r["pubchem_cid"] and r["inchikey"]
                                                 and r["isomeric_smiles"]):
            errs.append(f"drug {r['drug_id']}: resolved without CID/SMILES/InChIKey")
    if cache is not None:
        for r in drugs:
            try:
                exp = resolve_drug(r, cache)
            except CacheMissError as exc:
                errs.append(f"drug {r['drug_id']}: cache miss {exc}")
                continue
            for c in DRUG_OUTPUT_COLUMNS:
                if exp[c] != r[c]:
                    errs.append(f"drug {r['drug_id']}: {c} {r[c]!r} != recomputed {exp[c]!r}")
    # properties.csv
    pp = d / "properties.csv"
    if lic.get("data_redistribution") == "withheld":
        if pp.exists():
            errs.append("properties.csv committed although licence is withheld")
        props = None
    else:
        if not pp.exists():
            errs.append("properties.csv missing")
            return [f"{d.name}: {e}" for e in errs]
        with open(pp, newline="", encoding="utf-8") as fh:
            if next(csv.reader(fh)) != PROPERTY_COLUMNS:
                errs.append("properties.csv columns")
        props = read_csv(pp)
        if len(props) != len(drugs) * len(codes):
            errs.append(f"properties.csv has {len(props)} rows, expected "
                        f"{len(drugs)} x {len(codes)}")
        # rows must belong to their drugs and come in drug order, then property-column order
        order = [(dr["drug_id"], dr["row_label_verbatim"], c["property"])
                 for dr in drugs for c in pt.get("columns", [])]
        got = [(r["drug_id"], r["row_label_verbatim"], r["property"]) for r in props]
        if got != order:
            bad = next((g for g, o in zip(got, order) if g != o), None)
            errs.append(f"properties.csv rows not in drugs.csv order/labels (first: {bad})")
        if properties_sha256(props) != s["source"].get("properties_sha256"):
            errs.append("properties.csv differs from source.properties_sha256 "
                        "(edited after extraction?)")
        seen = set()
        ids_ok = {r["drug_id"] for r in drugs}
        for r in props:
            key = (r["drug_id"], r["property"])
            if key in seen:
                errs.append(f"duplicate property row {key}")
            seen.add(key)
            if r["drug_id"] not in ids_ok:
                errs.append(f"unknown drug_id {r['drug_id']}")
            if r["property"] not in codes:
                errs.append(f"unknown property {r['property']}")
            v, why = parse_value(r["value_verbatim"])
            if fmt_float(v) != r["value"] or why != r["missing_reason"]:
                errs.append(f"{key}: value/missing_reason inconsistent with value_verbatim")
    if cache is not None and not any("cache miss" in e for e in errs):
        if automated_checks(s, drugs, props) != s["automated_checks"]:
            errs.append("automated_checks out of date (run resolve)")
    # 'computed_pubchem' must be confirmed by the automated comparison
    per = s["automated_checks"].get("properties")
    for p in s["properties"]:
        if p.get("source_label") != "computed_pubchem" or not isinstance(per, dict):
            continue
        c = per.get(p.get("code"), {})
        n, k = c.get("pubchem_compared", 0), c.get("pubchem_match", 0)
        if not n or k / n < COMPUTED_PUBCHEM_MIN_MATCH:
            errs.append(f"property {p.get('code')}: source_label computed_pubchem but only "
                        f"{k}/{n} values match PubChem (< {COMPUTED_PUBCHEM_MIN_MATCH:.0%})")
    return [f"{d.name}: {e}" for e in errs]


def _ion_name_for_acid(name: str, title: str) -> bool:
    """A printed anion name ('valproate', 'meclofenamate') whose PubChem record is the acid."""
    n, t = name.strip().lower(), title.strip().lower()
    return n.endswith("ate") and t.endswith(" acid") and n != t


# --------------------------------------------------------------------------- commands

def cmd_resolve(keys: list[str] | None, offline: bool) -> int:
    cache = PubChemCache(offline=offline)
    try:
        for d in study_dirs(keys):
            s = load_study(d)
            drugs = read_csv(d / "drugs.csv")
            for r in drugs:
                r.update(resolve_drug(r, cache))
            write_csv(d / "drugs.csv", drugs, DRUG_COLUMNS)
            s["automated_checks"] = automated_checks(s, drugs, load_properties(d, s))
            dump_study(d, s)
            print(f"{d.name}: {s['automated_checks']['resolution_status_counts']}")
    finally:
        cache.save()
    return 0


def _xml_path(s: dict) -> Path:
    return XML_CACHE / s["study_key"] / f"{s['citation']['pmcid']}.xml"


def cmd_fetch_xml(keys: list[str] | None) -> int:  # pragma: no cover - network
    import requests
    for d in study_dirs(keys):
        s = load_study(d)
        p = _xml_path(s)
        if not p.exists():
            r = requests.get(EPMC.format(pmcid=s["citation"]["pmcid"]), timeout=600)
            r.raise_for_status()
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(r.content)
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        flag = "ok" if h == s["source"]["xml_sha256"] else "SHA-256 DIFFERS from study.json"
        print(f"{d.name}: {p} {flag}")
    return 0


def cmd_extract(keys: list[str] | None, record: bool = False) -> int:
    """Re-extract properties from the cached JATS XML and compare with properties.csv
    (or write it to the local cache for licence-withheld studies), compare its SHA-256 with
    ``source.properties_sha256`` (``--record`` writes it) and check that every
    ``equation_verbatim`` occurs verbatim in the XML."""
    bad = 0
    for d in study_dirs(keys):
        s = load_study(d)
        p = _xml_path(s)
        if not p.exists():
            print(f"{d.name}: no cached XML ({p}); run fetch-xml")
            bad += 1
            continue
        xml = p.read_bytes()
        drugs = read_csv(d / "drugs.csv")
        rows = [{c: r[c] for c in PROPERTY_COLUMNS}
                for r in extract_properties(xml, s["property_table"], drugs)]
        h = properties_sha256(rows)
        if record and s["source"].get("properties_sha256") != h:
            s["source"]["properties_sha256"] = h
            dump_study(d, s)
        if s["source"].get("properties_sha256") != h:
            print(f"{d.name}: extracted properties DIFFER FROM source.properties_sha256")
            bad += 1
        eq_bad = equation_problems(xml, s)
        if eq_bad:
            print(f"{d.name}: equation_verbatim not found in the XML for {eq_bad}")
            bad += 1
        if s["licence"]["data_redistribution"] == "withheld":
            q = XML_CACHE / s["study_key"] / "properties.csv"
            write_csv(q, rows, PROPERTY_COLUMNS)
            print(f"{d.name}: {len(rows)} rows written to {q} (not committed: licence)")
            continue
        same = read_csv(d / "properties.csv") == rows
        print(f"{d.name}: properties.csv {'matches' if same else 'DIFFERS FROM'} the XML")
        bad += not same
    return 1 if bad else 0


def cmd_validate(keys: list[str] | None) -> int:
    cache = PubChemCache(offline=True)
    n = 0
    dirs = study_dirs(keys)
    if not dirs:
        print("no study folders found")
        return 1
    for d in dirs:
        errs = validate_study(d, cache)
        for e in errs:
            print(e)
        n += len(errs)
        print(f"{d.name}: {'OK' if not errs else f'{len(errs)} problem(s)'}")
    return 1 if n else 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("validate", "resolve", "fetch-xml", "extract"):
        sp = sub.add_parser(name)
        sp.add_argument("studies", nargs="*")
        if name == "resolve":
            sp.add_argument("--offline", action="store_true",
                            help="use only the committed PubChem cache")
        if name == "extract":
            sp.add_argument("--record", action="store_true",
                            help="write source.properties_sha256 from the extracted rows")
    a = ap.parse_args(argv)
    keys = a.studies or None
    try:
        study_dirs(keys)
    except KeyError as exc:
        print(f"error: {exc.args[0]}", file=sys.stderr)
        return 2
    if a.cmd == "validate":
        return cmd_validate(keys)
    if a.cmd == "resolve":
        return cmd_resolve(keys, a.offline)
    if a.cmd == "fetch-xml":
        return cmd_fetch_xml(keys)
    return cmd_extract(keys, a.record)


if __name__ == "__main__":
    sys.exit(main())
