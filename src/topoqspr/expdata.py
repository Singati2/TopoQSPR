"""Measured (experimental) property targets for the headline analysis.

Class-3 targets of ``docs/EXPERIMENT_PLAN.md`` section 4: open, licence-checked datasets
of *measured* small-molecule properties, matched to drugs by the full InChIKey of the
standardised parent (:func:`topoqspr.graph.standardize`).  Sources, licences, curation
rules and drug coverage are documented in ``docs/EXPERIMENTAL_DATA.md``.

Interface (what the modelling pipeline uses)
--------------------------------------------
``get_target(name)`` returns one row per standardised parent::

    inchikey | value | n_measurements | source | sd | replicates_disagree

sorted by ``inchikey``; ``df.attrs`` carries ``target``, ``unit``, ``source`` and the
SHA-256 of the raw file.  ``load_records(name)`` returns the record-level audit table
(one row per source record, with a ``status`` reason code and a ``detail`` text for every
record that was not used).  :data:`TARGETS` lists the target names.

Pipeline per record (first failure wins; nothing is imputed)
------------------------------------------------------------
1. Raw value parsed from the pinned file.  Non-numeric or missing -> ``missing_value``;
   qualified values (``>``, ``<``, ``~``) and values at or beyond a target's
   ``assay_range`` (AstraZeneca assay limits) -> ``censored_value``; records the source itself
   marks as bad (Bradley ``donotuse``, ChEMBL ``data_validity_comment``) ->
   ``flagged_by_source``; an OPERA record listing several alternative structures
   (``a|b``) -> ``ambiguous_structure``.
2. Unit normalisation to the target unit (log10 mol/L for solubility, degC for melting
   point, log units for logP / logD).  A conversion that needs a molar mass uses the
   recorded structure, i.e. the form that was weighed; failure -> ``unit_conversion_failed``.
3. Structure: :func:`topoqspr.graph.standardize` on the recorded SMILES; its rejection
   codes are used unchanged (``unparsable``, ``disallowed_element``, ``mixture`` ...).
   A standardised parent without an InChIKey -> ``no_inchikey``.
4. Solid-state targets (solubility, melting point) only: a record whose recorded
   structure has more than one component (salt, hydrate, solvate, co-crystal) ->
   ``salt_or_solvate_form``, because the measured value belongs to that solid form, not
   to the parent.  ``keep_salt_forms=True`` disables this rule (sensitivity analysis).
5. Aggregation by full InChIKey of the standardised parent: ``value`` is the **median**
   of the record values; ``n_measurements`` is the number of underlying measurements
   (a record that is itself an aggregate, e.g. an AqSolDB row with ``Ocurrences`` = 3,
   counts as its number of occurrences, an OPERA record as its number of merged PHYSPROP
   IDs); ``sd`` is the larger of the sample SD (ddof=1)
   across records and any SD the source reports for a record (NaN if neither exists).
   ``replicates_disagree`` is ``sd > disagreement_sd`` of the target (0.5 log units for
   logS / logP / logD - the AqSolDB reliability cut-off; 10 degC for melting point).
   ``get_target`` drops disagreeing parents by default (``drop_disagreeing=True``).

Downloads
---------
Every file source has a pinned URL and SHA-256.  Files are stored under
``<cache_dir>/expdata/<source>/`` (default cache: :func:`topoqspr.chembl.default_cache_dir`),
never in the repository.  A download whose hash does not match is deleted and
:class:`ChecksumError` is raised; a cached file whose hash does not match also raises
(delete it by hand to re-download).  ChEMBL-hosted sources go through
:class:`topoqspr.chembl.ChemblClient` (release-keyed cache) and the canonical JSON of the
fetched records is checked against a SHA-256 pinned per ChEMBL release.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import zipfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd
import requests
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors

from .chembl import CacheMiss, ChemblClient, default_cache_dir
from .graph import REJECTION_CODES, standardize

__all__ = [
    "SOURCES",
    "TARGETS",
    "STATUS_CODES",
    "TARGET_COLUMNS",
    "RECORD_COLUMNS",
    "ChecksumError",
    "Source",
    "TargetSpec",
    "aggregate",
    "approved_drug_pool",
    "convert_unit",
    "coverage",
    "fetch",
    "get_target",
    "load_records",
    "process_records",
    "processing_fingerprint",
    "sha256_bytes",
    "standardised_keys",
    "chembl_records_sha256",
    "parse_aqsoldb",
    "parse_biogen",
    "parse_bradley",
    "parse_chembl_activities",
    "parse_opera",
    "read_sdf_properties",
]

#: Columns of :func:`get_target` (the first four are the pipeline contract).
TARGET_COLUMNS = ("inchikey", "value", "n_measurements", "source", "sd", "replicates_disagree")

#: Columns of :func:`load_records` / :func:`process_records`.
RECORD_COLUMNS = (
    "source",
    "record_id",
    "smiles",
    "raw_value",
    "raw_unit",
    "value",
    "n_reported",
    "sd_reported",
    "inchikey",
    "canonical_smiles",
    "status",
    "detail",
)

#: Record status codes: ``ok`` plus every reason a record is not used.
STATUS_CODES: tuple[str, ...] = (
    "ok",
    "missing_structure",
    "ambiguous_structure",
    "missing_value",
    "censored_value",
    "flagged_by_source",
    "unit_conversion_failed",
    *REJECTION_CODES,
    "no_inchikey",
    "salt_or_solvate_form",
)


class ChecksumError(RuntimeError):
    """A downloaded or cached source file does not have its pinned SHA-256."""


# ----------------------------------------------------------------------------- registry
@dataclass(frozen=True)
class Source:
    """A pinned, licence-checked data source (see ``docs/EXPERIMENTAL_DATA.md``)."""

    key: str
    title: str
    citation: str
    doi: str
    version: str
    licence: str
    redistribution: bool
    url: str = ""  # file sources: pinned download URL
    filename: str = ""
    sha256: str = ""
    chembl_assay: str = ""  # ChEMBL sources: assay ChEMBL ID
    chembl_sha256: Mapping[str, str] | None = None  # release -> sha256 of canonical records


@dataclass(frozen=True)
class TargetSpec:
    """A measured endpoint taken from one source."""

    name: str
    source: str
    endpoint: str
    unit: str
    solid_state: bool
    disagreement_sd: float
    description: str
    member: str = ""  # zip member (OPERA) or column name (Biogen)
    #: Finite measuring range of the assay, in ``unit``.  A value at or beyond either limit
    #: is a limit value, not a measurement, and becomes ``censored_value`` even when the
    #: source stores it with ``=``.
    assay_range: tuple[float, float] | None = None


SOURCES: dict[str, Source] = {
    s.key: s
    for s in (
        Source(
            key="aqsoldb",
            title="AqSolDB: curated reference set of aqueous solubility",
            citation="Sorkun MC, Khetan A, Er S. Sci. Data 6 (2019) 143",
            doi="10.1038/s41597-019-0151-1",
            version="Harvard Dataverse doi:10.7910/DVN/OVHAW8 V1 (2019-07-18), "
            "file curated-solubility-dataset (original CSV, id 3407241)",
            licence="CC0 1.0",
            redistribution=True,
            url="https://dataverse.harvard.edu/api/access/datafile/3407241?format=original",
            filename="curated-solubility-dataset.csv",
            sha256="7c47e021d03a5720f53cfddb7aef90c8911ea276a36306c839798733d8630332",
        ),
        Source(
            key="bradley",
            title="Jean-Claude Bradley Open Melting Point Dataset",
            citation="Bradley J-C, Williams A, Lang A. figshare dataset (2014)",
            doi="10.6084/m9.figshare.1031637.v2",
            version="figshare article 1031637 v2 (2014-05-20), file BradleyMeltingPointDataset.xlsx",
            licence="CC0",
            redistribution=True,
            url="https://ndownloader.figshare.com/files/1503990",
            filename="BradleyMeltingPointDataset.xlsx",
            sha256="e8a66818609bbad6b9a67844cd98cf5d18098724b9f5fd8e0413aaf97775aa12",
        ),
        Source(
            key="opera",
            title="OPERA curated PHYSPROP training/test sets (QSAR-ready)",
            citation="Mansouri K et al. J. Cheminform. 10 (2018) 10; "
            "Mansouri K et al. SAR QSAR Environ. Res. 27 (2016) 911-937",
            doi="10.1186/s13321-018-0263-1",
            version="github.com/kmansouri/OPERA commit 22f4db60e385 (v2.9.5, 2026-08-20), OPERA_Data.zip",
            licence="MIT (repository LICENSE)",
            redistribution=True,
            url="https://github.com/kmansouri/OPERA/raw/22f4db60e38548f03921c44aecc84eb1b792e955/OPERA_Data.zip",
            filename="OPERA_Data.zip",
            sha256="9dc5d6387201df0f66e9472b034aeaa09ef7a47f705b9b5f50d7c9af8a1aadc3",
        ),
        Source(
            key="biogen",
            title="Biogen public ADME set (3521 compounds)",
            citation="Fang C et al. J. Chem. Inf. Model. 63 (2023) 3263-3274",
            doi="10.1021/acs.jcim.3c00160",
            version="github.com/molecularinformatics/Computational-ADME commit 685c9c828b23 (2023-03-31)",
            licence="MIT (repository LICENSE)",
            redistribution=True,
            url="https://raw.githubusercontent.com/molecularinformatics/Computational-ADME/"
            "685c9c828b23d8672572c21f82bb7e0780dfba2f/ADME_public_set_3521.csv",
            filename="ADME_public_set_3521.csv",
            sha256="2cfabc2667740c224487876c33b23124159ef43294e0f9e4d926cb6276c95a3b",
        ),
        Source(
            key="az_logd",
            title="AstraZeneca deposited dataset: octanol/water logD7.4 (shake flask)",
            citation="Wenlock M, Tomkinson N. ChEMBL deposited dataset CHEMBL3301361 (ChEMBL 20, 2015)",
            doi="10.6019/CHEMBL3301361",
            version="ChEMBL assay CHEMBL3301363 (records pinned per ChEMBL release)",
            licence="CC BY-SA 3.0 (ChEMBL)",
            redistribution=True,
            chembl_assay="CHEMBL3301363",
            chembl_sha256={"ChEMBL_37": "71e5e47c44b7578912be28d75120dff49e498c5fea0fed5c38e2c4661a0af33a"},
        ),
        Source(
            key="az_solubility",
            title="AstraZeneca deposited dataset: solubility in pH 7.4 buffer (solid starting material)",
            citation="Wenlock M, Tomkinson N. ChEMBL deposited dataset CHEMBL3301361 (ChEMBL 20, 2015)",
            doi="10.6019/CHEMBL3301361",
            version="ChEMBL assay CHEMBL3301364 (records pinned per ChEMBL release)",
            licence="CC BY-SA 3.0 (ChEMBL)",
            redistribution=True,
            chembl_assay="CHEMBL3301364",
            chembl_sha256={"ChEMBL_37": "6f63966724685f214fbd5ce49e26b1155d51df622ceb2b2841b0a87582fde256"},
        ),
    )
}

TARGETS: dict[str, TargetSpec] = {
    t.name: t
    for t in (
        TargetSpec("logs_aqsoldb", "aqsoldb", "aqueous solubility (intrinsic/water, ~25 degC)",
                   "log10(mol/L)", True, 0.5, "AqSolDB 'Solubility' (LogS)"),
        TargetSpec("logs_opera", "opera", "water solubility (PHYSPROP)",
                   "log10(mol/L)", True, 0.5, "OPERA WS_QR.sdf 'LogMolar'", member="WS_QR.sdf"),
        TargetSpec("logs74_az", "az_solubility", "solubility, pH 7.4 buffer",
                   "log10(mol/L)", True, 0.5, "ChEMBL standard_value (nM) -> log10(mol/L)",
                   assay_range=(-7.0, math.log10(1513600e-9))),  # 100 nM .. 1.5136 mM
        TargetSpec("logs68_biogen", "biogen", "kinetic solubility, pH 6.8",
                   "log10(mol/L)", True, 0.5, "Biogen log10(ug/mL) -> log10(mol/L) with the recorded structure's MW",
                   member="LOG SOLUBILITY PH 6.8 (ug/mL)"),
        TargetSpec("mp_bradley", "bradley", "melting point", "degC", True, 10.0,
                   "Bradley 'mpC'; rows marked donotuse are excluded"),
        TargetSpec("mp_opera", "opera", "melting point (PHYSPROP)", "degC", True, 10.0,
                   "OPERA MP_QR.sdf 'MP'", member="MP_QR.sdf"),
        TargetSpec("logp_opera", "opera", "octanol/water partition coefficient (PHYSPROP)",
                   "log10", False, 0.5, "OPERA LogP_QR.sdf 'LogP'", member="LogP_QR.sdf"),
        TargetSpec("logd74_az", "az_logd", "octanol/water distribution coefficient, pH 7.4",
                   "log10", False, 0.5, "ChEMBL standard_value (LogD7.4)", assay_range=(-1.5, 4.5)),
    )
}

#: Bump when record processing changes; part of the processed-records cache key (which
#: also carries :func:`processing_fingerprint`, so a code change is caught even unbumped).
PROCESSING_VERSION = "4"

#: Tolerance (target units) for "at the assay limit" (:attr:`TargetSpec.assay_range`).
_LIMIT_TOL = 1e-6


# ----------------------------------------------------------------------------- download
def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _expdata_dir(cache_dir: str | os.PathLike | None) -> Path:
    return Path(cache_dir) if cache_dir is not None else default_cache_dir() / "expdata"


def fetch(
    source: str | Source,
    cache_dir: str | os.PathLike | None = None,
    *,
    offline: bool = False,
    session: requests.Session | None = None,
    timeout: float = 300.0,
) -> Path:
    """Path of the verified local copy of a file source, downloading it if needed.

    The file lives at ``<cache_dir>/<source.key>/<source.filename>`` (default cache
    ``default_cache_dir()/expdata``).  A cached file with the wrong SHA-256 raises
    :class:`ChecksumError` (it is not silently replaced).  A download with the wrong
    SHA-256 is deleted, then :class:`ChecksumError` is raised.  ``offline=True`` raises
    :class:`topoqspr.chembl.CacheMiss` if the file is not cached.
    """
    src = SOURCES[source] if isinstance(source, str) else source
    if not src.url or not src.sha256:
        raise ValueError(f"source {src.key!r} is not a file source")
    path = _expdata_dir(cache_dir) / src.key / src.filename
    if path.exists():
        got = _sha256_file(path)
        if got != src.sha256:
            raise ChecksumError(
                f"cached {path} has sha256 {got}, expected {src.sha256}; "
                "delete the file to download it again"
            )
        return path
    if offline:
        raise CacheMiss(f"{src.key}: {path} not cached (offline)")
    sess = session if session is not None else requests.Session()
    resp = sess.get(src.url, timeout=timeout)
    resp.raise_for_status()
    data = resp.content
    got = sha256_bytes(data)
    if got != src.sha256:
        raise ChecksumError(f"download of {src.url} has sha256 {got}, expected {src.sha256}; not stored")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp{os.getpid()}")
    tmp.write_bytes(data)
    os.replace(tmp, path)
    return path


# ----------------------------------------------------------------------------- units
_LN10 = math.log(10.0)


def convert_unit(value: float, unit_from: str, unit_to: str, *, mw: float | None = None) -> float:
    """Convert one value between the units used here; NaN in -> NaN out.

    Concentrations: ``mol/L``, ``mM``, ``uM``, ``nM``, ``g/L``, ``mg/L``, ``ug/mL`` and
    their ``log10(...)`` forms -> ``log10(mol/L)`` (mass units need ``mw`` in g/mol).
    Temperatures: ``degC``, ``K``, ``degF`` -> ``degC``.  Identity for equal units.
    Raises ``ValueError`` for an unknown pair, a non-positive concentration, or a missing
    or non-positive ``mw``.
    """
    v = float(value)
    if math.isnan(v):
        return math.nan
    if unit_from == unit_to:
        return v
    if unit_to == "degC":
        if unit_from == "K":
            return v - 273.15
        if unit_from == "degF":
            return (v - 32.0) * 5.0 / 9.0
        raise ValueError(f"cannot convert {unit_from!r} to degC")
    if unit_to != "log10(mol/L)":
        raise ValueError(f"unsupported target unit {unit_to!r}")
    m = re.fullmatch(r"log10\((.+)\)", unit_from)
    base = m.group(1) if m else unit_from
    molar = {"mol/L": 1.0, "mM": 1e-3, "uM": 1e-6, "nM": 1e-9}
    mass = {"g/L": 1.0, "mg/L": 1e-3, "ug/mL": 1e-3}  # -> g/L
    if base in molar:
        factor = molar[base]
    elif base in mass:
        if mw is None or not mw > 0:
            raise ValueError(f"converting {unit_from!r} needs a positive molar mass")
        factor = mass[base] / float(mw)
    else:
        raise ValueError(f"cannot convert {unit_from!r} to log10(mol/L)")
    if m:
        return v + math.log10(factor)
    if v <= 0:
        raise ValueError(f"non-positive concentration {v} {unit_from} has no logarithm")
    return math.log10(v * factor)


# ----------------------------------------------------------------------------- parsers
# Every parser returns a raw-record frame with columns
#   record_id, smiles, raw_value, raw_unit, value, n_reported, sd_reported, status, detail
# where status is "" (not yet decided) or an early reason code.
_RAW = ("record_id", "smiles", "raw_value", "raw_unit", "value", "n_reported", "sd_reported", "status", "detail")


def _num(x: Any) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return math.nan
    return v if math.isfinite(v) else math.nan


def _frame(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=list(_RAW))


def _blank(x: Any) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x)) or str(x).strip() == ""


def parse_aqsoldb(df: pd.DataFrame) -> pd.DataFrame:
    """AqSolDB CSV: ``ID, SMILES, Solubility`` (LogS, mol/L), ``SD``, ``Ocurrences``."""
    rows = []
    for r in df.itertuples(index=False):
        v = _num(r.Solubility)
        rows.append(
            dict(
                record_id=str(r.ID), smiles="" if _blank(r.SMILES) else str(r.SMILES),
                raw_value=r.Solubility, raw_unit="log10(mol/L)", value=v,
                n_reported=int(r.Ocurrences) if not _blank(r.Ocurrences) else 1,
                sd_reported=_num(r.SD), status="" if not math.isnan(v) else "missing_value",
                detail="",
            )
        )
    return _frame(rows)


def parse_bradley(df: pd.DataFrame) -> pd.DataFrame:
    """Bradley xlsx sheet: ``key, smiles, mpC, donotuse, donotusebecause``."""
    rows = []
    for r in df.itertuples(index=False):
        v = _num(r.mpC)
        status, detail = "", ""
        if math.isnan(v):
            status = "missing_value"
        elif not _blank(r.donotuse):
            status, detail = "flagged_by_source", f"donotuse: {'' if _blank(r.donotusebecause) else r.donotusebecause}"
        rows.append(
            dict(
                record_id=str(r.key), smiles="" if _blank(r.smiles) else str(r.smiles),
                raw_value=r.mpC, raw_unit="degC", value=v, n_reported=1, sd_reported=math.nan,
                status=status, detail=detail,
            )
        )
    return _frame(rows)


def read_sdf_properties(text: str) -> pd.DataFrame:
    """Data fields of an SD file as strings (one row per record; molblocks ignored)."""
    recs: list[dict] = []
    cur: dict[str, list[str]] = {}
    key: str | None = None
    for line in text.splitlines():
        if line.startswith("$$$$"):
            recs.append({k: "\n".join(v) for k, v in cur.items()})
            cur, key = {}, None
            continue
        m = re.match(r">\s*<([^>]+)>", line)
        if m:
            key = m.group(1)
            cur[key] = []
            continue
        if key is not None:
            if line.strip() == "":
                key = None
            else:
                cur[key].append(line)
    return pd.DataFrame(recs).fillna("")


_OPERA_VALUE = {"WS_QR.sdf": ("LogMolar", "log10(mol/L)"), "MP_QR.sdf": ("MP", "degC"), "LogP_QR.sdf": ("LogP", "log10")}


def parse_opera(df: pd.DataFrame, member: str) -> pd.DataFrame:
    """OPERA QSAR-ready SDF fields.  The structure is ``Original_SMILES`` (the recorded
    substance, salts included) when present, else ``Canonical_QSARr``.  OPERA records whose
    ``Original_SMILES`` lists several structures separated by ``|`` (stereoisomers merged
    by the QSAR-ready workflow) -> ``ambiguous_structure``.

    OPERA values are pre-aggregated: a record whose ``MPID`` lists several PHYSPROP entries
    (``a|b|c``) merges that many entries, so ``n_reported`` is the number of distinct
    PHYSPROP IDs (1 when the field is blank).  OPERA publishes no spread for a merged
    value, so ``sd_reported`` stays NaN."""
    col, unit = _OPERA_VALUE[member]
    rows = []
    for i, r in enumerate(df.to_dict("records")):
        smi = r.get("Original_SMILES") or r.get("Canonical_QSARr") or ""
        v = _num(r.get(col))
        rid = r.get("ChemID") or r.get("MPID") or str(i)
        status, detail = ("" if not math.isnan(v) else "missing_value"), ""
        if not status and "|" in smi:
            # OPERA merged several stereoisomers / tautomers into one record: the value
            # cannot be assigned to one full InChIKey
            status = "ambiguous_structure"
            detail = f"{smi.count('|') + 1} alternative structures in Original_SMILES"
        mpids = {x.strip() for x in str(r.get("MPID") or "").split("|") if x.strip()}
        rows.append(
            dict(
                record_id=f"{rid}#{i}", smiles=smi, raw_value=r.get(col), raw_unit=unit, value=v,
                n_reported=max(1, len(mpids)), sd_reported=math.nan, status=status, detail=detail,
            )
        )
    return _frame(rows)


@lru_cache(maxsize=200_000)
def _input_mw(smiles: str) -> float:
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles) if smiles else None
    return Descriptors.MolWt(mol) if mol is not None else math.nan


def parse_biogen(df: pd.DataFrame, column: str) -> pd.DataFrame:
    """Biogen public ADME CSV; ``column`` holds log10(ug/mL), converted to log10(mol/L)
    with the molar mass of the recorded SMILES."""
    rows = []
    for r in df.to_dict("records"):
        smi = "" if _blank(r.get("SMILES")) else str(r["SMILES"])
        raw = _num(r.get(column))
        status, detail, v = "", "", math.nan
        if math.isnan(raw):
            status = "missing_value"
        else:
            mw = _input_mw(smi)
            try:
                v = convert_unit(raw, "log10(ug/mL)", "log10(mol/L)", mw=mw)
            except ValueError as exc:
                status, detail = "unit_conversion_failed", str(exc)
        rows.append(
            dict(
                record_id=str(r.get("Internal ID")), smiles=smi, raw_value=r.get(column),
                raw_unit="log10(ug/mL)", value=v, n_reported=1, sd_reported=math.nan,
                status=status, detail=detail,
            )
        )
    return _frame(rows)


#: Fields kept from ChEMBL ``activity`` records.
ACTIVITY_FIELDS = (
    "activity_id",
    "molecule_chembl_id",
    "canonical_smiles",
    "standard_relation",
    "standard_value",
    "standard_units",
    "standard_type",
    "data_validity_comment",
)


def parse_chembl_activities(records: list[dict], unit_to: str) -> pd.DataFrame:
    """ChEMBL activity rows -> raw records.  ``standard_relation`` other than ``=`` (or
    empty) -> ``censored_value``; any ``data_validity_comment`` -> ``flagged_by_source``."""
    rows = []
    for a in sorted(records, key=lambda x: int(x["activity_id"])):
        raw = a.get("standard_value")
        unit = a.get("standard_units") or ("log10" if unit_to == "log10" else "")
        v = _num(raw)
        status, detail = "", ""
        rel = (a.get("standard_relation") or "=").strip()
        if math.isnan(v):
            status = "missing_value"
        elif rel != "=":
            status, detail = "censored_value", f"standard_relation {rel!r}"
        elif a.get("data_validity_comment"):
            status, detail = "flagged_by_source", f"data_validity_comment: {a['data_validity_comment']}"
        else:
            try:
                v = convert_unit(v, unit, unit_to)
            except ValueError as exc:
                status, detail, v = "unit_conversion_failed", str(exc), math.nan
        rows.append(
            dict(
                record_id=str(a["activity_id"]), smiles=a.get("canonical_smiles") or "",
                raw_value=raw, raw_unit=a.get("standard_units") or "", value=v if not status else math.nan,
                n_reported=1, sd_reported=math.nan, status=status, detail=detail,
            )
        )
    return _frame(rows)


def chembl_records_sha256(records: list[dict]) -> str:
    """SHA-256 of the canonical JSON of activity records (projected, sorted by activity_id)."""
    proj = [{k: r.get(k) for k in ACTIVITY_FIELDS} for r in sorted(records, key=lambda x: int(x["activity_id"]))]
    return sha256_bytes(json.dumps(proj, sort_keys=True, separators=(",", ":")).encode())


# ----------------------------------------------------------------------------- processing
def _std_summary(smiles: str) -> tuple[str, str, str, str, str, int]:
    """Picklable summary of :func:`standardize`: ``(code, detail, inchikey, canonical_smiles,
    parent_note, n_input_components)``; ``code`` is ``''`` on success."""
    res = standardize(smiles)
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles)
    n_comp = len(Chem.GetMolFrags(mol)) if mol is not None else 0
    code, detail = (res.rejection.code, res.rejection.detail) if res.rejection is not None else ("", "")
    return code, detail, res.inchikey or "", res.canonical_smiles or "", res.parent_note or "", n_comp


#: Probe structures covering every branch of :func:`topoqspr.graph.standardize` (pipeline
#: salt stripping, each parent-fallback rule, every rejection code, isotopes, stereo).
_STD_PROBES: tuple[str, ...] = (
    "CC(=O)Oc1ccccc1C(=O)O",
    "C[C@H](N)C(=O)O",
    "[2H]C([2H])([2H])c1ccccc1",
    "CN(C)C(=N)N=C(N)N.Cl",
    "O=C([O-])c1ccccc1.[Na+]",
    "[H+].[O-]C(=O)c1ccccc1",
    "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1.O=C(O)c1cc2ccccc2cc1O",
    "CCCCN.CCCCN",
    "C[N+](C)(C)CCO.O=C([O-])c1ccccc1O",
    "CC(=O)Nc1ccc(O)cc1.CCCCCCCCO",
    "CCO.c1ccccc1",
    "[Li+].[Li+].[O-]C([O-])=O",
    "C[Hg]C",
    "OB(O)c1ccccc1",
    "C[Si](C)(C)c1ccccc1",
    "C[Se]c1ccccc1",
    "O=S(=O)(O)O",
    "CC",
    "c1cccc2ncnc12",
    "CCO junk",
    "C[n+]1ccccc1",
)


def processing_fingerprint() -> str:
    """SHA-256 identifying the structure processing in effect: the ``topoqspr`` version,
    the source of ``graph.py`` and ``expdata.py``, and the output of the current
    :func:`standardize` on :data:`_STD_PROBES`.  Part of the processed-records cache key,
    so editing (or monkeypatching) the standardiser invalidates cached record tables.

    The probe output (about 2 s of standardisation) is memoised per ``standardize``
    function object and source hash, so it is recomputed whenever either changes."""
    from . import __version__, graph

    parts: list[Any] = [__version__]
    for f in (graph.__file__, __file__):
        parts.append(_sha256_file(Path(f)) if f and Path(f).exists() else "?")
    memo_key = (standardize, tuple(parts))
    probe = _PROBE_MEMO.get(memo_key)
    if probe is None:
        probe = [list(_std_summary(s)) for s in _STD_PROBES]
        _PROBE_MEMO.clear()
        _PROBE_MEMO[memo_key] = probe
    return sha256_bytes(json.dumps([*parts, probe], separators=(",", ":")).encode())


_PROBE_MEMO: dict[tuple, list] = {}  # (standardize, source hashes) -> probe output


_STD_MEMO: dict[str, tuple] = {}  # per-process memo of _std_summary
_STD_MEMO_FP: list[str] = [""]  # processing_fingerprint the memo was filled under
_POOL_MIN = 500  # below this many distinct SMILES a pool costs more than it saves


def _std_many(smiles: Iterable[str], n_jobs: int | None = None) -> dict[str, tuple]:
    """:func:`_std_summary` for each distinct SMILES, in a process pool when worthwhile.

    The result does not depend on ``n_jobs`` (a pure function mapped over sorted inputs).
    ``None`` -> ``min(8, os.cpu_count())``.
    """
    want = sorted({s for s in smiles if s})
    fp = processing_fingerprint()
    if fp != _STD_MEMO_FP[0]:  # standardiser changed: memo entries are stale
        _STD_MEMO.clear()
        _STD_MEMO_FP[0] = fp
    uniq = [s for s in want if s not in _STD_MEMO]
    n_jobs = min(8, os.cpu_count() or 1) if n_jobs is None else int(n_jobs)
    if n_jobs <= 1 or len(uniq) < _POOL_MIN:
        _STD_MEMO.update((s, _std_summary(s)) for s in uniq)
    else:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            _STD_MEMO.update(zip(uniq, ex.map(_std_summary, uniq, chunksize=64)))
    return {s: _STD_MEMO[s] for s in want}


def process_records(
    raw: pd.DataFrame, spec: TargetSpec, *, keep_salt_forms: bool = False, n_jobs: int | None = None
) -> pd.DataFrame:
    """Standardise raw records (module docstring, steps 3-4) -> record table (:data:`RECORD_COLUMNS`).

    ``n_jobs`` worker processes standardise the distinct SMILES (default ``min(8, cpus)``);
    the output does not depend on it.
    """
    recs = raw.to_dict("records")
    todo = [(r["smiles"] or "").strip() for r in recs if not (r["status"] or "")]
    std = _std_many(todo, n_jobs)
    out = []
    for r in recs:
        status, detail = r["status"] or "", r["detail"] or ""
        value = r["value"]
        key, can = "", ""
        smi = (r["smiles"] or "").strip()
        if not status and spec.assay_range is not None:
            lo, hi = spec.assay_range
            if value <= lo + _LIMIT_TOL or value >= hi - _LIMIT_TOL:
                status = "censored_value"
                detail = f"{value:.6g} at or beyond the assay range [{lo:.6g}, {hi:.6g}] {spec.unit}"
        if not status and not smi:
            status = "missing_structure"
        if not status:
            code, det, key, can, note, n_comp = std[smi]
            if code:
                status, detail = code, det
            elif not key:
                status = "no_inchikey"
            elif spec.solid_state and not keep_salt_forms and n_comp > 1:
                status = "salt_or_solvate_form"
                detail = f"{n_comp} components in the recorded structure" + (f"; {note}" if note else "")
        if status:
            value = math.nan
        out.append(
            dict(
                source=spec.source, record_id=r["record_id"], smiles=smi, raw_value=r["raw_value"],
                raw_unit=r["raw_unit"], value=value, n_reported=int(r["n_reported"]),
                sd_reported=r["sd_reported"], inchikey=key, canonical_smiles=can,
                status=status or "ok", detail=detail,
            )
        )
    df = pd.DataFrame(out, columns=list(RECORD_COLUMNS))
    df["raw_value"] = df["raw_value"].map(lambda x: "" if _blank(x) else str(x))
    df["value"] = df["value"].astype(float)
    df["sd_reported"] = df["sd_reported"].astype(float)
    df["n_reported"] = df["n_reported"].astype(int)
    bad = set(df["status"]) - set(STATUS_CODES)
    if bad:  # pragma: no cover - guard against a new code slipping in
        raise RuntimeError(f"unknown status codes {sorted(bad)}")
    return df


def aggregate(records: pd.DataFrame, spec: TargetSpec) -> pd.DataFrame:
    """Records with ``status == 'ok'`` -> one row per InChIKey (module docstring, step 5)."""
    ok = records[records["status"] == "ok"]
    rows = []
    for key, g in ok.groupby("inchikey", sort=True):
        vals = g["value"].to_numpy(dtype=float)
        between = float(np.std(vals, ddof=1)) if len(vals) > 1 else math.nan
        within = g["sd_reported"].to_numpy(dtype=float)
        within = float(np.nanmax(within)) if np.isfinite(within).any() else math.nan
        sd = float(np.nanmax([between, within])) if not (math.isnan(between) and math.isnan(within)) else math.nan
        rows.append(
            dict(
                inchikey=key,
                value=float(np.median(vals)),
                n_measurements=int(g["n_reported"].sum()),
                source=spec.source,
                sd=sd,
                replicates_disagree=bool(not math.isnan(sd) and sd > spec.disagreement_sd),
            )
        )
    df = pd.DataFrame(rows, columns=list(TARGET_COLUMNS))
    df["value"] = df["value"].astype(float)
    df["sd"] = df["sd"].astype(float)
    df["n_measurements"] = df["n_measurements"].astype(int)
    df["replicates_disagree"] = df["replicates_disagree"].astype(bool)
    return df.reset_index(drop=True)


# ----------------------------------------------------------------------------- loading
def _raw_records(
    spec: TargetSpec,
    cache_dir: str | os.PathLike | None,
    offline: bool,
    client: ChemblClient | None,
    session: requests.Session | None,
) -> tuple[pd.DataFrame, str]:
    src = SOURCES[spec.source]
    if src.chembl_assay:
        own = client is None
        cl = client if client is not None else ChemblClient(offline=offline)
        recs = cl.fetch_all(
            "activity",
            {"assay_chembl_id": src.chembl_assay, "order_by": "activity_id", "only": ",".join(ACTIVITY_FIELDS)},
        )
        if own and not offline:
            cl.verify_release()
        digest = chembl_records_sha256(recs)
        pinned = (src.chembl_sha256 or {}).get(cl.release)
        if pinned is None:
            raise ChecksumError(f"{src.key}: no pinned sha256 for {cl.release}; records hash {digest}")
        if digest != pinned:
            raise ChecksumError(f"{src.key}: {cl.release} records hash {digest}, expected {pinned}")
        return parse_chembl_activities(recs, spec.unit), f"{cl.release}:{digest}"
    path = fetch(src, cache_dir, offline=offline, session=session)
    if src.key == "aqsoldb":
        raw = parse_aqsoldb(pd.read_csv(path, dtype={"ID": str}, keep_default_na=True))
    elif src.key == "bradley":
        raw = parse_bradley(pd.read_excel(path, sheet_name="Sheet1", engine="openpyxl"))
    elif src.key == "opera":
        with zipfile.ZipFile(path) as zf:
            text = zf.read(spec.member).decode("latin-1")
        raw = parse_opera(read_sdf_properties(text), spec.member)
    elif src.key == "biogen":
        raw = parse_biogen(pd.read_csv(path), spec.member)
    else:  # pragma: no cover
        raise KeyError(src.key)
    return raw, src.sha256


def load_records(
    name: str,
    cache_dir: str | os.PathLike | None = None,
    *,
    offline: bool = False,
    keep_salt_forms: bool = False,
    client: ChemblClient | None = None,
    session: requests.Session | None = None,
    n_jobs: int | None = None,
) -> pd.DataFrame:
    """Record-level audit table of target ``name`` (:data:`RECORD_COLUMNS`), in source order.

    Processed tables are cached as ``<cache_dir>/processed/<name>-<hash>.csv.gz`` keyed by the
    raw data hash, ``keep_salt_forms``, :data:`PROCESSING_VERSION`, the RDKit / ChEMBL
    structure-pipeline versions and :func:`processing_fingerprint` (topoqspr version, source
    of ``graph.py`` / ``expdata.py``, standardiser output on fixed probes).
    ``df.attrs['processing_fingerprint']`` records the latter.
    """
    spec = TARGETS[name]
    raw, raw_hash = _raw_records(spec, cache_dir, offline, client, session)
    import chembl_structure_pipeline

    fp = processing_fingerprint()
    key = json.dumps(
        [name, raw_hash, keep_salt_forms, PROCESSING_VERSION, rdBase.rdkitVersion,
         getattr(chembl_structure_pipeline, "__version__", "?"), fp]
    )
    proc = _expdata_dir(cache_dir) / "processed" / f"{name}-{sha256_bytes(key.encode())[:16]}.csv.gz"
    if proc.exists():
        df = pd.read_csv(proc, dtype=str, keep_default_na=False)
        for c in ("value", "sd_reported"):
            df[c] = df[c].map(lambda x: float(x) if x != "" else math.nan).astype(float)
        df["n_reported"] = df["n_reported"].astype(int)
    else:
        df = process_records(raw, spec, keep_salt_forms=keep_salt_forms, n_jobs=n_jobs)
        proc.parent.mkdir(parents=True, exist_ok=True)
        tmp = proc.with_suffix(f".tmp{os.getpid()}")
        df.to_csv(tmp, index=False, compression={"method": "gzip", "mtime": 0})  # repr floats round-trip
        os.replace(tmp, proc)
    df.attrs.update(target=name, unit=spec.unit, source=spec.source, raw_sha256=raw_hash,
                    processing_fingerprint=fp)
    return df


def get_target(
    name: str,
    cache_dir: str | os.PathLike | None = None,
    *,
    offline: bool = False,
    drop_disagreeing: bool = True,
    keep_salt_forms: bool = False,
    client: ChemblClient | None = None,
    session: requests.Session | None = None,
    n_jobs: int | None = None,
) -> pd.DataFrame:
    """Measured target ``name`` as ``DataFrame[inchikey, value, n_measurements, source, sd,
    replicates_disagree]``, one row per standardised parent InChIKey, sorted by InChIKey.

    ``value`` is in ``TARGETS[name].unit`` (also in ``df.attrs['unit']``).  With
    ``drop_disagreeing`` (default) parents whose replicate SD exceeds the target's
    ``disagreement_sd`` are removed; ``df.attrs['n_dropped_disagreeing']`` records how many.
    """
    spec = TARGETS[name]
    recs = load_records(
        name, cache_dir, offline=offline, keep_salt_forms=keep_salt_forms, client=client, session=session,
        n_jobs=n_jobs,
    )
    agg = aggregate(recs, spec)
    n_dis = int(agg["replicates_disagree"].sum())
    if drop_disagreeing:
        agg = agg[~agg["replicates_disagree"]].reset_index(drop=True)
    agg.attrs.update(recs.attrs, n_dropped_disagreeing=n_dis if drop_disagreeing else 0,
                     n_disagreeing=n_dis, keep_salt_forms=keep_salt_forms)
    return agg


# ----------------------------------------------------------------------------- drugs
def approved_drug_pool(client: ChemblClient):
    """Pooled cohort of approved small-molecule parents for *any* indication.

    One scan of ``drug_indication`` with ``max_phase_for_ind=4`` (every MeSH heading),
    then exactly the curation of :func:`topoqspr.cohort.build_cohort` (parent collapse,
    small molecules, structure checks, InChIKey de-duplication).  This is the union of the
    controlled phase-4 cohorts of all headings, without an evidence filter.
    Returns a :class:`topoqspr.cohort.Cohort` with ``kind='pooled'``.
    """
    from .cohort import Cohort, _curate_families, _phase, _resolve_parents, _sort_manifest

    rows = client.fetch_all(
        "drug_indication",
        {
            "max_phase_for_ind": "4",
            "order_by": "drugind_id",
            "only": "drugind_id,molecule_chembl_id,parent_molecule_chembl_id,max_phase_for_ind,mesh_id",
        },
    )
    members: dict[str, dict] = {}
    fallback: dict[str, str] = {}
    for r in rows:
        cid = r["molecule_chembl_id"]
        m = members.setdefault(cid, {"phase": math.nan, "qualifying": True, "ref_types": set(), "efo_ids": set()})
        p = _phase(r.get("max_phase_for_ind"))
        if not math.isnan(p) and (math.isnan(m["phase"]) or p > m["phase"]):
            m["phase"] = p
        if r.get("parent_molecule_chembl_id"):
            fallback[cid] = r["parent_molecule_chembl_id"]
    records, parent_of = _resolve_parents(client, members, fallback)
    families: dict[str, dict] = {}
    for cid, m in members.items():
        families.setdefault(parent_of[cid], {})[cid] = m
    out = _curate_families(
        families, records, label="pooled approved (any MeSH heading)", mesh_col="*",
        release=client.release, min_phase=4.0, use_phase=True,
    )
    client.verify_release()
    return Cohort(
        kind="pooled",
        disease_label="pooled approved (any MeSH heading)",
        release=client.release,
        release_date=client.status().get("chembl_release_date", ""),
        parameters={
            "drug_indication_filter": {"max_phase_for_ind": 4},
            "min_phase": 4.0,
            "molecule_type": "Small molecule",
            "dedup": "standard_inchi_key, lowest parent ChEMBL ID kept",
        },
        manifest=_sort_manifest(out),
        extra={"indication_rows": len(rows), "n_mesh_ids": len({r.get("mesh_id") for r in rows})},
    )


def standardised_keys(included: pd.DataFrame, n_jobs: int | None = None) -> pd.DataFrame:
    """Standardise the ``canonical_smiles`` of included cohort parents with
    :func:`topoqspr.graph.standardize`: ``[parent_chembl_id, inchikey, status, detail]``
    (``inchikey`` empty and ``status`` the rejection code when standardisation fails)."""
    smiles = [r.canonical_smiles or "" for r in included.itertuples(index=False)]
    std = _std_many(smiles, n_jobs)
    rows = []
    for r, smi in zip(included.itertuples(index=False), smiles):
        code, det, key, _, _, _ = std[smi] if smi else ("unparsable", "empty SMILES", "", "", "", 0)
        st = code or ("ok" if key else "no_inchikey")
        rows.append(dict(parent_chembl_id=r.parent_chembl_id, inchikey=key if st == "ok" else "", status=st, detail=det))
    return pd.DataFrame(rows, columns=["parent_chembl_id", "inchikey", "status", "detail"])


def coverage(targets: Mapping[str, pd.DataFrame], inchikeys: Iterable[str]) -> dict[str, dict[str, int]]:
    """For each target, how many of ``inchikeys`` (full keys) have a value.

    Returns ``{target: {"n": drugs (distinct full keys) with a full-key match,
    "n_first_block": drugs whose first (connectivity) block matches, informational}}``.
    Both count drugs, i.e. distinct full keys; empty keys are ignored.
    """
    keys = sorted({k for k in inchikeys if k})
    out = {}
    for name, df in targets.items():
        have = set(df["inchikey"])
        have_b = {k.split("-")[0] for k in have}
        out[name] = {"n": sum(k in have for k in keys), "n_first_block": sum(k.split("-")[0] in have_b for k in keys)}
    return out
