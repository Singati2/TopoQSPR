"""Reproducible disease -> drug cohorts with a versioned manifest.

Three cohort definitions share one manifest format (see ``docs/COHORT_DEFINITION.md``):

* :func:`build_cohort` - the v1.0 *controlled* cohort.  The disease is an explicit list of
  MeSH descriptor IDs confirmed by the user (``drug_indication.mesh_id``); there is no
  substring retrieval and no keyword fallback.  Molecules are collapsed onto their ChEMBL
  parent (``molecule_hierarchy.parent_chembl_id``); the phase of a parent is

      phase(P) = max { max_phase_for_ind(r) : r in rows(selected mesh_ids),
                       molecule(r) in family(P), r passes the evidence filter,
                       max_phase_for_ind(r) >= 0 }

  (ChEMBL encodes "unknown" as -1 and early phase 1 as 0.5; Zdrazil et al. 2024,
  doi:10.1093/nar/gkad1004).  Parents are kept if ``phase >= min_phase``,
  ``molecule_type == 'Small molecule'`` and the parent structure passes the RDKit checks
  of :func:`structure_reason`; survivors are de-duplicated by standard InChIKey (lowest
  ChEMBL ID kept).
* :func:`build_atc_cohort` - alternative definition by WHO ATC code prefix (e.g. ``A10``,
  drugs used in diabetes), same parent / structure / de-duplication rules, no phase rule.
* :func:`naive_cohort` - the v0.x retrieval frozen for Experiment A
  (``mesh_heading__icontains``, all phases, no salt handling), but keeping every
  candidate and sorting, so it is deterministic.

Every candidate molecule appears in the manifest with ``decision`` (include/exclude) and a
``reason`` code.  Same inputs + same ChEMBL release => byte-identical CSV.
"""

from __future__ import annotations

import csv
import difflib
import hashlib
import io
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd
from rdkit import Chem, RDLogger

from . import __version__
from .chembl import ChemblClient, MeshClient, chembl_sort_key

__all__ = [
    "MANIFEST_COLUMNS",
    "REASONS",
    "Cohort",
    "build_atc_cohort",
    "build_cohort",
    "jaccard",
    "naive_cohort",
    "structure_reason",
    "suggest_headings",
]

MANIFEST_COLUMNS = (
    "disease_label",
    "mesh_ids",
    "chembl_release",
    "molecule_chembl_id",
    "parent_chembl_id",
    "pref_name",
    "canonical_smiles",
    "standard_inchi_key",
    "molecule_type",
    "phase_for_indication",
    "indication_ref_types",
    "decision",
    "reason",
)

#: Reason codes, in the order the checks are applied to a parent.  The first failing check
#: is reported.  ``duplicate_inchikey`` carries the kept parent as ``duplicate_inchikey:<ID>``.
REASONS = (
    "included",
    "collapsed_to_parent",
    "no_qualifying_evidence",
    "phase_unknown",
    "phase_below_min",
    "molecule_record_missing",
    "not_small_molecule",
    "no_smiles",
    "unparseable_smiles",
    "no_carbon",
    "metal_containing",
    "unsupported_element",
    "mixture",
    "duplicate_inchikey",
)

#: Elements allowed in a drug-like organic small molecule (H, B, C, N, O, F, Si, P, S, Cl,
#: Se, Br, I).  Boron (bortezomib, dutogliptin) and selenium (ebselen) occur in drugs.
ORGANIC_ELEMENTS = frozenset({1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53})
#: Non-metals and metalloids that are *not* in ORGANIC_ELEMENTS: He, Ne, Ar, Ge, As, Kr,
#: Sb, Te, Xe, At, Rn (and the dummy atom 0).  Anything else with Z > 0 is a metal.
_NONMETAL_OTHER = frozenset({0, 2, 10, 18, 32, 33, 36, 51, 52, 54, 85, 86})


def _phase(value: Any) -> float:
    """ChEMBL phase string/number -> float; ``None``/``-1`` (unknown) -> NaN."""
    if value is None or value == "":
        return math.nan
    x = float(value)
    return math.nan if x < 0 else x


def _fmt_phase(x: float) -> str:
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.1f}"


def structure_reason(smiles: str | None) -> str | None:
    """Return the structural exclusion reason for ``smiles`` or ``None`` if acceptable.

    Checks (first failure wins): ``no_smiles``; ``unparseable_smiles`` (RDKit sanitisation
    fails); ``no_carbon`` (inorganic); ``metal_containing`` (any atom outside
    :data:`ORGANIC_ELEMENTS` that is a metal); ``unsupported_element`` (other non-organic
    elements, e.g. As, Te, Xe, dummy atoms); ``mixture`` (more than one connected fragment
    in the *parent* structure, e.g. co-crystals or multi-component products).
    """
    if not smiles:
        return "no_smiles"
    RDLogger.DisableLog("rdApp.*")
    try:
        mol = Chem.MolFromSmiles(smiles)
    finally:
        RDLogger.EnableLog("rdApp.*")
    if mol is None:
        return "unparseable_smiles"
    zs = {a.GetAtomicNum() for a in mol.GetAtoms()}
    if 6 not in zs:
        return "no_carbon"
    foreign = zs - ORGANIC_ELEMENTS
    if foreign - _NONMETAL_OTHER:
        return "metal_containing"
    if foreign:
        return "unsupported_element"
    if len(Chem.GetMolFrags(mol)) > 1:
        return "mixture"
    return None


def _inchikey(rec: Mapping | None) -> str:
    """ChEMBL standard InChIKey; RDKit InChIKey of the canonical SMILES if ChEMBL has none."""
    structs = (rec or {}).get("molecule_structures") or {}
    key = structs.get("standard_inchi_key")
    if key:
        return key
    smi = structs.get("canonical_smiles")
    if smi:
        RDLogger.DisableLog("rdApp.*")
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                return Chem.MolToInchiKey(mol) or ""
        finally:
            RDLogger.EnableLog("rdApp.*")
    return ""


def _smiles(rec: Mapping | None) -> str:
    return ((rec or {}).get("molecule_structures") or {}).get("canonical_smiles") or ""


def _parent_of(cid: str, rec: Mapping | None, fallback: str | None) -> str:
    hier = (rec or {}).get("molecule_hierarchy") or {}
    return hier.get("parent_chembl_id") or fallback or cid


# ---------------------------------------------------------------------------------- result
@dataclass
class Cohort:
    """A cohort manifest plus the metadata written to its JSON sidecar."""

    kind: str
    disease_label: str
    release: str
    release_date: str
    parameters: dict
    manifest: pd.DataFrame
    extra: dict = field(default_factory=dict)

    @property
    def included(self) -> pd.DataFrame:
        """Included rows (one per curated parent), in manifest order."""
        return self.manifest[self.manifest["decision"] == "include"].reset_index(drop=True)

    def reason_counts(self) -> dict[str, int]:
        """Rows per reason code (detail after ``:`` stripped), sorted by code."""
        codes = self.manifest["reason"].str.split(":").str[0]
        return {k: int(v) for k, v in sorted(codes.value_counts().items())}

    def to_csv_bytes(self) -> bytes:
        """The manifest as UTF-8 CSV with ``\\n`` line endings (byte-stable)."""
        buf = io.StringIO()
        w = csv.writer(buf, lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
        w.writerow(MANIFEST_COLUMNS)
        for rec in self.manifest.itertuples(index=False):
            row = rec._asdict()
            row["phase_for_indication"] = _fmt_phase(row["phase_for_indication"])
            w.writerow(["" if row[c] is None else str(row[c]) for c in MANIFEST_COLUMNS])
        return buf.getvalue().encode("utf-8")

    def sidecar(self, timestamp_utc: str) -> dict:
        """JSON sidecar content.  ``timestamp_utc`` is supplied by the caller (ISO 8601)."""
        csv_bytes = self.to_csv_bytes()
        return {
            "kind": self.kind,
            "disease_label": self.disease_label,
            "chembl_release": self.release,
            "chembl_release_date": self.release_date,
            "parameters": self.parameters,
            "counts": {
                "manifest_rows": int(len(self.manifest)),
                "parents": int((self.manifest["molecule_chembl_id"] == self.manifest["parent_chembl_id"]).sum())
                if self.kind != "naive"
                else None,
                "included": int((self.manifest["decision"] == "include").sum()),
                "by_reason": self.reason_counts(),
            },
            **self.extra,
            "topoqspr_version": __version__,
            "created_utc": timestamp_utc,
            "manifest_sha256": hashlib.sha256(csv_bytes).hexdigest(),
        }

    def write(self, csv_path: str | Path, timestamp_utc: str, json_path: str | Path | None = None) -> tuple[Path, Path]:
        """Write ``<name>.csv`` and the JSON sidecar (default ``<name>.json``)."""
        csv_path = Path(csv_path)
        json_path = Path(json_path) if json_path is not None else csv_path.with_suffix(".json")
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        csv_path.write_bytes(self.to_csv_bytes())
        json_path.write_text(json.dumps(self.sidecar(timestamp_utc), sort_keys=True, indent=1) + "\n")
        return csv_path, json_path


def _sort_manifest(rows: list[dict]) -> pd.DataFrame:
    rows = sorted(
        rows,
        key=lambda r: (
            chembl_sort_key(r["parent_chembl_id"]),
            r["molecule_chembl_id"] != r["parent_chembl_id"],
            chembl_sort_key(r["molecule_chembl_id"]),
        ),
    )
    df = pd.DataFrame(rows, columns=list(MANIFEST_COLUMNS))
    df["phase_for_indication"] = df["phase_for_indication"].astype(float)
    return df


# ------------------------------------------------------------------------ shared curation
def _curate_families(
    families: Mapping[str, Mapping[str, dict]],
    records: Mapping[str, dict | None],
    *,
    label: str,
    mesh_col: str,
    release: str,
    min_phase: float | None,
    use_phase: bool,
) -> list[dict]:
    """Manifest rows for parent families.

    ``families[parent][member] = {"phase": float, "qualifying": bool, "ref_types": set}``
    where ``phase`` is the member's max phase over qualifying rows (NaN if none/unknown).
    """
    out: list[dict] = []
    kept: list[dict] = []
    for parent in sorted(families, key=chembl_sort_key):
        members = families[parent]
        prec = records.get(parent)
        fam_phase = [m["phase"] for m in members.values() if not math.isnan(m["phase"])]
        phase = max(fam_phase) if fam_phase else math.nan
        any_qual = any(m["qualifying"] for m in members.values())
        ref_types = sorted(set().union(*(m["ref_types"] for m in members.values())))
        base = {
            "disease_label": label,
            "mesh_ids": mesh_col,
            "chembl_release": release,
            "parent_chembl_id": parent,
        }
        reason = None
        if use_phase and not any_qual:
            reason = "no_qualifying_evidence"
        elif use_phase and min_phase is not None and math.isnan(phase):
            reason = "phase_unknown"
        elif use_phase and min_phase is not None and phase < min_phase:
            reason = "phase_below_min"
        elif prec is None:
            reason = "molecule_record_missing"
        elif prec.get("molecule_type") != "Small molecule":
            reason = "not_small_molecule"
        else:
            reason = structure_reason(_smiles(prec))
        row = {
            **base,
            "molecule_chembl_id": parent,
            "pref_name": (prec or {}).get("pref_name") or "",
            "canonical_smiles": _smiles(prec),
            "standard_inchi_key": _inchikey(prec),
            "molecule_type": (prec or {}).get("molecule_type") or "",
            "phase_for_indication": phase,
            "indication_ref_types": ";".join(ref_types),
            "decision": "exclude" if reason else "include",
            "reason": reason or "included",
        }
        out.append(row)
        if reason is None:
            kept.append(row)
        for member in sorted(members, key=chembl_sort_key):
            if member == parent:
                continue
            mrec = records.get(member)
            m = members[member]
            out.append(
                {
                    **base,
                    "molecule_chembl_id": member,
                    "pref_name": (mrec or {}).get("pref_name") or "",
                    "canonical_smiles": _smiles(mrec),
                    "standard_inchi_key": _inchikey(mrec),
                    "molecule_type": (mrec or {}).get("molecule_type") or "",
                    "phase_for_indication": m["phase"],
                    "indication_ref_types": ";".join(sorted(m["ref_types"])),
                    "decision": "exclude",
                    "reason": "collapsed_to_parent",
                }
            )
    # de-duplicate included parents by standard InChIKey; lowest ChEMBL ID wins (kept is sorted)
    seen: dict[str, str] = {}
    for row in kept:
        key = row["standard_inchi_key"]
        if not key:
            continue
        if key in seen:
            row["decision"] = "exclude"
            row["reason"] = f"duplicate_inchikey:{seen[key]}"
        else:
            seen[key] = row["molecule_chembl_id"]
    return out


def _resolve_parents(client: ChemblClient, ids: Iterable[str], fallback: Mapping[str, str]) -> tuple[dict, dict]:
    """Molecule records for ``ids`` and their parents, plus ``parent_of``.

    The parent is ``molecule_hierarchy.parent_chembl_id``; if the record has no hierarchy,
    ``fallback`` (``drug_indication.parent_molecule_chembl_id``) and then the molecule
    itself.  The mapping is followed to a fixed point (at most 5 hops), fetching each hop's
    records in one batch.
    """
    ids = sorted(set(ids), key=chembl_sort_key)
    records: dict[str, dict | None] = dict(client.molecules(ids))
    current = {cid: cid for cid in ids}
    done: set[str] = set()
    for _ in range(5):
        step: dict[str, str] = {}
        for cid in ids:
            if cid in done:
                continue
            cur = current[cid]
            nxt = _parent_of(cur, records.get(cur), fallback.get(cur))
            if nxt == cur:
                done.add(cid)
            else:
                step[cid] = nxt
        if not step:
            break
        need = sorted({p for p in step.values() if p not in records}, key=chembl_sort_key)
        if need:
            records.update(client.molecules(need))
        current.update(step)
    else:
        raise RuntimeError("molecule_hierarchy did not reach a fixed point in 5 hops")
    return records, current


# ------------------------------------------------------------------------ disease lookup
_TOKEN = re.compile(r"[a-z0-9]+")


def _tokens(text: str) -> list[str]:
    return _TOKEN.findall(text.lower())


def _match_score(query: str, name: str) -> float:
    """``0.6 * |Q n N| / |Q| + 0.4 * SequenceMatcher(sorted Q, sorted N).ratio()``.

    ``Q``/``N`` are the lower-cased alphanumeric token sets of query and name; sorting the
    tokens makes the fuzzy part order-insensitive ("type 2 diabetes" vs "Diabetes Mellitus,
    Type 2").
    """
    q, n = _tokens(query), _tokens(name)
    if not q or not n:
        return 0.0
    tok = len(set(q) & set(n)) / len(set(q))
    fuzzy = difflib.SequenceMatcher(None, " ".join(sorted(set(q))), " ".join(sorted(set(n)))).ratio()
    return 0.6 * tok + 0.4 * fuzzy


def suggest_headings(
    free_text: str,
    client: ChemblClient,
    *,
    limit: int = 10,
    mesh_client: MeshClient | None = None,
    min_score: float = 0.3,
) -> pd.DataFrame:
    """Rank the release's ``drug_indication`` MeSH headings against ``free_text``.

    Suggestions only - the cohort is built from the descriptor IDs the user confirms.
    Score: :func:`_match_score` against the heading; with ``mesh_client``, also against
    every MeSH term (entry term / synonym) of descriptors whose terms contain the text
    (NLM lookup API), taking the maximum.  Sorted by score desc, rows desc, mesh_id.
    Columns: mesh_id, mesh_heading, n_rows, score, matched_via.
    """
    synonyms: dict[str, list[str]] = {}
    if mesh_client is not None:
        synonyms = mesh_client.descriptors_for_text(free_text)
    out = []
    for h in client.mesh_headings():
        if not h["mesh_id"]:
            continue
        best, via = _match_score(free_text, h["mesh_heading"]), h["mesh_heading"]
        for term in synonyms.get(h["mesh_id"], []):
            s = _match_score(free_text, term)
            if s > best:
                best, via = s, f"MeSH term: {term}"
        if best >= min_score:
            out.append({**h, "score": round(best, 6), "matched_via": via})
    out.sort(key=lambda r: (-r["score"], -r["n_rows"], r["mesh_id"]))
    return pd.DataFrame(out[:limit], columns=["mesh_id", "mesh_heading", "n_rows", "score", "matched_via"])


_MESH_ID = re.compile(r"^D\d{6,9}$")


def _resolve_mesh_ids(
    mesh_ids: Iterable[str],
    client: ChemblClient,
    include_descendants: bool,
    mesh_client: MeshClient | None,
) -> tuple[list[str], dict[str, str], dict]:
    requested = sorted({str(m).strip() for m in mesh_ids})
    if not requested:
        raise ValueError("mesh_ids is empty: a cohort needs at least one confirmed MeSH descriptor ID")
    bad = [m for m in requested if not _MESH_ID.match(m)]
    if bad:
        raise ValueError(
            f"not MeSH descriptor IDs: {bad}. Use suggest_headings() to find IDs; "
            "free text is never used to retrieve drugs."
        )
    index = {h["mesh_id"]: h["mesh_heading"] for h in client.mesh_headings() if h["mesh_id"]}
    absent = [m for m in requested if m not in index]
    if absent:
        raise ValueError(f"MeSH IDs with no drug_indication rows in {client.release}: {absent}")
    info: dict[str, Any] = {"mesh_ids_requested": requested}
    used = set(requested)
    if include_descendants:
        if mesh_client is None:
            raise ValueError("include_descendants=True needs a MeshClient")
        desc_all: dict[str, str] = {}
        for m in requested:
            for d in mesh_client.descendants(m):
                desc_all[d["mesh_id"]] = d["mesh_heading"]
        added = sorted(d for d in desc_all if d in index and d not in used)
        info["mesh_descendants_added"] = {d: index[d] for d in added}
        info["mesh_descendants_without_rows"] = sorted(d for d in desc_all if d not in index)
        used |= set(added)
    used_sorted = sorted(used)
    return used_sorted, {m: index[m] for m in used_sorted}, info


# ------------------------------------------------------------------------ cohorts
def _ref_types(row: Mapping) -> set[str]:
    return {r.get("ref_type") or "" for r in (row.get("indication_refs") or [])} - {""}


def build_cohort(
    mesh_ids: Iterable[str],
    client: ChemblClient,
    *,
    min_phase: float | None = 4.0,
    evidence_sources: Iterable[str] | None = None,
    exclude_sole_evidence: Iterable[str] | None = None,
    include_descendants: bool = False,
    mesh_client: MeshClient | None = None,
    disease_label: str | None = None,
) -> Cohort:
    """Controlled cohort for an explicit list of confirmed MeSH descriptor IDs.

    Parameters
    ----------
    mesh_ids:
        MeSH descriptor IDs (``D######``) present in ChEMBL ``drug_indication``.
    min_phase:
        Keep parents whose indication phase (max over the family, see module docstring)
        is ``>= min_phase``.  ``None`` disables the phase rule (unknown phases kept).
    evidence_sources:
        If given, an indication row qualifies only if at least one of its
        ``indication_refs`` has a ``ref_type`` in this set.
    exclude_sole_evidence:
        If given, an indication row does not qualify when *all* its ref types are in this
        set, e.g. ``{"DailyMed"}`` drops indications supported only by text-mined label
        sections.  Rows with no refs count as sole-evidence rows.
    include_descendants:
        Add MeSH tree descendants of each ID that have rows in this release (NLM SPARQL;
        requires ``mesh_client``).  Off by default because MeSH is a polyhierarchy.
    """
    if min_phase is not None:
        min_phase = float(min_phase)
    used, headings, mesh_info = _resolve_mesh_ids(mesh_ids, client, include_descendants, mesh_client)
    allow = set(evidence_sources) if evidence_sources is not None else None
    sole = set(exclude_sole_evidence) if exclude_sole_evidence is not None else None
    rows = client.drug_indications(used)

    members: dict[str, dict] = {}
    fallback: dict[str, str] = {}
    for r in rows:
        cid = r["molecule_chembl_id"]
        types = _ref_types(r)
        qualifies = (allow is None or bool(types & allow)) and (sole is None or not types <= sole)
        m = members.setdefault(cid, {"phase": math.nan, "qualifying": False, "ref_types": set()})
        m["ref_types"] |= types
        if qualifies:
            m["qualifying"] = True
            p = _phase(r.get("max_phase_for_ind"))
            if not math.isnan(p) and (math.isnan(m["phase"]) or p > m["phase"]):
                m["phase"] = p
        if r.get("parent_molecule_chembl_id"):
            fallback[cid] = r["parent_molecule_chembl_id"]

    records, parent_of = _resolve_parents(client, members, fallback)
    families: dict[str, dict] = {}
    for cid, m in members.items():
        families.setdefault(parent_of[cid], {})[cid] = m
    n_mismatch = sum(1 for c in members if c in fallback and fallback[c] != parent_of[c])

    label = disease_label or " | ".join(headings[m] for m in used)
    mesh_col = ";".join(used)
    out = _curate_families(
        families,
        records,
        label=label,
        mesh_col=mesh_col,
        release=client.release,
        min_phase=min_phase,
        use_phase=True,
    )
    client.verify_release()
    return Cohort(
        kind="controlled",
        disease_label=label,
        release=client.release,
        release_date=client.status().get("chembl_release_date", ""),
        parameters={
            "mesh_ids": used,
            "min_phase": min_phase,
            "evidence_sources": sorted(allow) if allow is not None else None,
            "exclude_sole_evidence": sorted(sole) if sole is not None else None,
            "include_descendants": include_descendants,
            "molecule_type": "Small molecule",
            "dedup": "standard_inchi_key, lowest parent ChEMBL ID kept",
        },
        manifest=_sort_manifest(out),
        extra={
            "mesh_headings": headings,
            **mesh_info,
            "indication_rows": len(rows),
            "parent_mismatch_vs_drug_indication": n_mismatch,
        },
    )


_ATC_PREFIX = re.compile(r"[A-Z](\d{2}([A-Z]([A-Z](\d{2})?)?)?)?")
#: ATC level by code length: 1 -> level1 (anatomical main group), 3 -> level2, 4 -> level3,
#: 5 -> level4, 7 -> level5 (chemical substance).
_ATC_LEVEL = {1: "level1", 3: "level2", 4: "level3", 5: "level4", 7: "level5"}


def build_atc_cohort(
    atc_prefix: str | Iterable[str], client: ChemblClient, *, disease_label: str | None = None
) -> Cohort:
    """Cohort of molecules carrying a WHO ATC level-5 code starting with ``atc_prefix``.

    ``atc_prefix`` is one code prefix (``'A10'``) or several (``('C02', 'C03', 'C07')``);
    a molecule is a candidate if any of its codes starts with any prefix.  ATC codes in
    ChEMBL may sit on a salt record; every molecule is collapsed onto its parent exactly as
    in :func:`build_cohort`.  No phase rule is applied (an ATC code implies a marketed
    product somewhere), so ``phase_for_indication`` is empty and ``indication_ref_types``
    lists the matched codes as ``ATC:<code>``.
    """
    prefixes = [atc_prefix] if isinstance(atc_prefix, str) else list(atc_prefix)
    prefixes = sorted({p.strip().upper() for p in prefixes})
    if not prefixes:
        raise ValueError("atc_prefix is empty")
    for p in prefixes:
        if not _ATC_PREFIX.fullmatch(p):
            raise ValueError(f"not an ATC code prefix: {p!r}")
    members: dict[str, dict] = {}
    desc: dict[str, list[str]] = {}
    level5: set[str] = set()
    for prefix in prefixes:
        for rec in client.molecules_by_atc(prefix):
            codes = {c for c in (rec.get("atc_classifications") or []) if c.startswith(prefix)}
            m = members.setdefault(
                rec["molecule_chembl_id"], {"phase": math.nan, "qualifying": True, "ref_types": set()}
            )
            m["ref_types"] |= {f"ATC:{c}" for c in codes}
        classes = client.atc_classes(prefix)
        level5 |= {c["level5"] for c in classes}
        level = _ATC_LEVEL.get(len(prefix))
        name_key = "who_name" if level == "level5" else f"{level}_description"
        desc[prefix] = sorted({c.get(name_key) or "" for c in classes if level and c.get(level) == prefix} - {""})
    records, parent_of = _resolve_parents(client, members, {})
    families: dict[str, dict] = {}
    for cid, m in members.items():
        families.setdefault(parent_of[cid], {})[cid] = m
    label = disease_label or "ATC " + "+".join(
        p + (f" ({desc[p][0]})" if desc[p] else "") for p in prefixes
    )
    out = _curate_families(
        families, records, label=label, mesh_col="", release=client.release, min_phase=None, use_phase=False
    )
    client.verify_release()
    return Cohort(
        kind="atc",
        disease_label=label,
        release=client.release,
        release_date=client.status().get("chembl_release_date", ""),
        parameters={
            "atc_prefixes": prefixes,
            "molecule_type": "Small molecule",
            "dedup": "standard_inchi_key, lowest parent ChEMBL ID kept",
        },
        manifest=_sort_manifest(out),
        extra={"atc_description": desc, "atc_level5_codes": sorted(level5)},
    )


def naive_cohort(free_text: str, client: ChemblClient) -> Cohort:
    """The v0.x retrieval (``fetcher.fetch_chembl_drugs``) frozen for Experiment A.

    ``drug_indication`` rows with ``mesh_heading__icontains=free_text``; every distinct
    ``molecule_chembl_id`` (salt forms kept as separate molecules, all phases, all molecule
    types) whose ChEMBL ``canonical_smiles`` parses in RDKit is included.  Differences
    from v0.x, all removing nondeterminism: no ``max_drugs`` truncation of an unordered
    set and no early stop at ``2 * max_drugs`` IDs; rows sorted by ChEMBL ID.  The v0.x
    keyword fallbacks (``drug_indication.search``) are omitted because the endpoint is not
    searchable (HTTP 404 / client AssertionError), so they never returned rows.
    ``parent_chembl_id`` is reported for analysis only and is not used for curation.
    """
    if not free_text or not free_text.strip():
        raise ValueError("free_text is empty (it would match every drug_indication row)")
    rows = client.drug_indications_icontains(free_text)
    info: dict[str, dict] = {}
    headings: set[tuple[str, str]] = set()
    for r in rows:
        cid = r["molecule_chembl_id"]
        m = info.setdefault(cid, {"phase": math.nan, "ref_types": set()})
        m["ref_types"] |= _ref_types(r)
        p = _phase(r.get("max_phase_for_ind"))
        if not math.isnan(p) and (math.isnan(m["phase"]) or p > m["phase"]):
            m["phase"] = p
        headings.add((r.get("mesh_id") or "", r.get("mesh_heading") or ""))
    records = client.molecules(info) if info else {}
    mesh_col = ";".join(sorted({h[0] for h in headings if h[0]}))
    out = []
    for cid in sorted(info, key=chembl_sort_key):
        rec = records.get(cid)
        smi = _smiles(rec)
        reason = None
        if rec is None:
            reason = "molecule_record_missing"
        elif not smi:
            reason = "no_smiles"
        else:
            RDLogger.DisableLog("rdApp.*")
            try:
                if Chem.MolFromSmiles(smi) is None:
                    reason = "unparseable_smiles"
            finally:
                RDLogger.EnableLog("rdApp.*")
        out.append(
            {
                "disease_label": free_text,
                "mesh_ids": mesh_col,
                "chembl_release": client.release,
                "molecule_chembl_id": cid,
                "parent_chembl_id": _parent_of(cid, rec, None),
                "pref_name": (rec or {}).get("pref_name") or "",
                "canonical_smiles": smi,
                "standard_inchi_key": _inchikey(rec),
                "molecule_type": (rec or {}).get("molecule_type") or "",
                "phase_for_indication": info[cid]["phase"],
                "indication_ref_types": ";".join(sorted(info[cid]["ref_types"])),
                "decision": "exclude" if reason else "include",
                "reason": reason or "included",
            }
        )
    out.sort(key=lambda r: chembl_sort_key(r["molecule_chembl_id"]))
    df = pd.DataFrame(out, columns=list(MANIFEST_COLUMNS))
    df["phase_for_indication"] = df["phase_for_indication"].astype(float)
    client.verify_release()
    return Cohort(
        kind="naive",
        disease_label=free_text,
        release=client.release,
        release_date=client.status().get("chembl_release_date", ""),
        parameters={"mesh_heading__icontains": free_text, "min_phase": None, "salt_handling": None},
        manifest=df,
        extra={
            "indication_rows": len(rows),
            "mesh_headings_matched": {m: h for m, h in sorted(headings)},
        },
    )


def jaccard(a: Iterable[str], b: Iterable[str]) -> float:
    """``|A n B| / |A u B|`` of two ID sets (NaN if both empty)."""
    sa, sb = set(a) - {""}, set(b) - {""}
    union = sa | sb
    return len(sa & sb) / len(union) if union else math.nan
