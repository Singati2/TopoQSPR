"""End-to-end study pipeline: cohort -> structures -> indices -> features ->
targets -> nested-CV baseline ladder -> results directory + run record.

:func:`run_study` executes one study configuration (see ``docs/PIPELINE.md``
for the schema and every output file):

1. **cohort** - :func:`topoqspr.cohort.build_cohort` from confirmed MeSH IDs
   (ChEMBL, cached per release), or a cohort manifest written earlier, or a
   user SMILES CSV (``id``, ``smiles``);
2. **standardize** - :func:`topoqspr.graph.standardize` on every candidate;
   rejections are kept in ``molecules.csv`` with their code. Two candidates
   whose standardised parents share a full InChIKey are de-duplicated (the
   first in manifest order is kept, the other is excluded as
   ``duplicate_inchikey:<kept id>``), since a duplicate would put the same
   structure in a training and a test fold;
3. **indices** - :func:`topoqspr.indices.compute_indices` (NaN + reason);
4. **features** - the baseline ladder of :mod:`topoqspr.features` plus ``TI``
   and one ``TI_<family>`` set per index family; ``rdkit2d`` is subset per
   target with that target's exclusion rule;
5. **targets** - from :mod:`topoqspr.target_providers` (shipped: RDKit
   computed targets); rows with a missing target value are dropped for that
   target only, and counted;
6. **evaluate** - :func:`topoqspr.evaluate.compare_feature_sets` (repeated
   nested CV, identical outer folds, paired differences against ``size`` with
   corrected repeated-CV t intervals; docs/EVALUATION_PROTOCOL.md) once per
   *learner*: ``selected`` (family and hyperparameters chosen in the inner
   loop among all configured families) and each family on its own;
7. **write** - CSV outputs, ``run_record.json`` last.

Every result CSV is a deterministic function of (config, input data, ChEMBL
release, package versions): no clock, no hash-seed dependence, sorted rows,
``\\n`` line endings. :func:`rerun` re-executes a run record offline from the
ChEMBL cache and checks the CSVs byte for byte.
"""

from __future__ import annotations

import copy
import json
import logging
import math
import os
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np
import pandas as pd
from rdkit import rdBase

from . import __version__
from . import evaluate as ev
from . import features as ft
from . import provenance as pv
from .graph import standardize
from .indices import FAMILIES as INDEX_FAMILIES
from .indices import INDEX_NAMES, REGISTRY, compute_indices
from .models import DEFAULT_FAMILIES, FAMILIES as MODEL_FAMILIES
from .target_providers import TargetSpec, get_target_provider

__all__ = [
    "CONFIG_SCHEMA",
    "DEFAULT_FEATURE_SETS",
    "RESULT_CSVS",
    "StudyResult",
    "load_config",
    "resolve_config",
    "run_study",
    "rerun",
    "compute_index_table",
]

log = logging.getLogger("topoqspr.pipeline")

#: Version of the study-config layout.
CONFIG_SCHEMA = 1

#: Feature sets of EXPERIMENT_PLAN §5. ``TI_family`` expands to one
#: ``TI_<family>`` set per index family that has at least one index.
DEFAULT_FEATURE_SETS: tuple[str, ...] = (
    "null", "size", "size_rings", "composition", "TI", "TI_family", "rdkit2d",
)
_BASE_SETS = ("null", "size", "size_rings", "composition", "TI", "TI_family", "rdkit2d")

#: Result CSVs that a re-run must reproduce byte for byte.
RESULT_CSVS: tuple[str, ...] = (
    "manifest.csv", "molecules.csv", "features.csv", "targets.csv", "failures.csv",
    "folds.csv", "predictions.csv", "paired.csv", "inner.csv", "summary.csv",
)

_DEFAULTS: dict[str, Any] = {
    "schema": CONFIG_SCHEMA,
    "name": "study",
    "cohort": None,  # required
    "standardize": {"min_heavy_atoms": 3},
    "indices": {"names": None},
    "feature_sets": list(DEFAULT_FEATURE_SETS),
    "targets": [{"provider": "rdkit", "names": ["MolWt", "MolLogP", "TPSA", "MolMR"],
                 "options": {}}],
    "evaluation": {
        "learners": None,  # None -> ["selected", *families]
        "families": list(DEFAULT_FAMILIES),
        "outer_splits": 5,
        "seeds": list(ev.DEFAULT_SEEDS),
        "inner_splits": 3,
        "random_state": 0,
        "ad": "knn",
        "split": "random",
        "reference": "size",
        "alpha": 0.05,
        "min_n": 50,
        "min_n_test": 10,
    },
    "output": {"features_format": "csv"},
    "runtime": {"n_jobs": 1, "chembl_cache_dir": None, "offline": False},
}

_COHORT_DEFAULTS: dict[str, dict[str, Any]] = {
    "chembl": {"mesh_ids": None, "min_phase": 4.0, "evidence_sources": None,
               "exclude_sole_evidence": None, "exclude_diagnostics": False,
               "include_descendants": False, "disease_label": None,
               "chembl_release": None},
    "manifest": {"path": None},
    "smiles_csv": {"path": None, "id_column": "id", "smiles_column": "smiles"},
}


# ------------------------------------------------------------------ config
def load_config(path: str | os.PathLike) -> dict:
    """Read a study config from JSON, or YAML when PyYAML is installed."""
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() in (".yaml", ".yml"):
        try:
            import yaml  # type: ignore[import-not-found]
        except ImportError:
            raise RuntimeError("YAML configs need PyYAML, which is not installed; "
                               "use a JSON config") from None
        return yaml.safe_load(text)
    return json.loads(text)


def _merge(defaults: Mapping, given: Mapping, where: str) -> dict:
    unknown = sorted(set(given) - set(defaults))
    if unknown:
        raise ValueError(f"unknown key(s) {unknown} in {where}; allowed: {sorted(defaults)}")
    out = copy.deepcopy(dict(defaults))
    for k, v in given.items():
        out[k] = copy.deepcopy(v)
    return out


def _sorted_or_none(x):
    if x is None:
        return None
    if isinstance(x, str):
        return [x]
    return sorted(str(v) for v in x)


def resolve_config(config: Mapping) -> dict:
    """Validate ``config`` and fill every default (the result is what gets hashed).

    Unknown keys raise, so a typo cannot silently fall back to a default.
    """
    if not isinstance(config, Mapping):
        raise TypeError("config must be a mapping")
    cfg = _merge(_DEFAULTS, config, "config")
    if cfg["schema"] != CONFIG_SCHEMA:
        raise ValueError(f"config schema {cfg['schema']!r} not supported (expected {CONFIG_SCHEMA})")
    for sect in ("standardize", "indices", "evaluation", "output", "runtime"):
        cfg[sect] = _merge(_DEFAULTS[sect], cfg[sect] or {}, sect)

    coh = cfg["cohort"]
    if not isinstance(coh, Mapping) or "source" not in coh:
        raise ValueError("config.cohort must be a mapping with a 'source' "
                         f"({sorted(_COHORT_DEFAULTS)})")
    src = coh["source"]
    if src not in _COHORT_DEFAULTS:
        raise ValueError(f"cohort.source {src!r} not in {sorted(_COHORT_DEFAULTS)}")
    rest = {k: v for k, v in coh.items() if k != "source"}
    coh = {"source": src, **_merge(_COHORT_DEFAULTS[src], rest, f"cohort ({src})")}
    if src == "chembl":
        ids = coh["mesh_ids"]
        if not ids:
            raise ValueError("cohort.mesh_ids is required for source 'chembl'")
        coh["mesh_ids"] = sorted({str(m).strip().upper() for m in ([ids] if isinstance(ids, str) else ids)})
        if coh["min_phase"] is not None:
            coh["min_phase"] = float(coh["min_phase"])
        coh["evidence_sources"] = _sorted_or_none(coh["evidence_sources"])
        coh["exclude_sole_evidence"] = _sorted_or_none(coh["exclude_sole_evidence"])
        coh["exclude_diagnostics"] = bool(coh["exclude_diagnostics"])
        coh["include_descendants"] = bool(coh["include_descendants"])
        if coh["include_descendants"]:
            raise ValueError("include_descendants needs the NLM MeSH client and is not "
                             "supported by the pipeline yet; list the descendant IDs explicitly")
    elif not coh["path"]:
        raise ValueError(f"cohort.path is required for source {src!r}")
    cfg["cohort"] = coh

    cfg["standardize"]["min_heavy_atoms"] = int(cfg["standardize"]["min_heavy_atoms"])
    names = cfg["indices"]["names"]
    if names is not None:
        unknown = sorted(set(names) - set(INDEX_NAMES))
        if unknown:
            raise KeyError(f"unknown index names {unknown}")
        cfg["indices"]["names"] = [n for n in INDEX_NAMES if n in set(names)]

    fsets = list(cfg["feature_sets"])
    bad = [s for s in fsets if s not in _BASE_SETS]
    if bad:
        raise ValueError(f"unknown feature sets {bad}; allowed: {list(_BASE_SETS)}")
    if len(set(fsets)) != len(fsets):
        raise ValueError(f"duplicate feature sets in {fsets}")
    cfg["feature_sets"] = [s for s in _BASE_SETS if s in fsets]  # canonical order

    tg = cfg["targets"]
    if not tg:
        raise ValueError("config.targets is empty")
    cfg["targets"] = [
        {**_merge({"provider": None, "names": None, "options": {}}, t, "targets[]")}
        for t in tg
    ]
    for t in cfg["targets"]:
        if not t["provider"]:
            raise ValueError("every targets[] entry needs a 'provider'")
        t["options"] = dict(t["options"] or {})

    e = cfg["evaluation"]
    fams = list(e["families"])
    unknown = [f for f in fams if f not in MODEL_FAMILIES]
    if unknown:
        raise KeyError(f"unknown model families {unknown}; known: {sorted(MODEL_FAMILIES)}")
    e["families"] = fams
    learners = ["selected", *fams] if e["learners"] is None else list(e["learners"])
    bad = [x for x in learners if x != "selected" and x not in MODEL_FAMILIES]
    if bad:
        raise ValueError(f"unknown learners {bad}: use 'selected' or a model family name")
    if len(set(learners)) != len(learners):
        raise ValueError(f"duplicate learners in {learners}")
    e["learners"] = learners
    e["seeds"] = [int(s) for s in e["seeds"]]
    for k in ("outer_splits", "inner_splits", "random_state", "min_n", "min_n_test"):
        e[k] = int(e[k])
    e["alpha"] = float(e["alpha"])
    if e["split"] not in ("random", "scaffold"):
        raise ValueError("evaluation.split must be 'random' or 'scaffold'")
    if e["ad"] not in (None, "knn", "leverage"):
        raise ValueError("evaluation.ad must be 'knn', 'leverage' or null")
    if e["reference"] not in cfg["feature_sets"] or e["reference"] == "TI_family":
        raise ValueError(f"evaluation.reference {e['reference']!r} must be one of the feature sets")
    if cfg["output"]["features_format"] not in ("csv", "parquet"):
        raise ValueError("output.features_format must be 'csv' or 'parquet'")
    cfg["runtime"]["n_jobs"] = int(cfg["runtime"]["n_jobs"])
    cfg["runtime"]["offline"] = bool(cfg["runtime"]["offline"])
    cfg["name"] = str(cfg["name"])
    # round-trip through canonical JSON: only JSON types survive, hashing is exact
    return json.loads(pv.canonical_json(cfg))


# ----------------------------------------------------------------- helpers
def _write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False, lineterminator="\n", encoding="utf-8")


class _Timer:
    def __init__(self) -> None:
        self.t: dict[str, float] = {}

    def run(self, stage: str, fn: Callable, *a, **kw):
        t0 = time.perf_counter()
        out = fn(*a, **kw)
        self.t[stage] = self.t.get(stage, 0.0) + time.perf_counter() - t0
        log.info("stage %s: %.1f s", stage, self.t[stage])
        return out


@dataclass
class StudyResult:
    """Paths and in-memory tables of a finished study."""

    out_dir: Path
    run_record: dict
    summary: pd.DataFrame
    molecules: pd.DataFrame
    files: list[str] = field(default_factory=list)


# ------------------------------------------------------------------ cohort
def _read_text_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])


def _cohort_stage(cfg: dict, out: Path, base: Path, timestamp_utc: str, client):
    """Write manifest.csv/.json; return (candidates frame, chembl info, inputs)."""
    coh = cfg["cohort"]
    inputs: dict[str, dict] = {}
    if coh["source"] == "chembl":
        from .chembl import ChemblClient
        from .cohort import build_cohort

        if client is None:
            rt = cfg["runtime"]
            client = ChemblClient(rt["chembl_cache_dir"], offline=rt["offline"],
                                  release=coh["chembl_release"])
        if coh["chembl_release"] and client.release != coh["chembl_release"]:
            raise RuntimeError(f"config pins {coh['chembl_release']} but the client serves "
                               f"{client.release}")
        cohort = build_cohort(
            coh["mesh_ids"], client, min_phase=coh["min_phase"],
            evidence_sources=coh["evidence_sources"],
            exclude_sole_evidence=coh["exclude_sole_evidence"],
            exclude_diagnostics=coh["exclude_diagnostics"],
            include_descendants=False, disease_label=coh["disease_label"])
        cohort.write(out / "manifest.csv", timestamp_utc=timestamp_utc)
        inc = cohort.included
        cand = pd.DataFrame({"id": inc["parent_chembl_id"].astype(str),
                             "input_smiles": inc["canonical_smiles"].astype(str)})
        chembl = dict(release=cohort.release, release_date=cohort.release_date,
                      source="ChEMBL REST API via topoqspr.chembl.ChemblClient",
                      offline=bool(getattr(client, "offline", False)))
        return cand, chembl, inputs

    path = Path(coh["path"])
    path = path if path.is_absolute() else base / path
    if not path.exists():
        raise FileNotFoundError(f"cohort input {path} not found")
    inputs["cohort"] = dict(path=coh["path"], sha256=pv.sha256_file(path))
    if coh["source"] == "manifest":
        from .cohort import MANIFEST_COLUMNS

        df = _read_text_csv(path)
        missing = [c for c in MANIFEST_COLUMNS if c not in df.columns]
        if missing:
            raise ValueError(f"{path} is not a cohort manifest; missing columns {missing}")
        shutil.copyfile(path, out / "manifest.csv")
        side = path.with_suffix(".json")
        if side.exists():
            shutil.copyfile(side, out / "manifest.json")
            inputs["cohort_sidecar"] = dict(path=str(Path(coh["path"]).with_suffix(".json")),
                                            sha256=pv.sha256_file(side))
        inc = df[df["decision"] == "include"]
        cand = pd.DataFrame({"id": inc["parent_chembl_id"].to_numpy(),
                             "input_smiles": inc["canonical_smiles"].to_numpy()})
        releases = sorted(set(df["chembl_release"]) - {""})
        chembl = None
        if releases:
            if len(releases) > 1:
                raise ValueError(f"manifest mixes ChEMBL releases {releases}")
            date = ""
            if side.exists():
                date = json.loads(side.read_text()).get("chembl_release_date", "")
            chembl = dict(release=releases[0], release_date=date,
                          source=f"cohort manifest {coh['path']}", offline=True)
        return cand, chembl, inputs

    # smiles_csv
    df = _read_text_csv(path)
    for c in (coh["id_column"], coh["smiles_column"]):
        if c not in df.columns:
            raise ValueError(f"{path}: column {c!r} not found (columns: {list(df.columns)})")
    ids = df[coh["id_column"]].str.strip()
    if (ids == "").any():
        raise ValueError(f"{path}: empty id in rows {list(np.flatnonzero(ids == ''))}")
    dup = ids[ids.duplicated()].unique().tolist()
    if dup:
        raise ValueError(f"{path}: duplicate ids {dup[:10]}")
    cand = pd.DataFrame({"id": ids.to_numpy(), "input_smiles": df[coh["smiles_column"]].to_numpy()})
    _write_csv(cand.rename(columns={"input_smiles": "smiles"}), out / "manifest.csv")
    (out / "manifest.json").write_text(json.dumps(dict(
        kind="smiles_csv", source=coh["path"], input_sha256=inputs["cohort"]["sha256"],
        rows=int(len(cand)), id_column=coh["id_column"], smiles_column=coh["smiles_column"],
        topoqspr_version=__version__, created_utc=timestamp_utc,
        manifest_sha256=pv.sha256_file(out / "manifest.csv")), sort_keys=True, indent=1) + "\n")
    return cand, None, inputs


# ------------------------------------------------------------ standardize
def _standardize_stage(cand: pd.DataFrame, min_heavy_atoms: int):
    rows, mols, graphs = [], [], []
    seen: dict[str, str] = {}
    for rid, smi in zip(cand["id"], cand["input_smiles"]):
        res = standardize(smi, min_heavy_atoms=min_heavy_atoms)
        row = dict(id=rid, input_smiles=smi, status="ok", reason="", detail="",
                   smiles=res.canonical_smiles or "", inchikey=res.inchikey or "",
                   parent_note=res.parent_note or "", n_heavy=math.nan, n_edges=math.nan)
        if not res.ok:
            row.update(status="rejected", reason=res.rejection.code, detail=res.rejection.detail)
        elif res.inchikey and res.inchikey in seen:
            row.update(status="rejected", reason=f"duplicate_inchikey:{seen[res.inchikey]}",
                       detail="standardised parent identical to an earlier molecule")
        else:
            if res.inchikey:
                seen[res.inchikey] = rid
            row.update(n_heavy=res.graph.n, n_edges=res.graph.m)
            mols.append(res.mol)
            graphs.append(res.graph)
        rows.append(row)
    table = pd.DataFrame(rows, columns=["id", "input_smiles", "status", "reason", "detail",
                                        "smiles", "inchikey", "parent_note", "n_heavy", "n_edges"])
    table["n_heavy"] = table["n_heavy"].astype("Int64")
    table["n_edges"] = table["n_edges"].astype("Int64")
    return table, mols, graphs


# --------------------------------------------------------------- indices
def compute_index_table(graphs: Sequence, index, names: Sequence[str] | None = None,
                        exact: bool = False) -> tuple[pd.DataFrame, list[dict]]:
    """Index values (one row per graph) and ``{"row","feature","reason"}`` failures.

    Undefined indices are NaN with the registry's reason. An unexpected
    exception makes every index of that molecule NaN with the exception as
    reason (never 0). ``exact=True`` keeps Python ints (Hosoya Z can exceed
    2**53) in object columns; otherwise values are floats.
    """
    cols = list(INDEX_NAMES) if names is None else list(names)
    rows, failures = [], []
    for label, g in zip(index, graphs):
        try:
            vals, reasons = compute_indices(g, cols)
        except Exception as exc:  # noqa: BLE001 - recorded, never imputed
            rows.append({c: math.nan for c in cols})
            failures += [{"row": label, "feature": c, "reason": f"exception: {exc!r}"[:300]}
                         for c in cols]
            continue
        row = {}
        for c in cols:
            v = vals[c]
            if not exact:
                try:
                    v = float(v)
                except OverflowError:
                    failures.append({"row": label, "feature": c,
                                     "reason": "value exceeds float range"})
                    v = math.nan
            row[c] = v
        for c, why in reasons.items():
            failures.append({"row": label, "feature": c, "reason": why})
        rows.append(row)
    df = pd.DataFrame(rows, index=pd.Index(index), columns=cols,
                      dtype=None if exact else float)
    return df, failures


# ------------------------------------------------------------- features
def _feature_stage(cfg, ids, mols, graphs):
    """Wide feature table (prefixed columns), set -> columns map, failures."""
    idx = pd.Index(ids, name="id")
    blocks, failures, sets = [], [], {}
    fsets = cfg["feature_sets"]
    sr = ft.size_rings_features(mols, idx)
    blocks.append(sr.add_prefix("size."))
    failures += [dict(block="size", **f) for f in sr.attrs.get("failures", [])]
    sets["null"] = []
    sets["size"] = [f"size.{c}" for c in ft.SIZE_COLUMNS]
    sets["size_rings"] = [f"size.{c}" for c in ft.SIZE_RINGS_COLUMNS]
    comp = ft.composition_features(mols, idx)
    blocks.append(comp.add_prefix("composition."))
    failures += [dict(block="composition", **f) for f in comp.attrs.get("failures", [])]
    sets["composition"] = [f"composition.{c}" for c in ft.COMPOSITION_COLUMNS]
    if "TI" in fsets or "TI_family" in fsets:
        names = cfg["indices"]["names"]
        ti, tfail = compute_index_table(graphs, idx, names)
        blocks.append(ti.add_prefix("ti."))
        failures += [dict(block="ti", **f) for f in tfail]
        cols = list(ti.columns)
        sets["TI"] = [f"ti.{c}" for c in cols]
        fam_of = {s.name: s.family for s in REGISTRY}
        for fam in INDEX_FAMILIES:
            members = [f"ti.{c}" for c in cols if fam_of[c] == fam]
            if members:
                sets[f"TI_{fam}"] = members
    if "rdkit2d" in fsets:
        rd = ft.rdkit2d_features(mols, None, idx)
        blocks.append(rd.add_prefix("rdkit2d."))
        failures += [dict(block="rdkit2d", **f) for f in rd.attrs.get("failures", [])]
    wide = pd.concat(blocks, axis=1)
    wide.index.name = "id"
    selected = {}
    for s in fsets:
        if s == "TI_family":
            selected.update({k: v for k, v in sets.items() if k.startswith("TI_")})
        elif s != "rdkit2d":
            selected[s] = sets[s]
        else:
            selected["rdkit2d"] = None  # per target
    return wide, selected, failures


def _sets_for_target(wide: pd.DataFrame, sets: Mapping, spec: TargetSpec, rows) -> dict:
    out = {}
    for name, cols in sets.items():
        if cols is None:  # rdkit2d with this target's exclusions
            cols = [f"rdkit2d.{c}" for c in ft.rdkit2d_names(spec.rdkit2d_exclusion)]
        out[name] = wide.iloc[rows][cols] if cols else pd.DataFrame(index=wide.index[rows])
    return out


# --------------------------------------------------------------- targets
def _target_stage(cfg, mol_table, mols):
    specs: list[TargetSpec] = []
    frames, failures, prov = [], [], []
    for entry in cfg["targets"]:
        provider = get_target_provider(entry["provider"], entry["names"], entry["options"])
        sp = provider.specs()
        vals = provider.compute(mol_table, mols)
        got = list(vals.values.columns)
        if got != [s.name for s in sp] or not vals.values.index.equals(mol_table.index):
            raise ValueError(f"target provider {entry['provider']!r} returned columns {got} / "
                             "an index that does not match its specs and the molecule table")
        specs += sp
        frames.append(vals.values.astype(float))
        failures += [dict(block=f"target:{entry['provider']}", **f) for f in vals.failures]
        prov.append(provider.provenance())
    names = [s.name for s in specs]
    if len(set(names)) != len(names):
        raise ValueError(f"duplicate target names across providers: {names}")
    values = pd.concat(frames, axis=1)
    return specs, values, failures, prov


# -------------------------------------------------------------- evaluate
def _evaluate_stage(cfg, wide, sets, specs, values, mols, timer):
    e = cfg["evaluation"]
    groups_all = None
    if e["split"] == "scaffold":
        groups_all = ev.scaffold_groups(mols)
    folds, preds, paired, inner, summ = [], [], [], [], []
    n_jobs = cfg["runtime"]["n_jobs"]
    for spec in specs:
        y = values[spec.name]
        keep = np.isfinite(y.to_numpy(dtype=float))
        rows = np.flatnonzero(keep)
        y_t = y.iloc[rows]
        groups = None if groups_all is None else groups_all[rows]
        sets_t = _sets_for_target(wide, sets, spec, rows)
        for learner in e["learners"]:
            fams = tuple(e["families"]) if learner == "selected" else (learner,)
            t0 = time.perf_counter()
            cmp = ev.compare_feature_sets(
                sets_t, y_t, reference=e["reference"], families=fams, alpha=e["alpha"],
                outer_splits=e["outer_splits"], seeds=e["seeds"],
                inner_splits=e["inner_splits"], groups=groups,
                random_state=e["random_state"], ad=e["ad"], n_jobs=n_jobs,
                min_n=e["min_n"], min_n_test=e["min_n_test"])
            dt = time.perf_counter() - t0
            timer.t[f"evaluate:{spec.name}:{learner}"] = dt
            log.info("evaluate %s / %s: %.1f s", spec.name, learner, dt)
            key = dict(target=spec.name, learner=learner)
            for set_name, res in cmp.results.items():
                f = res.folds.copy()
                f.insert(0, "set", set_name)
                folds.append(_prefix(f, key))
                p = res.predictions.copy()
                p.insert(0, "set", set_name)
                preds.append(_prefix(p, key))
                inn = res.inner.copy()
                inn.insert(0, "set", set_name)
                inner.append(_prefix(inn, key))
            pr = cmp.paired[["set", "repeat", "fold", "d_r2", "d_mae", "d_rmse",
                             "n_train", "n_test"]].copy()
            paired.append(_prefix(pr, key))
            s = cmp.summary.copy()
            s["family_counts"] = [json.dumps(cmp.results[n].summary()["family_counts"],
                                             sort_keys=True) for n in s["set"]]
            s["n"] = [cmp.results[n].meta["n"] for n in s["set"]]
            s["min_n_test"] = [cmp.results[n].meta["min_n_test"] for n in s["set"]]
            s["underpowered_reason"] = [cmp.results[n].meta["underpowered_reason"]
                                        for n in s["set"]]
            s["reference"] = e["reference"]
            s["target_class"] = spec.target_class
            s["n_target_missing"] = int((~keep).sum())
            summ.append(_prefix(s, key))
    cat = lambda xs: pd.concat(xs, ignore_index=True) if xs else pd.DataFrame()  # noqa: E731
    summary = cat(summ)
    lead = ["target", "target_class", "learner", "set", "reference", "n", "n_target_missing", "p"]
    summary = summary[lead + [c for c in summary.columns if c not in lead]]
    return dict(folds=cat(folds), predictions=cat(preds), paired=cat(paired),
                inner=cat(inner), summary=summary)


def _prefix(df: pd.DataFrame, key: Mapping) -> pd.DataFrame:
    for i, (k, v) in enumerate(key.items()):
        df.insert(i, k, v)
    return df


# ------------------------------------------------------------------- run
def run_study(config: Mapping, out_dir: str | os.PathLike, *, timestamp_utc: str,
              base_dir: str | os.PathLike | None = None, client=None,
              overwrite: bool = False) -> StudyResult:
    """Run a study and write its results directory (see module docstring).

    Parameters
    ----------
    config : study configuration (validated and completed by :func:`resolve_config`).
    out_dir : results directory; must be empty or absent unless ``overwrite``.
    timestamp_utc : ISO-8601 time stamp recorded in the manifest sidecar and
        the run record (passed in, never read from the clock here).
    base_dir : directory that relative input paths in the config refer to
        (default: the current directory); recorded in the run record.
    client : optional :class:`topoqspr.chembl.ChemblClient` (e.g. offline,
        pinned to a release); built from ``config.runtime`` otherwise.
    """
    code = pv.git_info()  # the source tree as imported, before anything can change it
    cfg = resolve_config(config)
    out = Path(out_dir)
    if out.exists() and any(out.iterdir()):
        if not overwrite:
            raise FileExistsError(f"{out} is not empty (pass overwrite=True to replace)")
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    base = Path(base_dir).resolve() if base_dir is not None else Path.cwd().resolve()
    timer = _Timer()

    cand, chembl, inputs = timer.run("cohort", _cohort_stage, cfg, out, base, timestamp_utc,
                                     client)
    with rdBase.BlockLogs():
        mol_table_all, mols, graphs = timer.run(
            "standardize", _standardize_stage, cand, cfg["standardize"]["min_heavy_atoms"])
    ok = mol_table_all[mol_table_all["status"] == "ok"]
    ids = ok["id"].tolist()
    log.info("cohort: %d candidates, %d standardised", len(cand), len(ids))
    if len(ids) < 2:
        raise ValueError(f"only {len(ids)} usable molecules; nothing to evaluate")
    with rdBase.BlockLogs():
        wide, sets, ffail = timer.run("features", _feature_stage, cfg, ids, mols, graphs)
    mol_table = pd.DataFrame({"input_smiles": ok["input_smiles"].to_numpy(),
                              "smiles": ok["smiles"].to_numpy(),
                              "inchikey": ok["inchikey"].to_numpy()},
                             index=pd.Index(ids, name="id"))
    with rdBase.BlockLogs():
        specs, tvalues, tfail, tprov = timer.run("targets", _target_stage, cfg, mol_table, mols)
    t0 = time.perf_counter()
    res = _evaluate_stage(cfg, wide, sets, specs, tvalues, mols, timer)
    timer.t["evaluate"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    _write_csv(mol_table_all, out / "molecules.csv")
    feats = wide.reset_index()
    if cfg["output"]["features_format"] == "parquet":
        feats.to_parquet(out / "features.parquet", index=False)
    _write_csv(feats, out / "features.csv")
    _write_csv(tvalues.reset_index(), out / "targets.csv")
    fail = pd.DataFrame(ffail + tfail, columns=["block", "row", "feature", "reason"])
    _write_csv(fail, out / "failures.csv")
    for name in ("folds", "predictions", "paired", "inner", "summary"):
        _write_csv(res[name], out / f"{name}.csv")
    set_cols = {k: v for k, v in sets.items() if v is not None}
    set_cols["rdkit2d"] = {s.name: [f"rdkit2d.{c}" for c in ft.rdkit2d_names(s.rdkit2d_exclusion)]
                           for s in specs} if "rdkit2d" in sets else None
    meta = dict(feature_sets={k: v for k, v in set_cols.items() if v is not None},
                targets=[s.as_dict() for s in specs], target_providers=tprov)
    (out / "study.json").write_text(json.dumps(meta, sort_keys=True, indent=1) + "\n")
    timer.t["write"] = time.perf_counter() - t0

    written = sorted(p.name for p in out.iterdir() if p.is_file())
    manifest_sha = pv.sha256_file(out / "manifest.csv")
    rej = mol_table_all.loc[mol_table_all["status"] != "ok", "reason"].str.split(":").str[0]
    counts = dict(
        candidates=int(len(cand)), standardised=int(len(ids)),
        rejected={k: int(v) for k, v in sorted(rej.value_counts().items())},
        targets={s.name: int(np.isfinite(tvalues[s.name].to_numpy(float)).sum()) for s in specs},
        feature_set_sizes={k: len(v) for k, v in set_cols.items() if isinstance(v, list)},
        feature_failures=int(len(fail)),
    )
    e = cfg["evaluation"]
    seeds = dict(outer_seeds=e["seeds"], inner_seed_rule="1000 * outer_seed + fold",
                 estimator_random_state=e["random_state"],
                 outer_splitter="KFold(shuffle=True)" if e["split"] == "random"
                 else "balanced_group_kfold (Bemis-Murcko scaffolds)")
    record = pv.build_run_record(
        config=cfg, timestamp_utc=timestamp_utc, chembl=chembl, seeds=seeds, inputs=inputs,
        manifest_sha256=manifest_sha, outputs=pv.file_records(out, written),
        timings_s=timer.t, counts=counts, base_dir=str(base), code=code,
        extra=dict(result_csvs=[c for c in RESULT_CSVS if c in written]))
    pv.write_run_record(record, out / "run_record.json")
    return StudyResult(out, record, res["summary"], mol_table_all, written + ["run_record.json"])


# ----------------------------------------------------------------- rerun
def rerun(record_path: str | os.PathLike, out_dir: str | os.PathLike, *, timestamp_utc: str,
          cache_dir: str | os.PathLike | None = None, base_dir: str | os.PathLike | None = None,
          n_jobs: int | None = None, offline: bool = True) -> dict:
    """Re-execute a run record and compare its result CSVs byte for byte.

    The ChEMBL client is offline (by default) and pinned to the recorded
    release, so the cohort is rebuilt from the cache of that release only.
    Input files must hash to the recorded values. Returns a report (also
    written to ``<out_dir>/rerun_report.json``) with ``identical`` /
    ``different`` / ``missing`` CSV names and the package versions that
    differ from the recorded environment; ``ok`` is True iff every result
    CSV is byte-identical.
    """
    rec = pv.read_run_record(record_path)
    cfg = copy.deepcopy(rec["config"])
    if n_jobs is not None:
        cfg["runtime"]["n_jobs"] = int(n_jobs)
    if cache_dir is not None:
        cfg["runtime"]["chembl_cache_dir"] = str(cache_dir)
    cfg["runtime"]["offline"] = bool(offline)
    base = Path(base_dir) if base_dir is not None else Path(rec.get("base_dir") or ".")
    for name, info in rec.get("inputs", {}).items():
        p = Path(info["path"])
        p = p if p.is_absolute() else base / p
        if not p.exists():
            raise FileNotFoundError(f"input {name} ({p}) not found; pass base_dir")
        if pv.sha256_file(p) != info["sha256"]:
            raise ValueError(f"input {name} ({p}) differs from the recorded file (sha256)")
    client = None
    if cfg["cohort"]["source"] == "chembl":
        from .chembl import ChemblClient

        client = ChemblClient(cfg["runtime"]["chembl_cache_dir"], offline=cfg["runtime"]["offline"],
                              release=rec["chembl"]["release"])
    res = run_study(cfg, out_dir, timestamp_utc=timestamp_utc, base_dir=base, client=client)
    new = res.run_record
    report = dict(record=str(record_path), config_sha256=rec["config_sha256"],
                  config_sha256_rerun=new["config_sha256"], identical=[], different=[],
                  missing=[])
    for name in rec.get("result_csvs", RESULT_CSVS):
        old = rec["outputs"].get(name)
        cur = new["outputs"].get(name)
        if old is None or cur is None:
            report["missing"].append(name)
        elif old["sha256"] == cur["sha256"]:
            report["identical"].append(name)
        else:
            report["different"].append(name)
    old_pk, new_pk = rec["environment"]["packages"], new["environment"]["packages"]
    report["package_differences"] = {
        k: [old_pk.get(k), new_pk.get(k)] for k in sorted(set(old_pk) | set(new_pk))
        if old_pk.get(k) != new_pk.get(k)}
    report["code"] = dict(recorded=rec["code"].get("commit"), rerun=new["code"].get("commit"),
                          rerun_dirty=new["code"].get("dirty"))
    report["ok"] = (not report["different"] and not report["missing"]
                    and report["config_sha256"] == report["config_sha256_rerun"])
    (Path(out_dir) / "rerun_report.json").write_text(
        json.dumps(report, sort_keys=True, indent=1) + "\n")
    return report
