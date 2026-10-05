"""End-to-end study pipeline: cohort -> structures -> indices -> features ->
targets -> nested-CV baseline ladder -> results directory + run record.

:func:`run_study` executes one study configuration (see ``docs/PIPELINE.md``
for the schema and every output file):

1. **cohort** - :func:`topoqspr.cohort.build_cohort` from confirmed MeSH IDs
   (ChEMBL, cached per release), or the pooled approved-drug set
   (:func:`topoqspr.expdata.approved_drug_pool`), or a cohort manifest
   written earlier, or a user SMILES CSV (``id``, ``smiles``);
2. **standardize** - :func:`topoqspr.graph.standardize` on every candidate;
   rejections are kept in ``molecules.csv`` with their code. Two candidates
   whose standardised parents share a full InChIKey are de-duplicated (the
   first in manifest order is kept, the other is excluded as
   ``duplicate_inchikey:<kept id>``), since a duplicate would put the same
   structure in a training and a test fold. Parents that share only the
   InChIKey first block (stereoisomers, racemate and enantiomer,
   isotopologues: the same graph, so identical 2D features) are all kept and
   labelled with a ``connectivity_group``; the evaluation keeps every such
   group inside one fold;
3. **indices** - :func:`topoqspr.indices.compute_indices` (NaN + reason);
4. **features** - the baseline ladder of :mod:`topoqspr.features` plus ``TI``
   and one ``TI_<family>`` set per index family; ``rdkit2d`` is subset per
   target with that target's exclusion rule. The ``null``, ``size`` and
   ``composition`` baselines and the ``size`` reference are mandatory
   (:func:`resolve_config` rejects a config without them unless
   ``evaluation.allow_custom_baselines`` is true, which the run record shows);
5. **targets** - from :mod:`topoqspr.target_providers` (shipped: RDKit
   computed targets); rows with a missing target value are dropped for that
   target only, and counted;
6. **evaluate** - :func:`topoqspr.evaluate.compare_feature_sets` (repeated
   nested CV, identical outer folds, paired differences against ``size`` with
   corrected repeated-CV t intervals; docs/EVALUATION_PROTOCOL.md) once per
   *learner*: ``selected`` (family and hyperparameters chosen in the inner
   loop among all configured families) and each family on its own;
7. **protocol** - the confirmatory protocol of EXPERIMENT_PLAN §7 on the
   confirmatory family named in ``config.protocol.family``: the
   Y-randomised control (:func:`permutation_control`), the whole-procedure
   paired bootstrap (:func:`paired_bootstrap`; by default only for the pooled
   cohort, other cohorts need ``protocol.bootstrap.enabled = true``) and the
   confirmatory table
   (corrected-t CI, BH and BY q, practical-significance flag, AD share,
   permutation p); plus applicability-domain coverage and in-domain metrics
   for every evaluated (target, learner, set) (:func:`ad_report`);
8. **write** - CSV outputs, ``run_record.json`` last. All files are written
   to a staging directory and moved into the results directory only when the
   run has succeeded, so a failed run leaves nothing behind.

Every result CSV is a deterministic function of (config, input data, ChEMBL
release, package versions): no clock, no hash-seed dependence, sorted rows,
``\\n`` line endings. :func:`rerun` re-executes a run record offline from the
ChEMBL cache and checks the CSVs byte for byte.
"""

from __future__ import annotations

import copy
import importlib.util
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
from . import stats
from .graph import standardize
from .indices import FAMILIES as INDEX_FAMILIES
from .chembl import CHEMBL_DATA_LICENCE
from .indices import INDEX_NAMES, REGISTRY, compute_indices
from .models import DEFAULT_FAMILIES, FAMILIES as MODEL_FAMILIES
from .target_providers import TargetSpec, get_target_provider

__all__ = [
    "CONFIG_SCHEMA",
    "DEFAULT_FEATURE_SETS",
    "MANDATORY_BASELINES",
    "MANDATORY_REFERENCE",
    "RESULT_CSVS",
    "StudyResult",
    "load_config",
    "resolve_config",
    "run_study",
    "rerun",
    "compute_index_table",
    "read_results",
    "permutation_control",
    "paired_bootstrap",
    "ad_report",
    "confirmatory_table",
    "permutation_indices",
    "PROTOCOL_CSVS",
    "AD_CSVS",
]

log = logging.getLogger("topoqspr.pipeline")

#: Version of the study-config layout.
CONFIG_SCHEMA = 1

#: Feature sets of EXPERIMENT_PLAN §5. ``TI_family`` expands to one
#: ``TI_<family>`` set per index family that has at least one index.
DEFAULT_FEATURE_SETS: tuple[str, ...] = (
    "null", "size", "size_rings", "composition", "TI", "TI_family", "rdkit2d",
)
_BASE_SETS = ("null", "size", "size_rings", "composition", "TI", "TI_family", "TI_composition", "rdkit2d")

#: Result CSVs that a re-run must reproduce byte for byte.
RESULT_CSVS: tuple[str, ...] = (
    "manifest.csv", "molecules.csv", "features.csv", "targets.csv", "failures.csv",
    "folds.csv", "predictions.csv", "paired.csv", "inner.csv", "summary.csv",
)

#: Applicability-domain report, written for every run with ``evaluation.ad`` set.
AD_CSVS: tuple[str, ...] = ("ad.csv", "ad_summary.csv")

#: Files of the protocol stage (written only when ``config.protocol`` is set).
PROTOCOL_CSVS: tuple[str, ...] = (
    "permutation.csv", "permutation_summary.csv", "bootstrap.csv", "confirmatory.csv",
)

_PROTOCOL_DEFAULTS: dict[str, Any] = {
    # null: no confirmatory analysis; "auto": every class-3 target x (TI vs
    # evaluation.reference); or a list of {"target", "set", "reference"} entries
    # (or "target:SET-REFERENCE" strings)
    "family": None,
    "learner": "selected",          # learner whose results form the confirmatory table
    "practical_threshold": 0.05,    # EXPERIMENT_PLAN §7 minimal meaningful Delta R2
    "q": 0.05,                      # FDR level for BH / BY
    "permutation": {"enabled": True, "n_permutations": 200, "seed": 0,
                    "learner": None, "feature_sets": None, "quantile": 0.95},
    # enabled null (default): on for cohort.source "pooled" only (EXPERIMENT_PLAN §7: at
    # per-disease n most resamples are undefined); any other cohort needs an explicit
    # true, which the resolved config (and so the run record) then carries
    "bootstrap": {"enabled": None, "n_boot": 200, "seed": 0, "families": None},
}

#: Baselines every study must evaluate (EXPERIMENT_PLAN §5), and the reference
#: every paired difference is taken against. Overridable only with
#: ``evaluation.allow_custom_baselines = true``, which is kept in the resolved
#: config and written to the run record (``baselines``).
MANDATORY_BASELINES: tuple[str, ...] = ("null", "size", "composition")
MANDATORY_REFERENCE = "size"

#: Default of protocol.bootstrap.enabled (written to run_record.protocol.bootstrap).
ENABLED_RULE_BOOT = ("protocol.bootstrap.enabled null (default) = on for cohort.source "
                     "'pooled' only (EXPERIMENT_PLAN §7); other cohorts run it only with an "
                     "explicit true")

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
        # true lets feature_sets omit a mandatory baseline or use another reference;
        # dropped from the resolved config when false so existing config hashes hold
        "allow_custom_baselines": False,
    },
    # null = no protocol stage (and no "protocol" key in the resolved config, so
    # configs written before the protocol stage keep their config_sha256)
    "protocol": None,
    "output": {"features_format": "csv"},
    "runtime": {"n_jobs": 1, "chembl_cache_dir": None, "offline": False,
                "expdata_cache_dir": None},
}

_COHORT_DEFAULTS: dict[str, dict[str, Any]] = {
    "chembl": {"mesh_ids": None, "min_phase": 4.0, "evidence_sources": None,
               "exclude_sole_evidence": None, "exclude_diagnostics": False,
               "include_descendants": False, "disease_label": None,
               "chembl_release": None},
    # pooled approved small-molecule drugs for any indication (expdata.approved_drug_pool)
    "pooled": {"chembl_release": None},
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
        cfg = yaml.safe_load(text)
    else:
        cfg = json.loads(text)
    if not isinstance(cfg, Mapping):
        raise ValueError(f"{path}: a study config must be a mapping (JSON object), "
                         f"got {type(cfg).__name__}")
    return cfg


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
        raise ValueError(f"config must be a mapping, got {type(config).__name__}")
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
    elif src == "pooled":
        pass
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
    if not tg or isinstance(tg, (str, Mapping)):
        raise ValueError("config.targets must be a non-empty list of provider entries")
    if not all(isinstance(t, Mapping) for t in tg):
        raise ValueError("every targets[] entry must be a mapping")
    cfg["targets"] = [
        {**_merge({"provider": None, "names": None, "options": {}}, t, "targets[]")}
        for t in tg
    ]
    for t in cfg["targets"]:
        if not t["provider"]:
            raise ValueError("every targets[] entry needs a 'provider'")
        if t["names"] is not None:
            names_t = [t["names"]] if isinstance(t["names"], str) else list(t["names"])
            if not names_t:
                raise ValueError(f"targets[] entry {t['provider']!r}: 'names' is empty "
                                 "(null = all targets of the provider)")
            if len(set(names_t)) != len(names_t):
                raise ValueError(f"targets[] entry {t['provider']!r}: duplicate names {names_t}")
            t["names"] = [str(n) for n in names_t]
        t["options"] = dict(t["options"] or {})

    e = cfg["evaluation"]
    fams = list(e["families"])
    unknown = [f for f in fams if f not in MODEL_FAMILIES]
    if unknown:
        raise KeyError(f"unknown model families {unknown}; known: {sorted(MODEL_FAMILIES)}")
    if not fams:
        raise ValueError("evaluation.families is empty")
    if len(set(fams)) != len(fams):
        raise ValueError(f"duplicate model families in {fams}")
    e["families"] = fams
    learners = ["selected", *fams] if e["learners"] is None else list(e["learners"])
    bad = [x for x in learners if x != "selected" and x not in MODEL_FAMILIES]
    if bad:
        raise ValueError(f"unknown learners {bad}: use 'selected' or a model family name")
    if len(set(learners)) != len(learners):
        raise ValueError(f"duplicate learners in {learners}")
    if not learners:
        raise ValueError("evaluation.learners is empty (null = ['selected', *families])")
    e["learners"] = learners
    seeds = e["seeds"]
    if isinstance(seeds, (str, int)) or not seeds:
        raise ValueError("evaluation.seeds must be a non-empty list of integers")
    e["seeds"] = [int(s) for s in seeds]
    if len(set(e["seeds"])) != len(e["seeds"]):
        # identical repeats would be counted as independent and narrow every CI
        raise ValueError(f"duplicate outer seeds in {e['seeds']}")
    for k in ("outer_splits", "inner_splits", "random_state", "min_n", "min_n_test"):
        e[k] = int(e[k])
    for k in ("outer_splits", "inner_splits"):
        if e[k] < 2:
            raise ValueError(f"evaluation.{k} must be >= 2, got {e[k]}")
    e["alpha"] = float(e["alpha"])
    if e["split"] not in ("random", "scaffold"):
        raise ValueError("evaluation.split must be 'random' or 'scaffold'")
    if e["ad"] not in (None, "knn", "leverage"):
        raise ValueError("evaluation.ad must be 'knn', 'leverage' or null")
    if e["reference"] not in cfg["feature_sets"] or e["reference"] == "TI_family":
        raise ValueError(f"evaluation.reference {e['reference']!r} must be one of the feature sets")
    if not isinstance(e["allow_custom_baselines"], bool):
        raise ValueError("evaluation.allow_custom_baselines must be true or false")
    if e["allow_custom_baselines"]:
        e["allow_custom_baselines"] = True
    else:
        missing = [b for b in MANDATORY_BASELINES if b not in cfg["feature_sets"]]
        if missing:
            raise ValueError(f"feature_sets {cfg['feature_sets']} lack the mandatory baseline(s) "
                             f"{missing} (EXPERIMENT_PLAN §5); set "
                             "evaluation.allow_custom_baselines = true to override (recorded "
                             "in the run record)")
        if e["reference"] != MANDATORY_REFERENCE:
            raise ValueError(f"evaluation.reference {e['reference']!r}: the mandatory reference "
                             f"is {MANDATORY_REFERENCE!r}; set evaluation.allow_custom_baselines "
                             "= true to override (recorded in the run record)")
        del e["allow_custom_baselines"]
    if cfg["output"]["features_format"] not in ("csv", "parquet"):
        raise ValueError("output.features_format must be 'csv' or 'parquet'")
    if cfg["output"]["features_format"] == "parquet" and not _parquet_engine():
        # checked here, not at the write stage, so that a missing engine cannot
        # discard a finished evaluation
        raise ValueError("output.features_format 'parquet' needs pyarrow or fastparquet, "
                         "neither of which is installed (not in requirements-lock.txt); "
                         "use 'csv'")
    cfg["runtime"]["n_jobs"] = int(cfg["runtime"]["n_jobs"])
    cfg["runtime"]["offline"] = bool(cfg["runtime"]["offline"])
    cfg["name"] = str(cfg["name"])
    if cfg["protocol"] is None:
        del cfg["protocol"]
    else:
        cfg["protocol"] = _resolve_protocol(cfg["protocol"], cfg)
    # round-trip through canonical JSON: only JSON types survive, hashing is exact
    return json.loads(pv.canonical_json(cfg))


def _valid_set_name(name: str, fsets: Sequence[str]) -> bool:
    if name == "TI_family":
        return False
    if name in fsets:
        return True
    return ("TI_family" in fsets and name.startswith("TI_")
            and name[3:] in INDEX_FAMILIES)


def _parse_family_entry(entry, reference: str) -> dict:
    """``{"target", "set", "reference"}`` or ``"target:SET-REFERENCE"`` / ``"target:SET"``."""
    if isinstance(entry, str):
        target, sep, comp = entry.partition(":")
        if not sep or not target or not comp:
            raise ValueError(f"protocol.family entry {entry!r}: use 'target:SET-REFERENCE' "
                             "or a mapping with target, set, reference")
        st, _, ref = comp.partition("-")
        entry = {"target": target, "set": st, "reference": ref or reference}
    if not isinstance(entry, Mapping):
        raise ValueError(f"protocol.family entry {entry!r} must be a mapping or a string")
    e = _merge({"target": None, "set": None, "reference": reference}, entry, "protocol.family[]")
    if not e["target"] or not e["set"]:
        raise ValueError(f"protocol.family entry {dict(entry)!r} needs 'target' and 'set'")
    return {k: str(e[k]) for k in ("target", "set", "reference")}


def _resolve_protocol(given, cfg: dict) -> dict:
    """Validate and complete the ``protocol`` section (see :data:`_PROTOCOL_DEFAULTS`)."""
    if not isinstance(given, Mapping):
        raise ValueError("config.protocol must be a mapping or null")
    p = _merge(_PROTOCOL_DEFAULTS, given, "protocol")
    for sub in ("permutation", "bootstrap"):
        p[sub] = _merge(_PROTOCOL_DEFAULTS[sub], p[sub] or {}, f"protocol.{sub}")
    e = cfg["evaluation"]
    fsets = cfg["feature_sets"]
    fam = p["family"]
    if fam is None or fam == "auto":
        pass
    elif isinstance(fam, (str, Mapping)) or not fam:
        raise ValueError("protocol.family must be null, 'auto' or a non-empty list of entries")
    else:
        rows = [_parse_family_entry(x, e["reference"]) for x in fam]
        keys = [(r["target"], r["set"], r["reference"]) for r in rows]
        if len(set(keys)) != len(keys):
            raise ValueError(f"duplicate entries in protocol.family: {keys}")
        for r in rows:
            for k in ("set", "reference"):
                if not _valid_set_name(r[k], fsets):
                    raise ValueError(f"protocol.family: {k} {r[k]!r} is not an evaluated feature "
                                     f"set ({fsets}; TI_<family> needs 'TI_family')")
            if r["reference"] != e["reference"]:
                raise ValueError(f"protocol.family: reference {r['reference']!r} differs from "
                                 f"evaluation.reference {e['reference']!r}; paired differences "
                                 "are computed against evaluation.reference only")
            if r["set"] == r["reference"]:
                raise ValueError(f"protocol.family: set and reference are both {r['set']!r}")
        p["family"] = rows
    if p["learner"] not in e["learners"]:
        raise ValueError(f"protocol.learner {p['learner']!r} is not in evaluation.learners "
                         f"{e['learners']}")
    p["practical_threshold"] = float(p["practical_threshold"])
    p["q"] = float(p["q"])
    if not 0.0 < p["q"] < 1.0:
        raise ValueError("protocol.q must lie in (0, 1)")
    perm = p["permutation"]
    perm["enabled"] = bool(perm["enabled"])
    perm["n_permutations"] = int(perm["n_permutations"])
    perm["seed"] = int(perm["seed"])
    perm["quantile"] = float(perm["quantile"])
    if perm["enabled"] and perm["n_permutations"] < 1:
        raise ValueError("protocol.permutation.n_permutations must be >= 1")
    if perm["learner"] is not None and perm["learner"] not in e["learners"]:
        raise ValueError(f"protocol.permutation.learner {perm['learner']!r} is not in "
                         "evaluation.learners (its observed R2 comes from the main evaluation)")
    if perm["feature_sets"] is not None:
        fs = [perm["feature_sets"]] if isinstance(perm["feature_sets"], str) \
            else list(perm["feature_sets"])
        bad = [s for s in fs if not _valid_set_name(s, fsets)]
        if bad or not fs or len(set(fs)) != len(fs):
            raise ValueError(f"protocol.permutation.feature_sets {fs}: unknown, empty or duplicated")
        perm["feature_sets"] = fs
    boot = p["bootstrap"]
    # null = default: pooled cohort only; other cohorts need an explicit opt-in
    boot["enabled"] = (cfg["cohort"]["source"] == "pooled") if boot["enabled"] is None \
        else bool(boot["enabled"])
    boot["n_boot"] = int(boot["n_boot"])
    boot["seed"] = int(boot["seed"])
    if boot["enabled"] and boot["n_boot"] < 2:
        raise ValueError("protocol.bootstrap.n_boot must be >= 2")
    if boot["families"] is not None:
        bf = [boot["families"]] if isinstance(boot["families"], str) else list(boot["families"])
        bad = [f for f in bf if f not in MODEL_FAMILIES]
        if bad or not bf or len(set(bf)) != len(bf):
            raise ValueError(f"protocol.bootstrap.families {bf}: unknown, empty or duplicated")
        boot["families"] = bf
    return p


# ----------------------------------------------------------------- helpers
def _parquet_engine() -> str | None:
    """Name of an installed pandas parquet engine, or None."""
    for mod in ("pyarrow", "fastparquet"):
        if importlib.util.find_spec(mod) is not None:
            return mod
    return None


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


def read_results(out_dir: str | os.PathLike, name: str = "summary", **kwargs) -> pd.DataFrame:
    """Read ``<out_dir>/<name>.csv`` of a results directory with the right NA rule.

    Missing values are written as empty cells, and only empty cells are NaN.
    A default ``pandas.read_csv`` also turns the strings ``"null"``, ``"NA"``,
    ``"None"``, ``"nan"`` ... into NaN, which silently loses the ``null``
    baseline feature set (``set == "null"``) and could hit IDs and SMILES.
    Extra keyword arguments go to ``pandas.read_csv``.
    """
    path = Path(out_dir) / (name if name.endswith(".csv") else f"{name}.csv")
    kwargs.setdefault("keep_default_na", False)
    kwargs.setdefault("na_values", [""])
    return pd.read_csv(path, **kwargs)


def _cohort_stage(cfg: dict, out: Path, base: Path, timestamp_utc: str, client):
    """Write manifest.csv/.json; return (candidates frame, chembl info, inputs, client).

    ``client`` is the ChEMBL client used for a ``chembl`` or ``pooled`` cohort
    (None otherwise); it is handed on to target providers (ChEMBL-hosted
    measured targets must come from the same release)."""
    coh = cfg["cohort"]
    inputs: dict[str, dict] = {}
    if coh["source"] in ("chembl", "pooled"):
        from .chembl import ChemblClient

        if client is None:
            rt = cfg["runtime"]
            client = ChemblClient(rt["chembl_cache_dir"], offline=rt["offline"],
                                  release=coh["chembl_release"])
        if coh["chembl_release"] and client.release != coh["chembl_release"]:
            raise RuntimeError(f"config pins {coh['chembl_release']} but the client serves "
                               f"{client.release}")
        if coh["source"] == "pooled":
            from .expdata import approved_drug_pool

            cohort = approved_drug_pool(client)
        else:
            from .cohort import build_cohort

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
                      offline=bool(getattr(client, "offline", False)),
                      data_licence=CHEMBL_DATA_LICENCE)
        return cand, chembl, inputs, client

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
        side = path.with_suffix(".json")
        side_meta: dict = {}
        if side.exists():
            side_meta = json.loads(side.read_text(encoding="utf-8"))
            want = side_meta.get("manifest_sha256")
            if want and want != inputs["cohort"]["sha256"]:
                raise ValueError(
                    f"{path} does not match its sidecar {side.name} (manifest_sha256 "
                    f"{want[:12]}..., file {inputs['cohort']['sha256'][:12]}...): the manifest "
                    "was edited after it was built. Rebuild the cohort, or remove the sidecar "
                    "to run the edited manifest as unverified")
        shutil.copyfile(path, out / "manifest.csv")
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
            verified = bool(side_meta.get("manifest_sha256"))
            chembl = dict(release=releases[0],
                          release_date=side_meta.get("chembl_release_date", ""),
                          source=f"cohort manifest {coh['path']}", offline=True,
                          data_licence=CHEMBL_DATA_LICENCE,
                          # True: the CSV hashes to its sidecar's manifest_sha256;
                          # False: no sidecar, so the release label is unverified
                          manifest_verified=verified)
            if not verified:
                log.warning("%s has no sidecar with manifest_sha256; its ChEMBL release "
                            "label %s is recorded as unverified", path, releases[0])
        return cand, chembl, inputs, None

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
    return cand, None, inputs, None


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
    # Stereoisomers, racemate/enantiomer pairs and isotopologues have different
    # full InChIKeys but the same first block (connectivity + H layer), hence the
    # same molecular graph and identical 2D features. They are kept as separate
    # molecules but labelled with a common group (the first member's id) so that
    # the evaluation keeps each group inside one fold.
    first: dict[str, str] = {}
    members: dict[str, int] = {}
    for row in rows:
        if row["status"] == "ok" and row["inchikey"]:
            b1 = row["inchikey"].split("-")[0]
            first.setdefault(b1, row["id"])
            members[b1] = members.get(b1, 0) + 1
    for row in rows:
        b1 = row["inchikey"].split("-")[0] if row["inchikey"] else ""
        row["connectivity_group"] = (first[b1] if row["status"] == "ok" and b1
                                     and members.get(b1, 0) > 1 else "")
    table = pd.DataFrame(rows, columns=["id", "input_smiles", "status", "reason", "detail",
                                        "smiles", "inchikey", "parent_note", "n_heavy", "n_edges",
                                        "connectivity_group"])
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
    if "TI" in fsets or "TI_family" in fsets or "TI_composition" in fsets:
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
    if "TI_composition" in fsets:
        # EXPERIMENT_PLAN §7 secondary: what the indices add on top of element composition
        sets["TI_composition"] = sets["TI"] + sets["composition"]
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
def _make_providers(cfg) -> list[tuple[dict, Any, list[TargetSpec]]]:
    """Instantiate every configured provider and read its specs.

    Runs before the cohort stage, so a provider that cannot be imported or
    built, or target names that collide, fail in seconds, not after hours.
    """
    out, names = [], []
    for entry in cfg["targets"]:
        try:
            provider = get_target_provider(entry["provider"], entry["names"], entry["options"])
        except (ImportError, AttributeError) as exc:
            raise ValueError(f"target provider {entry['provider']!r} could not be loaded: "
                             f"{exc}") from exc
        sp = list(provider.specs())
        if not sp:
            raise ValueError(f"target provider {entry['provider']!r} has no targets")
        names += [s.name for s in sp]
        out.append((entry, provider, sp))
    if len(set(names)) != len(names):
        raise ValueError(f"duplicate target names across providers: {names}")
    return out


def _check_target_contract(label: str, values: pd.DataFrame, failures: list[dict]) -> None:
    """Every missing value must come with a reason; nothing may be +/-inf."""
    arr = values.to_numpy(dtype=float)
    if np.isinf(arr).any():
        r, c = np.argwhere(np.isinf(arr))[0]
        raise ValueError(f"target provider {label!r} returned an infinite value "
                         f"({values.index[r]!r}, {values.columns[c]!r}); return NaN with a "
                         "failure reason instead")
    have = {(str(f.get("row")), str(f.get("feature"))) for f in failures
            if str(f.get("reason") or "").strip()}
    rows, cols = np.nonzero(np.isnan(arr))
    missing = [(values.index[r], values.columns[c]) for r, c in zip(rows, cols)
               if (str(values.index[r]), str(values.columns[c])) not in have]
    if missing:
        raise ValueError(f"target provider {label!r} returned {len(missing)} NaN value(s) "
                         f"without a failure reason, e.g. {missing[:5]}; every missing value "
                         "needs a {'row', 'feature', 'reason'} record (NaN + reason)")


def _target_stage(providers, mol_table, mols):
    specs: list[TargetSpec] = []
    frames, failures, prov = [], [], []
    for entry, provider, sp in providers:
        vals = provider.compute(mol_table, mols)
        got = list(vals.values.columns)
        if got != [s.name for s in sp] or not vals.values.index.equals(mol_table.index):
            raise ValueError(f"target provider {entry['provider']!r} returned columns {got} / "
                             "an index that does not match its specs and the molecule table")
        v = vals.values.astype(float)
        _check_target_contract(entry["provider"], v, vals.failures)
        specs += sp
        frames.append(v)
        failures += [dict(block=f"target:{entry['provider']}", **f) for f in vals.failures]
        prov.append(provider.provenance())
    values = pd.concat(frames, axis=1)
    return specs, values, failures, prov


# -------------------------------------------------------------- evaluate
def _cv_groups(split: str, mols, conn_groups: Sequence[str], ids: Sequence[str]):
    """Outer/inner CV group labels, or None for plain shuffled KFold.

    ``conn_groups`` is the ``connectivity_group`` column of the standardised
    molecules ("" = no other molecule shares its InChIKey first block). Members
    of a connectivity group have identical 2D features and are always kept in
    one fold: in ``scaffold`` mode by giving them their first member's scaffold
    group (this matters for acyclic molecules, which are otherwise singletons);
    in ``random`` mode by switching to grouped folds with one group per
    connectivity group and singletons elsewhere. Without any such group,
    ``random`` stays KFold(shuffle=True).
    """
    conn = [g or i for g, i in zip(conn_groups, ids)]
    if split == "scaffold":
        scaf = ev.scaffold_groups(mols)
        pos = {i: k for k, i in enumerate(ids)}
        return np.asarray([scaf[pos[g]] for g in conn], dtype=object)
    if not any(conn_groups):
        return None
    return np.asarray(conn, dtype=object)


def _too_few_rows(n: int, groups, outer: int, inner: int) -> str:
    """Reason why nested CV cannot be split at all ('' if it can)."""
    units = n if groups is None else len(set(groups))
    if units < outer:
        return f"{units} CV units (molecules or groups) with a value < outer_splits={outer}"
    if units - math.ceil(units / outer) < inner:
        return (f"{units} CV units leave fewer than inner_splits={inner} units in an "
                "outer training part")
    return ""


def _evaluate_stage(cfg, wide, sets, specs, values, groups_all, timer):
    """Run the baseline ladder; returns the result frames plus ``context``:
    per target the rows, y, CV groups, feature sets and the
    :class:`~topoqspr.evaluate.ComparisonResult` of every learner (used by the
    protocol stage), and ``skipped`` targets with the reason."""
    e = cfg["evaluation"]
    folds, preds, paired, inner, summ = [], [], [], [], []
    n_jobs = cfg["runtime"]["n_jobs"]
    context: dict[str, dict] = {}
    skipped: dict[str, str] = {}
    n_values: dict[str, int] = {}
    for spec in specs:
        y = values[spec.name]
        keep = np.isfinite(y.to_numpy(dtype=float))
        rows = np.flatnonzero(keep)
        n_values[spec.name] = int(rows.size)
        y_t = y.iloc[rows]
        groups = None if groups_all is None else groups_all[rows]
        why = _too_few_rows(len(rows), groups, e["outer_splits"], e["inner_splits"])
        if why:
            # a measured target may cover only a handful of a small cohort's drugs:
            # reported (run record counts.skipped_targets), never a crash or a 0
            skipped[spec.name] = why
            log.warning("target %s not evaluated: %s", spec.name, why)
            continue
        sets_t = _sets_for_target(wide, sets, spec, rows)
        context[spec.name] = dict(spec=spec, rows=rows, y=y_t, groups=groups, sets=sets_t,
                                  n_missing=int((~keep).sum()), cmp={})
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
            context[spec.name]["cmp"][learner] = cmp
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
    if len(summary):
        summary = summary[lead + [c for c in summary.columns if c not in lead]]
    else:
        summary = pd.DataFrame(columns=lead)
    return dict(folds=cat(folds), predictions=cat(preds), paired=cat(paired),
                inner=cat(inner), summary=summary, context=context, skipped=skipped,
                n_values=n_values)


def _prefix(df: pd.DataFrame, key: Mapping) -> pd.DataFrame:
    for i, (k, v) in enumerate(key.items()):
        df.insert(i, k, v)
    return df


# --------------------------------------------------------------- protocol
def _learner_families(learner: str, families: Sequence[str]) -> tuple[str, ...]:
    return tuple(families) if learner == "selected" else (learner,)


def _nested_stat(res: ev.NestedCVResult) -> dict:
    """Mean pooled R^2 over repeats (the Y-randomisation statistic; NaN when any
    outer fold failed, as :func:`topoqspr.evaluate.nested_r2_statistic`)."""
    s = res.summary()
    return dict(r2_pooled_mean=s["r2_pooled_mean"] if s["n_failed_folds"] == 0 else math.nan,
                r2_fold_mean=s["r2_fold_mean"], mae_fold_mean=s["mae_fold_mean"],
                n_failed_folds=int(s["n_failed_folds"]), reason=s["reason"])


#: How permutation_control permutes y (written to run_record.protocol.permutation).
PERMUTATION_SCHEME = (
    "row-wise (default_rng([seed, b]).permutation(n)) when the CV is ungrouped; restricted "
    "(block) permutation when CV groups are given: with the same generator, the rows of "
    "singleton groups and of groups whose size no other group shares are permuted row-wise "
    "among themselves (first draw), then, for each group size s >= 2 shared by >= 2 groups "
    "(ascending s; groups in sorted label order), whole groups' y blocks are exchanged by a "
    "random permutation of those groups (row order inside a block kept). Only groups of a "
    "size shared with another group keep their within-group y structure under the null; a "
    "group whose size is unique (e.g. the only 3-member group) has its rows permuted "
    "row-wise with the singletons, so its rows can draw y from different groups. Every "
    "draw is a permutation of y, so the p-value stays valid.")


def permutation_indices(n: int, groups, rng: np.random.Generator) -> np.ndarray:
    """Index vector ``perm`` such that ``y[perm]`` is one Y-randomised target.

    ``groups=None``: ``rng.permutation(n)`` (exchangeable rows). With CV group
    labels the permutation is **restricted** (:data:`PERMUTATION_SCHEME`):
    groups of size >= 2 that share their size with at least one other group
    swap their y values whole, so a target that is shared within such a group
    (stereoisomers, scaffold series) keeps that structure under the null.
    Groups whose size no other group has cannot be exchanged as a block; their
    rows are permuted row-wise together with the singletons, so for them (and
    only them) the within-group structure is *not* kept and the null is the
    row-wise one. The result is always a permutation of ``range(n)``.
    """
    if groups is None:
        return rng.permutation(n)
    g = np.asarray(groups, dtype=object).astype(str)
    if g.size != n:
        raise ValueError(f"groups has {g.size} labels for {n} rows")
    uniq, inv = np.unique(g, return_inverse=True)
    members = [np.flatnonzero(inv == i) for i in range(uniq.size)]
    sizes = np.array([m.size for m in members])
    by_size: dict[int, list[int]] = {}
    for i, sz in enumerate(sizes):
        by_size.setdefault(int(sz), []).append(i)
    perm = np.arange(n)
    rest = [i for sz, gi in by_size.items() if sz == 1 or len(gi) == 1 for i in gi]
    if rest:
        rows = np.sort(np.concatenate([members[i] for i in rest]))
        perm[rows] = rows[rng.permutation(rows.size)]
    for sz in sorted(by_size):
        gi = by_size[sz]
        if sz == 1 or len(gi) == 1:
            continue
        pi = rng.permutation(len(gi))
        for k, j in enumerate(pi):
            perm[members[gi[k]]] = members[gi[j]]
    return perm


def _one_permutation(feature_sets, y, b, seed, groups, families, nested_kwargs):
    perm = permutation_indices(y.size, groups, np.random.default_rng([seed, b]))
    yp = y[perm]
    out = []
    for name, X in feature_sets.items():
        res = ev.nested_cv(X, yp, families, groups=groups, ad=None, **nested_kwargs)
        out.append(dict(set=name, permutation=b, rng=f"default_rng([{seed}, {b}])",
                        **_nested_stat(res)))
    return out


def permutation_control(feature_sets: Mapping[str, pd.DataFrame], y, *,
                        families: Sequence[str], n_permutations: int = 200, seed: int = 0,
                        groups=None, n_jobs: int = 1, **nested_kwargs) -> pd.DataFrame:
    """Y-randomised control (EXPERIMENT_PLAN §7; Ruecker et al. 2007).

    For b = 0..K-1 the target is permuted with
    :func:`permutation_indices` and ``numpy.random.default_rng([seed, b])``:
    ``rng.permutation(n)`` without groups, a restricted block permutation of
    equal-size groups' y values with ``groups`` (the same permutation for
    every feature set, so sets are compared on identical nulls) and the
    **whole** nested procedure of :func:`topoqspr.evaluate.nested_cv` is rerun
    with the same ``families``, outer seeds, inner rule and ``groups``. The
    rows of X and their group labels are not moved, so the outer and inner
    folds are exactly those of the real analysis. The statistic is the mean
    pooled out-of-fold R^2 over repeats (NaN if any outer fold failed).

    ``nested_kwargs`` go to :func:`~topoqspr.evaluate.nested_cv` (outer_splits,
    seeds, inner_splits, random_state, min_n, min_n_test); the AD is not fitted
    (it does not change predictions). Permutations run in parallel with
    ``n_jobs`` (identical results).

    Returns one row per (set, permutation): set, permutation, rng,
    r2_pooled_mean, r2_fold_mean, mae_fold_mean, n_failed_folds, reason.
    """
    y = np.asarray(y, dtype=float).ravel()
    families = tuple(families)
    args = (feature_sets, y)
    if n_jobs == 1:
        outs = [_one_permutation(*args, b, seed, groups, families, nested_kwargs)
                for b in range(n_permutations)]
    else:
        from joblib import Parallel, delayed
        outs = Parallel(n_jobs=n_jobs)(
            delayed(_one_permutation)(*args, b, seed, groups, families, nested_kwargs)
            for b in range(n_permutations))
    rows = [r for o in outs for r in o]
    df = pd.DataFrame(rows, columns=["set", "permutation", "rng", "r2_pooled_mean", "r2_fold_mean",
                                     "mae_fold_mean", "n_failed_folds", "reason"])
    order = {s: i for i, s in enumerate(feature_sets)}
    return df.sort_values(["set", "permutation"], key=lambda c: c.map(order) if c.name == "set"
                          else c, kind="mergesort").reset_index(drop=True)


def _one_resample(feature_sets, y, set_name, reference, units, cv_groups, b, seed, families,
                  nested_kwargs):
    labels = feature_sets[reference].index
    if b < 0:  # the original data, same procedure and families
        idx = np.arange(y.size)
        groups = cv_groups
        rng_txt = "original"
    else:
        idx, _ = stats.group_bootstrap(units, np.random.default_rng([seed, b]))
        base = units if cv_groups is None else cv_groups
        # every copy of an original molecule keeps its original group label, so
        # copies (and connectivity-group partners) never straddle a fold boundary
        groups = np.asarray(base, dtype=object)[idx]
        rng_txt = f"default_rng([{seed}, {b}])"
    new_index = pd.Index([f"{labels[i]}#{k}" for k, i in enumerate(idx)], name=labels.name)
    fs = {nm: feature_sets[nm].iloc[idx].set_axis(new_index, axis=0)
          for nm in (set_name, reference)}
    yb = pd.Series(y[idx], index=new_index)
    n_units = int(len(set(np.asarray(units, dtype=object)[idx].tolist())))
    row = dict(resample=b, rng=rng_txt, n_rows=int(idx.size), n_units=n_units,
               d_r2=math.nan, d_mae=math.nan, r2_set=math.nan, r2_reference=math.nan,
               mae_set=math.nan, mae_reference=math.nan, n_failed_folds=0, reason="")
    why = _too_few_rows(idx.size, groups, nested_kwargs["outer_splits"],
                        nested_kwargs["inner_splits"])
    if why:
        row["reason"] = why
        return row
    cmp = ev.compare_feature_sets(fs, yb, reference=reference, families=families,
                                  groups=groups, ad=None, **nested_kwargs)
    d = cmp.paired[cmp.paired["set"] == set_name]
    s = cmp.summary.set_index("set")
    row.update(d_r2=float(d["d_r2"].mean(skipna=False)), d_mae=float(d["d_mae"].mean(skipna=False)),
               r2_set=float(s.loc[set_name, "r2_fold_mean"]),
               r2_reference=float(s.loc[reference, "r2_fold_mean"]),
               mae_set=float(s.loc[set_name, "mae_fold_mean"]),
               mae_reference=float(s.loc[reference, "mae_fold_mean"]),
               n_failed_folds=int(s["n_failed_folds"].sum()),
               reason="; ".join(dict.fromkeys(r for r in s["reason"].astype(str) if r)))
    return row


def paired_bootstrap(feature_sets: Mapping[str, pd.DataFrame], y, *, set_name: str,
                     reference: str, families: Sequence[str], n_boot: int = 200, seed: int = 0,
                     units=None, cv_groups=None, alpha: float = 0.05, n_jobs: int = 1,
                     **nested_kwargs) -> tuple[pd.DataFrame, dict]:
    """Whole-procedure paired bootstrap of Delta R^2 and Delta MAE (EXPERIMENT_PLAN §7
    sensitivity; CV intervals can undercover, Bates, Hastie & Tibshirani 2023).

    Resample b = 0..B-1 draws the CV *units* with replacement
    (:func:`topoqspr.stats.group_bootstrap` with ``default_rng([seed, b])``):
    ``units`` are the connectivity-group labels (a molecule without a
    stereo/isotope partner is its own unit; default: every row its own unit),
    so stereoisomers are drawn together. On each resample the identical
    repeated nested CV is rerun for ``set_name`` and ``reference`` (same outer
    seeds, inner rule and ``families``) with grouped folds in which every copy
    of an original unit carries its original label (``cv_groups`` if given,
    e.g. scaffold groups, else the unit) - otherwise a duplicated molecule
    would sit in a training and a test fold at once. The statistics are the
    fold means of d_R2 = R2_set - R2_ref and d_MAE = MAE_ref - MAE_set, i.e.
    the same estimands as the corrected-t analysis.

    A row ``resample = -1`` ("original") reruns the procedure once on the
    original data with these families (it equals the main analysis when the
    families are the confirmatory learner's). The percentile interval
    (:func:`topoqspr.stats.percentile_interval`) is NaN with a reason if any
    resample failed. ``families`` may be a faster subset than the main
    learner (cost-aware sensitivity analysis; recorded).

    Returns (rows per resample, summary dict).
    """
    y = np.asarray(y, dtype=float).ravel()
    labels = feature_sets[reference].index
    units = np.asarray([str(i) for i in labels] if units is None else units, dtype=object)
    if units.size != y.size:
        raise ValueError(f"units has {units.size} labels for {y.size} rows")
    families = tuple(families)
    jobs = list(range(-1, n_boot))
    args = (feature_sets, y, set_name, reference, units, cv_groups)
    if n_jobs == 1:
        rows = [_one_resample(*args, b, seed, families, nested_kwargs) for b in jobs]
    else:
        from joblib import Parallel, delayed
        rows = Parallel(n_jobs=n_jobs)(
            delayed(_one_resample)(*args, b, seed, families, nested_kwargs) for b in jobs)
    df = pd.DataFrame(rows)
    boot = df[df["resample"] >= 0]
    ci_r2 = stats.percentile_interval(boot["d_r2"].to_numpy(), alpha)
    ci_mae = stats.percentile_interval(boot["d_mae"].to_numpy(), alpha)
    orig = df[df["resample"] == -1].iloc[0]
    summary = dict(n_boot=int(n_boot), seed=int(seed), families=list(families),
                   d_r2_original=float(orig["d_r2"]), d_mae_original=float(orig["d_mae"]),
                   d_r2_ci_low=ci_r2["ci_low"], d_r2_ci_high=ci_r2["ci_high"],
                   d_r2_boot_mean=ci_r2["mean"], d_r2_boot_sd=ci_r2["sd"],
                   d_mae_ci_low=ci_mae["ci_low"], d_mae_ci_high=ci_mae["ci_high"],
                   n_valid=ci_r2["n_valid"],
                   reason="; ".join(x for x in (ci_r2["reason"], str(orig["reason"])) if x))
    return df, summary


_AD_KEYS = ["target", "learner", "set"]


def ad_report(predictions: pd.DataFrame, folds: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Applicability-domain coverage and in-domain metrics (EXPERIMENT_PLAN §7).

    Per outer fold (``ad``): n_test, n_in_domain, n_outside, share_outside =
    n_outside / n_test (NaN with ``ad_reason`` when the AD was not fitted, e.g.
    the zero-column ``null`` set), and R2/MAE on the in-domain and on the
    outside test rows (secondary; NaN when fewer than 2 / 1 rows or constant
    y). Per (target, learner, set) (``ad_summary``): mean share outside over
    folds, pooled share over all outer-test rows, and for in-domain and
    outside rows the mean over repeats of the pooled per-repeat R2 and MAE,
    next to the all-row pooled values. No NaN is skipped: a fold without an
    AD makes the summary shares NaN and ``reason`` says so.
    """
    keys = _AD_KEYS + ["repeat", "fold"]
    fcols = keys + ["n_test", "ad_reason"]
    if predictions is None or not len(predictions):
        return (pd.DataFrame(columns=fcols), pd.DataFrame(columns=_AD_KEYS))
    reasons = folds.set_index(keys)["ad_reason"]
    rows = []
    for k, g in predictions.groupby(keys, sort=False):
        dom = g["in_domain"].to_numpy(dtype=float)
        yt, yp = g["y_true"].to_numpy(float), g["y_pred"].to_numpy(float)
        n = int(len(g))
        defined = bool(np.all(np.isfinite(dom)))
        m_in = ev._metrics(yt[dom == 1], yp[dom == 1]) if defined else dict(r2=math.nan, mae=math.nan)
        m_out = ev._metrics(yt[dom == 0], yp[dom == 0]) if defined else dict(r2=math.nan, mae=math.nan)
        n_in = int(np.sum(dom == 1)) if defined else math.nan
        n_out = int(np.sum(dom == 0)) if defined else math.nan
        rows.append(dict(zip(keys, k), n_test=n, n_in_domain=n_in, n_outside=n_out,
                         share_outside=(n_out / n) if defined else math.nan,
                         r2_in_domain=m_in["r2"], mae_in_domain=m_in["mae"],
                         r2_outside=m_out["r2"], mae_outside=m_out["mae"],
                         ad_reason=str(reasons.get(k, "")) if not defined else ""))
    ad = pd.DataFrame(rows)
    summ = []
    for k, g in predictions.groupby(_AD_KEYS, sort=False):
        fa = ad[(ad["target"] == k[0]) & (ad["learner"] == k[1]) & (ad["set"] == k[2])]
        n_undef = int(fa["share_outside"].isna().sum())
        rep = []
        for r, gr in g.groupby("repeat", sort=True):
            dom = gr["in_domain"].to_numpy(dtype=float)
            yt, yp = gr["y_true"].to_numpy(float), gr["y_pred"].to_numpy(float)
            ok = bool(np.all(np.isfinite(dom)))
            m_all = ev._metrics(yt, yp)
            m_in = ev._metrics(yt[dom == 1], yp[dom == 1]) if ok else dict(r2=math.nan, mae=math.nan)
            m_out = ev._metrics(yt[dom == 0], yp[dom == 0]) if ok else dict(r2=math.nan, mae=math.nan)
            rep.append((m_all["r2"], m_all["mae"], m_in["r2"], m_in["mae"], m_out["r2"], m_out["mae"]))
        rep = pd.DataFrame(rep, columns=["r2", "mae", "r2_in", "mae_in", "r2_out", "mae_out"])
        dom_all = g["in_domain"].to_numpy(dtype=float)
        pooled = float(np.mean(dom_all == 0)) if np.all(np.isfinite(dom_all)) else math.nan
        why = []
        if n_undef:
            first = fa.loc[fa["share_outside"].isna(), "ad_reason"].iloc[0]
            why.append(f"AD not available in {n_undef} of {len(fa)} folds ({first})")
        summ.append(dict(zip(_AD_KEYS, k), n_folds=int(len(fa)),
                         share_outside_fold_mean=float(fa["share_outside"].mean(skipna=False)),
                         share_outside_fold_max=float(fa["share_outside"].max(skipna=False)),
                         share_outside_pooled=pooled,
                         r2_pooled_mean=float(rep["r2"].mean(skipna=False)),
                         mae_pooled_mean=float(rep["mae"].mean(skipna=False)),
                         r2_in_domain_pooled_mean=float(rep["r2_in"].mean(skipna=False)),
                         mae_in_domain_pooled_mean=float(rep["mae_in"].mean(skipna=False)),
                         r2_outside_pooled_mean=float(rep["r2_out"].mean(skipna=False)),
                         mae_outside_pooled_mean=float(rep["mae_out"].mean(skipna=False)),
                         reason="; ".join(why)))
    return ad, pd.DataFrame(summ)


def _expand_family(protocol: Mapping, specs: Sequence[TargetSpec], reference: str) -> list[dict]:
    fam = protocol.get("family")
    if fam is None:
        return []
    if fam == "auto":
        return [dict(target=s.name, set="TI", reference=reference)
                for s in specs if s.target_class == 3]
    return [dict(r) for r in fam]


CONFIRMATORY_COLUMNS = [
    "target", "target_class", "comparison", "set", "reference", "learner", "n", "underpowered",
    "in_family", "family_m",
    "d_r2_mean", "d_r2_ci_low", "d_r2_ci_high", "d_r2_p", "d_r2_q_bh", "d_r2_q_by",
    "informative_bh", "informative_by",
    "d_mae_mean", "d_mae_ci_low", "d_mae_ci_high", "d_mae_p",
    "practical_threshold", "practically_significant",
    "ad_share_outside_set", "ad_share_outside_reference",
    "perm_k", "perm_p_set", "perm_uninformative_set", "perm_p_reference",
    "perm_uninformative_reference",
    "boot_b", "boot_families", "boot_d_r2_original", "boot_d_r2_ci_low", "boot_d_r2_ci_high",
    "boot_d_mae_original", "boot_d_mae_ci_low", "boot_d_mae_ci_high",
    "reason",
]

#: Membership rule of the multiplicity family (written to run_record.protocol).
FAMILY_RULE = (
    "m = number of pre-specified family members that are not underpowered "
    "(EXPERIMENT_PLAN 7: >= min_n measured values; EVALUATION_PROTOCOL: an underpowered "
    "result carries no ML claim). Underpowered members (n < evaluation.min_n or smallest "
    "outer test fold < evaluation.min_n_test) and members whose target was not evaluated "
    "with fewer than min_n values are left out of the family: q NaN, informative_* and "
    "practically_significant None. Every other member counts in m even when its p is "
    "undefined or it has no result; such a member is never rejected (as p = 1) and keeps "
    "NaN q with the reason, so an evaluation failure cannot loosen the correction.")


def confirmatory_table(family: Sequence[Mapping], summary: pd.DataFrame, ad_summary: pd.DataFrame,
                       *, learner: str, practical_threshold: float = 0.05, q: float = 0.05,
                       perm_summary: pd.DataFrame | None = None,
                       boot: Mapping | None = None, target_classes: Mapping | None = None,
                       skipped: Mapping | None = None, n_values: Mapping | None = None,
                       min_n: int | None = None, permutation_note: str = "permutation not run",
                       bootstrap_note: str = "bootstrap not run",
                       ad_note: str = "AD not fitted") -> pd.DataFrame:
    """One row per (target x comparison) of the confirmatory family.

    From the main analysis (``summary`` of ``learner``): n, mean d_R2 with the
    corrected repeated-CV t CI and p, d_MAE with CI and p.

    **Multiplicity family** (:data:`FAMILY_RULE`): a member is in the family
    (``in_family``) unless its result is underpowered, or its target was not
    evaluated (``skipped``) and has fewer than ``min_n`` measured values
    (``n_values``; unknown n counts as in the family). Excluded members get
    NaN q, None informative_* / practically_significant and the reason
    "underpowered: not in confirmatory family". Across the family, BH and BY
    adjusted q of d_r2_p (:func:`topoqspr.stats.benjamini_hochberg`,
    :func:`~topoqspr.stats.benjamini_yekutieli`) with m = ``family_m`` = the
    number of members in the family, **including** members whose p is NaN or
    that have no result (never rejected; NaN q with the reason).
    ``informative_bh`` = q_bh <= q and d_r2_mean > 0 (EXPERIMENT_PLAN §7);
    ``informative_by`` likewise with BY. ``practically_significant`` =
    d_r2_mean >= ``practical_threshold``.

    AD: mean share of outer-test molecules outside the domain for the set and
    the reference. Permutation: K and the empirical p (+1 corrected) and
    95th-percentile flag for the set and for the reference. Bootstrap: B,
    families, original-data estimate and percentile CIs. Anything unavailable
    is NaN, with ``reason`` (``permutation_note`` / ``bootstrap_note`` /
    ``ad_note`` say why a component was not run at all).
    """
    skipped = dict(skipped or {})
    target_classes = dict(target_classes or {})
    n_values = dict(n_values or {})
    s = summary[summary["learner"] == learner].set_index(["target", "set"]) if len(summary) \
        else pd.DataFrame()
    ads = ad_summary[ad_summary["learner"] == learner].set_index(["target", "set"]) \
        if ad_summary is not None and len(ad_summary) else pd.DataFrame()
    ps = perm_summary.set_index(["target", "set"]) \
        if perm_summary is not None and len(perm_summary) else pd.DataFrame()
    rows = []
    for fam in family:
        t, st, ref = fam["target"], fam["set"], fam["reference"]
        row = {c: math.nan for c in CONFIRMATORY_COLUMNS}
        row.update(target=t, target_class=target_classes.get(t, math.nan),
                   comparison=f"{st}-{ref}", set=st, reference=ref, learner=learner,
                   practical_threshold=practical_threshold, reason="", in_family=True)
        why = []
        evaluated = False
        if t in skipped:
            why.append(f"target not evaluated: {skipped[t]}")
            nv = n_values.get(t)
            if nv is not None:
                row["n"] = int(nv)
                if min_n is not None and int(nv) < int(min_n):
                    row["underpowered"] = True
                    row["in_family"] = False
        elif (t, st) not in s.index:
            why.append(f"no result for set {st!r} (feature set not produced for this run)")
        else:
            evaluated = True
            r = s.loc[(t, st)]
            row.update(n=int(r["n"]), underpowered=bool(r["underpowered"]),
                       **{c: float(r[c]) for c in ("d_r2_mean", "d_r2_ci_low", "d_r2_ci_high",
                                                   "d_r2_p", "d_mae_mean", "d_mae_ci_low",
                                                   "d_mae_ci_high", "d_mae_p")})
            if isinstance(r["reason"], str) and r["reason"]:
                why.append(r["reason"])
            if bool(r["underpowered"]):
                why.append(f"underpowered: {r['underpowered_reason']}")
                row["in_family"] = False
            if not math.isfinite(row["d_r2_p"]):
                why.append("d_r2_p undefined (never rejected"
                           + ("; counted in m)" if row["in_family"] else ")"))
        if not row["in_family"]:
            why.append("underpowered: not in confirmatory family")
        if evaluated:
            for nm, col, sfx in ((st, "ad_share_outside_set", "set"),
                                 (ref, "ad_share_outside_reference", "reference")):
                if (t, nm) in ads.index:
                    row[col] = float(ads.loc[(t, nm), "share_outside_fold_mean"])
                    ar = ads.loc[(t, nm), "reason"]
                    if not math.isfinite(row[col]):
                        why.append(f"AD {sfx}: {ar or 'share undefined'}")
                else:
                    why.append(f"AD {sfx}: {ad_note}")
        for nm, sfx in ((st, "set"), (ref, "reference")):
            if (t, nm) in ps.index:
                pr = ps.loc[(t, nm)]
                row["perm_k"] = int(pr["n_permutations"])
                row[f"perm_p_{sfx}"] = float(pr["p"])
                row[f"perm_uninformative_{sfx}"] = pr["uninformative"]
                if isinstance(pr["reason"], str) and pr["reason"]:
                    why.append(f"permutation {sfx}: {pr['reason']}")
            elif evaluated:
                why.append(f"permutation {sfx} ({nm}): {permutation_note}")
        b = (boot or {}).get((t, st, ref))
        if b is not None:
            row.update(boot_b=b["n_boot"], boot_families=json.dumps(b["families"]),
                       boot_d_r2_original=b["d_r2_original"], boot_d_r2_ci_low=b["d_r2_ci_low"],
                       boot_d_r2_ci_high=b["d_r2_ci_high"], boot_d_mae_original=b["d_mae_original"],
                       boot_d_mae_ci_low=b["d_mae_ci_low"], boot_d_mae_ci_high=b["d_mae_ci_high"])
            if b["reason"]:
                # the bootstrap's original-data row repeats the main analysis's fold
                # reasons: report only what is new, each clause once
                seen = {x for w in why for x in w.split("; ")}
                extra = [x for x in dict.fromkeys(b["reason"].split("; ")) if x and x not in seen]
                if extra:
                    why.append("bootstrap: " + "; ".join(extra))
        elif evaluated:
            why.append(bootstrap_note)
        row["reason"] = "; ".join(dict.fromkeys(why))
        rows.append(row)
    df = pd.DataFrame(rows, columns=CONFIRMATORY_COLUMNS)
    if len(df):
        member = df["in_family"].to_numpy(dtype=bool)
        m = int(member.sum())
        p = df["d_r2_p"].to_numpy(dtype=float)
        p_fam = np.where(member, p, np.nan)
        q_bh, rej_bh = stats.benjamini_hochberg(p_fam, q, m=m)
        q_by, rej_by = stats.benjamini_yekutieli(p_fam, q, m=m)
        df["d_r2_q_bh"], df["d_r2_q_by"] = q_bh, q_by
        df["family_m"] = m
        pos = df["d_r2_mean"].to_numpy(dtype=float) > 0
        defined = np.isfinite(p_fam)
        df["informative_bh"] = np.where(defined, rej_bh & pos, None)
        df["informative_by"] = np.where(defined, rej_by & pos, None)
        d = df["d_r2_mean"].to_numpy(dtype=float)
        df["practically_significant"] = np.where(np.isfinite(d) & member,
                                                 d >= practical_threshold, None)
        df["underpowered"] = df["underpowered"].astype(object)
        df["in_family"] = df["in_family"].astype(bool)
    return df


def _protocol_stage(cfg, res, specs, units_all, timer):
    """Permutation control, bootstrap and the confirmatory table (+ the AD report).

    Returns (frames to write, protocol record for the run record)."""
    e = cfg["evaluation"]
    n_jobs = cfg["runtime"]["n_jobs"]
    out: dict[str, pd.DataFrame] = {}
    t0 = time.perf_counter()
    out["ad"], out["ad_summary"] = ad_report(res["predictions"], res["folds"]) \
        if e["ad"] is not None else (None, None)
    timer.t["protocol:ad"] = time.perf_counter() - t0
    p = cfg.get("protocol")
    if p is None:
        return out, None
    family = _expand_family(p, specs, e["reference"])
    nested_kwargs = dict(outer_splits=e["outer_splits"], seeds=e["seeds"],
                         inner_splits=e["inner_splits"], random_state=e["random_state"],
                         min_n=e["min_n"], min_n_test=e["min_n_test"])
    ctx = res["context"]
    learner = p["learner"]
    record = dict(family=family, family_rule=("auto: every class-3 target x (TI vs reference)"
                                              if p["family"] == "auto" else "explicit list"),
                  learner=learner, learner_families=list(_learner_families(learner, e["families"])),
                  practical_threshold=p["practical_threshold"], q=p["q"], alpha=e["alpha"],
                  multiplicity="Benjamini-Hochberg over the family's d_r2_p (BY as sensitivity)",
                  family_membership=FAMILY_RULE, min_n=e["min_n"], min_n_test=e["min_n_test"],
                  informative_rule="q <= q-level and d_r2_mean > 0",
                  practical_rule="d_r2_mean >= practical_threshold")

    # --- Y-randomised control
    perm = p["permutation"]
    perm_rows, perm_summ = [], []
    t0 = time.perf_counter()
    if perm["enabled"] and family:
        pl_learner = perm["learner"] or learner
        fams = _learner_families(pl_learner, e["families"])
        targets = list(dict.fromkeys(f["target"] for f in family))
        for t in targets:
            if t not in ctx:
                continue
            t1 = time.perf_counter()
            c = ctx[t]
            sets_t = perm["feature_sets"] or list(dict.fromkeys(
                x for f in family if f["target"] == t for x in (f["set"], f["reference"])))
            sets_t = [x for x in sets_t if x in c["sets"]]
            fs = {x: c["sets"][x] for x in sets_t}
            df = permutation_control(fs, c["y"].to_numpy(float), families=fams,
                                     n_permutations=perm["n_permutations"], seed=perm["seed"],
                                     groups=c["groups"], n_jobs=n_jobs, **nested_kwargs)
            df.insert(0, "learner", pl_learner)
            df.insert(0, "target", t)
            perm_rows.append(df)
            obs = c["cmp"][pl_learner].results
            for x in sets_t:
                sobs = _nested_stat(obs[x])["r2_pooled_mean"]
                pv_ = stats.permutation_pvalue(sobs, df.loc[df["set"] == x, "r2_pooled_mean"],
                                               quantile=perm["quantile"])
                perm_summ.append(dict(target=t, learner=pl_learner, set=x,
                                      n_permutations=perm["n_permutations"], seed=perm["seed"],
                                      observed_r2_pooled_mean=pv_["observed"], p=pv_["p"],
                                      null_mean=pv_["null_mean"], null_sd=pv_["null_sd"],
                                      null_q95=pv_["null_q95"], null_max=pv_["null_max"],
                                      n_valid=pv_["n_valid"], n_failed=pv_["n_failed"],
                                      uninformative=pv_["uninformative"], reason=pv_["reason"]))
            timer.t[f"protocol:permutation:{t}"] = time.perf_counter() - t1
            log.info("protocol permutations %s (K=%d, %d sets): %.1f s", t,
                     perm["n_permutations"], len(sets_t), timer.t[f"protocol:permutation:{t}"])
        record["permutation"] = dict(
            n_permutations=perm["n_permutations"], seed=perm["seed"], learner=pl_learner,
            families=list(fams), quantile=perm["quantile"],
            rng_rule=("permutation b of a target: permutation_indices(n, groups, "
                      "numpy.random.default_rng([seed, b])), i.e. rng.permutation(n) for a "
                      "row-wise target and the restricted block draw of `scheme` for a "
                      "block target (see scheme_by_target)"),
            statistic="mean over repeats of pooled out-of-fold R2 (NaN if any outer fold failed)",
            p_rule="(1 + #{null >= observed}) / (1 + K_valid)",
            uninformative_rule="observed <= quantile of the permutation distribution",
            folds="identical to the main analysis (rows, groups and seeds unchanged)",
            scheme=PERMUTATION_SCHEME,
            scheme_by_target={t: ("block" if ctx[t]["groups"] is not None else "row-wise")
                              for t in targets if t in ctx},
            enabled=True)
    else:
        record["permutation"] = dict(enabled=False, reason=(
            "protocol.permutation.enabled = false" if not perm["enabled"] else "empty family"))
    timer.t["protocol:permutation"] = time.perf_counter() - t0
    out["permutation"] = pd.concat(perm_rows, ignore_index=True) if perm_rows else \
        pd.DataFrame(columns=["target", "learner", "set", "permutation", "rng", "r2_pooled_mean",
                              "r2_fold_mean", "mae_fold_mean", "n_failed_folds", "reason"])
    perm_df = pd.DataFrame(perm_summ, columns=[
        "target", "learner", "set", "n_permutations", "seed", "observed_r2_pooled_mean", "p",
        "null_mean", "null_sd", "null_q95", "null_max", "n_valid", "n_failed", "uninformative",
        "reason"])
    out["permutation_summary"] = perm_df

    # --- whole-procedure paired bootstrap
    boot = p["bootstrap"]
    boot_rows, boot_summ = [], {}
    t0 = time.perf_counter()
    if boot["enabled"] and family:
        fams = tuple(boot["families"]) if boot["families"] else \
            _learner_families(learner, e["families"])
        for f in family:
            t = f["target"]
            if t not in ctx or f["set"] not in ctx[t]["sets"]:
                continue
            t1 = time.perf_counter()
            c = ctx[t]
            df, sm = paired_bootstrap(
                c["sets"], c["y"].to_numpy(float), set_name=f["set"], reference=f["reference"],
                families=fams, n_boot=boot["n_boot"], seed=boot["seed"],
                units=units_all[c["rows"]], cv_groups=c["groups"], alpha=e["alpha"],
                n_jobs=n_jobs, **nested_kwargs)
            df.insert(0, "reference", f["reference"])
            df.insert(0, "set", f["set"])
            df.insert(0, "target", t)
            boot_rows.append(df)
            boot_summ[(t, f["set"], f["reference"])] = sm
            timer.t[f"protocol:bootstrap:{t}:{f['set']}"] = time.perf_counter() - t1
            log.info("protocol bootstrap %s / %s (B=%d): %.1f s", t, f["set"], boot["n_boot"],
                     timer.t[f"protocol:bootstrap:{t}:{f['set']}"])
        record["bootstrap"] = dict(
            n_boot=boot["n_boot"], seed=boot["seed"], families=list(fams),
            families_rule=("fast subset (protocol.bootstrap.families), cost-aware sensitivity"
                           if boot["families"] else "the confirmatory learner's families"),
            rng_rule="resample b: numpy.random.default_rng([seed, b]); stats.group_bootstrap",
            enabled_rule=ENABLED_RULE_BOOT,
            unit="connectivity group (stereoisomers/isotopologues drawn together); other "
                 "molecules singly",
            cv_groups="every copy keeps its original group label (main-analysis CV group, or "
                      "the unit), so copies never straddle a fold",
            interval="percentile, alpha = evaluation.alpha; NaN if any resample failed",
            enabled=True)
    else:
        record["bootstrap"] = dict(enabled=False, enabled_rule=ENABLED_RULE_BOOT, reason=(
            "protocol.bootstrap.enabled = false" if not boot["enabled"] else "empty family"))
    timer.t["protocol:bootstrap"] = time.perf_counter() - t0
    out["bootstrap"] = pd.concat(boot_rows, ignore_index=True) if boot_rows else pd.DataFrame(
        columns=["target", "set", "reference", "resample", "rng", "n_rows", "n_units", "d_r2",
                 "d_mae", "r2_set", "r2_reference", "mae_set", "mae_reference",
                 "n_failed_folds", "reason"])

    out["confirmatory"] = confirmatory_table(
        family, res["summary"], out["ad_summary"], learner=learner,
        practical_threshold=p["practical_threshold"], q=p["q"], perm_summary=perm_df,
        boot=boot_summ, target_classes={s.name: s.target_class for s in specs},
        skipped=res["skipped"], n_values=res.get("n_values"), min_n=e["min_n"],
        permutation_note=(
            "permutation disabled (protocol.permutation.enabled = false)" if not perm["enabled"]
            else "set not in protocol.permutation.feature_sets" if perm["feature_sets"]
            else "permutation not available"),
        bootstrap_note=("bootstrap disabled (protocol.bootstrap.enabled = false)"
                        if not boot["enabled"] else "bootstrap not available"),
        ad_note="AD not fitted (evaluation.ad = null)" if e["ad"] is None else "AD not available")
    return out, record


def _json_ready(x):
    """NaN/inf -> None and numpy scalars -> Python, for the run record (NaN forbidden)."""
    if isinstance(x, Mapping):
        return {str(k): _json_ready(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_json_ready(v) for v in x]
    if isinstance(x, np.generic):
        x = x.item()
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def _check_family_targets(cfg, providers) -> None:
    """Fail fast (before the cohort stage) if the family names an unknown target."""
    p = cfg.get("protocol")
    if not p or not isinstance(p.get("family"), list):
        return
    names = {s.name for _, _, sp in providers for s in sp}
    bad = sorted({f["target"] for f in p["family"]} - names)
    if bad:
        raise ValueError(f"protocol.family names targets {bad} that no provider supplies "
                         f"({sorted(names)})")


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
    final = Path(out_dir)
    if final.exists() and (not final.is_dir() or any(final.iterdir())) and not overwrite:
        raise FileExistsError(f"{final} is not empty (pass overwrite=True to replace)")
    providers = _make_providers(cfg)
    _check_family_targets(cfg, providers)
    base = Path(base_dir).resolve() if base_dir is not None else Path.cwd().resolve()
    # Everything is written into a staging directory next to the results
    # directory and moved there only when the run has finished, so a failed run
    # leaves no partial results (and an overwritten directory stays intact).
    final.parent.mkdir(parents=True, exist_ok=True)
    out = final.parent / f".{final.name}.partial"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir()
    try:
        record, summary, mol_table_all, written = _run_stages(
            cfg, providers, out, base, timestamp_utc, client, code)
        _publish(out, final, overwrite)
    finally:
        if out.exists():
            shutil.rmtree(out, ignore_errors=True)
    return StudyResult(final, record, summary, mol_table_all, written + ["run_record.json"])


def _publish(stage: Path, final: Path, overwrite: bool) -> None:
    """Move the staged files into ``final`` (replacing its contents if ``overwrite``)."""
    if final.exists() and not final.is_dir():
        if not overwrite:
            raise FileExistsError(f"{final} exists and is not a directory")
        final.unlink()
    final.mkdir(parents=True, exist_ok=True)
    old = list(final.iterdir())
    if old and not overwrite:
        raise FileExistsError(f"{final} is not empty (pass overwrite=True to replace)")
    for p in old:
        if p.is_dir() and not p.is_symlink():
            shutil.rmtree(p)
        else:
            p.unlink()
    for p in sorted(stage.iterdir()):
        shutil.move(str(p), str(final / p.name))


def _run_stages(cfg, providers, out: Path, base: Path, timestamp_utc: str, client, code):
    timer = _Timer()

    cand, chembl, inputs, client = timer.run("cohort", _cohort_stage, cfg, out, base,
                                             timestamp_utc, client)
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
    rt = cfg["runtime"]
    for _, provider, _ in providers:
        if hasattr(provider, "set_runtime"):
            provider.set_runtime(cache_dir=rt["expdata_cache_dir"],
                                 chembl_cache_dir=rt["chembl_cache_dir"],
                                 offline=rt["offline"], client=client)
    with rdBase.BlockLogs():
        specs, tvalues, tfail, tprov = timer.run("targets", _target_stage, providers, mol_table,
                                                 mols)
    e = cfg["evaluation"]
    conn = ok["connectivity_group"].tolist()
    groups_all = _cv_groups(e["split"], mols, conn, ids)
    t0 = time.perf_counter()
    res = _evaluate_stage(cfg, wide, sets, specs, tvalues, groups_all, timer)
    timer.t["evaluate"] = time.perf_counter() - t0
    # bootstrap units: connectivity group, else the molecule itself
    units_all = np.asarray([g or i for g, i in zip(conn, ids)], dtype=object)
    t0 = time.perf_counter()
    proto, proto_record = _protocol_stage(cfg, res, specs, units_all, timer)
    timer.t["protocol"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    _write_csv(mol_table_all, out / "molecules.csv")
    feats = wide.reset_index()
    if cfg["output"]["features_format"] == "parquet":
        feats.to_parquet(out / "features.parquet", index=False, engine=_parquet_engine())
    _write_csv(feats, out / "features.csv")
    _write_csv(tvalues.reset_index(), out / "targets.csv")
    fail = pd.DataFrame(ffail + tfail, columns=["block", "row", "feature", "reason"])
    _write_csv(fail, out / "failures.csv")
    for name in ("folds", "predictions", "paired", "inner", "summary"):
        _write_csv(res[name], out / f"{name}.csv")
    for name, df in proto.items():
        if df is not None:
            _write_csv(df, out / f"{name}.csv")
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
        connectivity_groups=dict(
            groups=int(len({g for g in conn if g})),
            molecules=int(sum(1 for g in conn if g)),
            rule="standardised parents sharing the InChIKey first block (stereoisomers, "
                 "racemate/enantiomer, isotopologues): all kept, always in the same CV fold"),
        skipped_targets=dict(sorted(res["skipped"].items())),
    )
    if e["split"] == "scaffold":
        splitter = "balanced_group_kfold (Bemis-Murcko scaffolds; connectivity groups merged)"
    elif groups_all is not None:
        splitter = "balanced_group_kfold (connectivity groups; other molecules singletons)"
    else:
        splitter = "KFold(shuffle=True)"
    seeds = dict(outer_seeds=e["seeds"], inner_seed_rule="1000 * outer_seed + fold",
                 estimator_random_state=e["random_state"], outer_splitter=splitter)
    record = pv.build_run_record(
        config=cfg, timestamp_utc=timestamp_utc, chembl=chembl, seeds=seeds, inputs=inputs,
        manifest_sha256=manifest_sha, outputs=pv.file_records(out, written),
        timings_s=timer.t, counts=counts, base_dir=str(base), code=code,
        extra=dict(result_csvs=[c for c in RESULT_CSVS + AD_CSVS + PROTOCOL_CSVS
                                if c in written],
                   targets=[s.as_dict() for s in specs],
                   target_providers=_json_ready(tprov),
                   protocol=_json_ready(proto_record),
                   baselines=_baselines_record(cfg)))
    pv.write_run_record(record, out / "run_record.json")
    return record, res["summary"], mol_table_all, written


def _baselines_record(cfg: Mapping) -> dict:
    """Run-record entry: were the mandatory baselines and reference used?"""
    e = cfg["evaluation"]
    custom = bool(e.get("allow_custom_baselines", False))
    return dict(mandatory=list(MANDATORY_BASELINES), mandatory_reference=MANDATORY_REFERENCE,
                allow_custom_baselines=custom,
                missing=[b for b in MANDATORY_BASELINES if b not in cfg["feature_sets"]],
                reference=e["reference"],
                rule=("custom baselines allowed (evaluation.allow_custom_baselines = true)"
                      if custom else "enforced by resolve_config"))


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
    if cfg["cohort"]["source"] in ("chembl", "pooled"):
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
