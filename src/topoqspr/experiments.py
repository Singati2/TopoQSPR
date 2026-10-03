"""Experiment runners for EXPERIMENT_PLAN §8 (A, C, D) and the §7 size-dependence analysis.

Each runner returns an :class:`ExperimentResult` (a ``summary`` dict plus named tables)
that :meth:`ExperimentResult.write` stores as CSV + ``summary.json``. Everything is
deterministic (fixed seeds, sorted rows, no clock in any table); undefined values are NaN
with a ``reason``, never 0. The specification of every output column is in
``docs/EXPERIMENTS.md``.

* :func:`experiment_a` - cohort sensitivity: naive (v0.x retrieval) vs controlled cohort
  of one disease. Sizes, Jaccard overlaps, distributions (two-sample KS, BH), top-k
  indices by |Spearman rho| per target (BH), and optionally the nested
  ΔR²(TI − size) of both cohorts (``model=True``).
* :func:`experiment_c` - the optimism chain on one dataset: training R² →
  ordinary CV R² (v0.x protocol) → repeated nested CV R² (random folds, connectivity
  groups kept together) → scaffold-grouped nested CV R²; B_train, B_CV, B_scaffold.
* :func:`experiment_d` - the replication harness for a published single-index
  TI-QSPR study (study spec in JSON/YAML): index recomputation, in-sample fit,
  LOOCV, nested CV of the genre's own selection procedure and of the baseline ladder,
  descriptor-pool-randomised and y-permutation controls (Rücker et al. 2007),
  Q²_F3 and CCC.
* :func:`size_dependence` - cross-validated R² of every descriptor regressed on the size
  features (heavy atoms, bonds) and on size + rings, plus a per-family summary.

Quarantine (EXPERIMENT_PLAN §2, §11): until the plan is frozen, modelling on ChEMBL disease
cohorts is allowed only for the pilot (Migraine Disorders, D008881).
:func:`check_modelling_allowed` enforces this for every runner that takes MeSH IDs; the
command line (``tools/run_experiment.py``) calls it before any network access.

Only public APIs of :mod:`topoqspr.cohort`, :mod:`topoqspr.evaluate`,
:mod:`topoqspr.stats`, :mod:`topoqspr.features` and :mod:`topoqspr.models` are used,
plus :func:`topoqspr.pipeline.compute_index_table`. Standardisation, InChIKey
de-duplication and connectivity grouping follow ``pipeline`` (docs/PIPELINE.md §1) but
are re-implemented in :func:`prepare_structures` because the pipeline's stage function
is private.

References
----------
Rücker C, Rücker G, Meringer M (2007) J Chem Inf Model 47:2345-2357. doi:10.1021/ci700157b
Consonni V, Ballabio D, Todeschini R (2009) J Chem Inf Model 49:1669-1678. doi:10.1021/ci900115y
Lin LI-K (1989) Biometrics 45:255-268. doi:10.2307/2532051
Nadeau C, Bengio Y (2003) Machine Learning 52:239-281. doi:10.1023/A:1024068626366
Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300. doi:10.1111/j.2517-6161.1995.tb02031.x
Phipson B, Smyth GK (2010) Stat Appl Genet Mol Biol 9:39. doi:10.2202/1544-6115.1585
"""

from __future__ import annotations

import json
import math
import re
import unicodedata
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats as _st
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

from . import evaluate as ev
from . import features as ft
from . import stats
from .graph import standardize
from .indices import FAMILIES as INDEX_FAMILIES
from .indices import INDEX_NAMES, REGISTRY
from .models import FAST_FAMILIES, DropConstantAndDuplicate, ModelFamily
from .pipeline import compute_index_table
from .target_providers import get_target_provider

__all__ = [
    "PILOT_MESH_IDS",
    "PILOT_NAIVE_MESH_IDS",
    "QuarantineError",
    "check_modelling_allowed",
    "cohort_mesh_ids",
    "Structures",
    "prepare_structures",
    "index_table",
    "size_table",
    "cv_groups",
    "feature_set",
    "FEATURE_SET_NAMES",
    "compute_targets",
    "ExperimentResult",
    "experiment_a",
    "experiment_c",
    "experiment_d",
    "load_study_spec",
    "map_index_name",
    "SingleIndexPolynomial",
    "single_index_family",
    "POLY2_FAMILY",
    "SIZE_DEPENDENCE_FAMILIES",
    "size_dependence",
    "q2_f3",
    "ccc",
]

#: MeSH IDs on which modelling is allowed before the plan is frozen (EXPERIMENT_PLAN §2).
PILOT_MESH_IDS = frozenset({"D008881"})

#: MeSH headings the pilot's *naive* (v0.x free-text) cohort may carry: Migraine Disorders
#: and its two MeSH children that have ChEMBL drug_indication rows (ChEMBL_37 heading list:
#: D020325 Migraine with Aura, D020326 Migraine without Aura). Used only for the naive
#: cohort of Experiment A; the controlled cohort must be exactly :data:`PILOT_MESH_IDS`.
PILOT_NAIVE_MESH_IDS = frozenset({"D008881", "D020325", "D020326"})


class QuarantineError(RuntimeError):
    """Modelling was requested on a cohort that is quarantined until the plan is frozen."""


def check_modelling_allowed(mesh_ids: Iterable[str] | None, plan_frozen: bool = False, *,
                            allowed: Iterable[str] = PILOT_MESH_IDS) -> None:
    """Raise :class:`QuarantineError` unless every MeSH ID is the pilot's.

    EXPERIMENT_PLAN §2/§3.3: no model performance for any disease other than the
    quarantined pilot may be computed before the plan is frozen. ``plan_frozen=True``
    lifts the guard (pass it only once the plan's header records the frozen commit).
    ``None`` (no disease, e.g. a SMILES file) is not checked: the guard is a seat belt for
    ChEMBL disease cohorts, not a lock on arbitrary inputs. An empty or blank list is
    *not* ``None``: it means a cohort whose disease is unknown (e.g. an ATC / pooled
    approved-drug manifest with an empty ``mesh_ids`` column) and is refused.
    ``allowed`` is the permitted set (default :data:`PILOT_MESH_IDS`).
    """
    if plan_frozen or mesh_ids is None:
        return
    if isinstance(mesh_ids, str):
        mesh_ids = [mesh_ids]
    allowed = {str(a).strip().upper() for a in allowed}
    ids = sorted({str(m).strip().upper() for m in mesh_ids if str(m).strip()})
    bad = [m for m in ids if m not in allowed] if ids else ["<none: disease unknown>"]
    if bad:
        raise QuarantineError(
            f"modelling on MeSH {bad} is quarantined until EXPERIMENT_PLAN.md is frozen "
            f"(only the pilot {sorted(PILOT_MESH_IDS)} is allowed); before the freeze, "
            "non-pilot cohorts run in counts-only mode, or pass plan_frozen=True after "
            "freezing")


def cohort_mesh_ids(c) -> list[str]:
    """MeSH IDs of a cohort: ``parameters['mesh_ids']`` (if any) united with the
    ``mesh_ids`` cells (``;``-separated) of the included manifest rows. ``[]`` if none
    (disease unknown, refused by :func:`check_modelling_allowed`)."""
    ids: set[str] = set()
    params = getattr(c, "parameters", None) or {}
    pm = params.get("mesh_ids") if isinstance(params, Mapping) else None
    if pm:
        ids |= {str(m).strip().upper() for m in ([pm] if isinstance(pm, str) else pm)
                if str(m).strip()}
    man = c.manifest if hasattr(c, "manifest") else c
    if isinstance(man, pd.DataFrame) and "mesh_ids" in man.columns:
        inc = man[man["decision"] == "include"] if "decision" in man.columns else man
        for v in inc["mesh_ids"]:
            if isinstance(v, str):
                ids |= {x.strip().upper() for x in v.split(";") if x.strip()}
    return sorted(ids)


# =============================================================== structures
@dataclass
class Structures:
    """Standardised structures of a candidate list.

    ``table``: one row per candidate (``id, input_smiles, status, reason, detail, smiles,
    inchikey, parent_note, n_heavy, n_edges, connectivity_group``), as ``molecules.csv``
    of the pipeline. ``ids``/``mols``/``graphs``: the accepted molecules in input order.
    """

    table: pd.DataFrame
    ids: list[str]
    mols: list
    graphs: list

    @property
    def ok(self) -> pd.DataFrame:
        return self.table[self.table["status"] == "ok"].set_index("id").loc[self.ids]

    @property
    def inchikeys(self) -> list[str]:
        return list(self.ok["inchikey"])

    def connectivity_labels(self) -> np.ndarray:
        """Group label per accepted molecule: its connectivity group or its own id."""
        ok = self.ok
        return np.asarray([g or i for g, i in zip(ok["connectivity_group"], ok.index)],
                          dtype=object)

    @property
    def n_connectivity_groups(self) -> int:
        g = self.ok["connectivity_group"]
        return int(g[g != ""].nunique())


def prepare_structures(ids: Sequence[str], smiles: Sequence[str], *,
                       min_heavy_atoms: int = 3) -> Structures:
    """Standardise every candidate (:func:`topoqspr.graph.standardize`).

    Same rules as the pipeline's standardise stage (docs/PIPELINE.md §1): rejections are
    kept with their code; a parent whose full InChIKey equals an earlier one is rejected
    as ``duplicate_inchikey:<kept id>``; parents sharing only the InChIKey first block get
    a common ``connectivity_group`` (the first member's id).
    """
    ids = [str(i) for i in ids]
    if len(ids) != len(smiles):
        raise ValueError(f"{len(ids)} ids for {len(smiles)} SMILES")
    if len(set(ids)) != len(ids):
        dup = sorted({i for i in ids if ids.count(i) > 1})
        raise ValueError(f"duplicate ids {dup[:10]}")
    rows, kept_ids, mols, graphs = [], [], [], []
    seen: dict[str, str] = {}
    for rid, smi in zip(ids, smiles):
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
            kept_ids.append(rid)
            mols.append(res.mol)
            graphs.append(res.graph)
        rows.append(row)
    first: dict[str, str] = {}
    count: dict[str, int] = {}
    for row in rows:
        if row["status"] == "ok" and row["inchikey"]:
            b1 = row["inchikey"].split("-")[0]
            first.setdefault(b1, row["id"])
            count[b1] = count.get(b1, 0) + 1
    for row in rows:
        b1 = row["inchikey"].split("-")[0] if row["inchikey"] else ""
        row["connectivity_group"] = (first[b1] if row["status"] == "ok" and b1
                                     and count.get(b1, 0) > 1 else "")
    table = pd.DataFrame(rows, columns=["id", "input_smiles", "status", "reason", "detail",
                                        "smiles", "inchikey", "parent_note", "n_heavy",
                                        "n_edges", "connectivity_group"])
    table["n_heavy"] = table["n_heavy"].astype("Int64")
    table["n_edges"] = table["n_edges"].astype("Int64")
    return Structures(table, kept_ids, mols, graphs)


def index_table(st: Structures, names: Sequence[str] | None = None) -> tuple[pd.DataFrame, list[dict]]:
    """All (or the named) registry indices of the accepted molecules (NaN + reason)."""
    return compute_index_table(st.graphs, pd.Index(st.ids, name="id"), names)


def size_table(st: Structures) -> pd.DataFrame:
    """``n_heavy``, ``n_bonds``, ``n_rings`` of the hydrogen-suppressed graphs."""
    return ft.size_rings_features(st.mols, pd.Index(st.ids, name="id"))


FEATURE_SET_NAMES: tuple[str, ...] = (
    ("size", "size_rings", "composition", "TI", "rdkit2d")
    + tuple(f"TI_{f}" for f in INDEX_FAMILIES))


def feature_set(st: Structures, name: str, *, index_names: Sequence[str] | None = None,
                rdkit2d_exclusion: str | None = None) -> pd.DataFrame:
    """One feature set of EXPERIMENT_PLAN §5 for the accepted molecules.

    ``size``, ``size_rings``, ``composition``, ``TI`` (all or ``index_names``),
    ``TI_<family>`` (the registry indices of that family) or ``rdkit2d`` minus the
    exclusions of ``rdkit2d_exclusion`` (a key of ``features.TARGET_EXCLUSIONS``;
    required, because every target needs an explicit leakage decision).
    """
    idx = pd.Index(st.ids, name="id")
    if name == "size":
        return ft.size_features(st.mols, idx)
    if name == "size_rings":
        return ft.size_rings_features(st.mols, idx)
    if name == "composition":
        return ft.composition_features(st.mols, idx)
    if name == "TI":
        return index_table(st, index_names)[0]
    if name.startswith("TI_") and name[3:] in INDEX_FAMILIES:
        return index_table(st, [s.name for s in REGISTRY if s.family == name[3:]])[0]
    if name == "rdkit2d":
        if rdkit2d_exclusion is None:
            raise ValueError("rdkit2d needs the target's exclusion rule (rdkit2d_exclusion)")
        return ft.rdkit2d_features(st.mols, rdkit2d_exclusion, idx)
    raise ValueError(f"unknown feature set {name!r}; known: {list(FEATURE_SET_NAMES)}")


def cv_groups(st: Structures, split: str = "random") -> np.ndarray | None:
    """CV group labels as in the pipeline (docs/PIPELINE.md, *Connectivity groups*).

    ``random``: ``None`` (plain shuffled KFold) when no two molecules share an InChIKey
    first block, else one group per connectivity group and singletons elsewhere.
    ``scaffold``: Bemis–Murcko scaffold (acyclic → singleton), connectivity-group members
    taking their first member's scaffold.
    """
    conn = st.connectivity_labels()
    if split == "random":
        return None if st.n_connectivity_groups == 0 else conn
    if split == "scaffold":
        scaf = ev.scaffold_groups(st.mols)
        pos = {i: k for k, i in enumerate(st.ids)}
        return np.asarray([scaf[pos[g]] for g in conn], dtype=object)
    raise ValueError("split must be 'random' or 'scaffold'")


def compute_targets(st: Structures, provider: str = "rdkit", names: Sequence[str] | None = None,
                    options: Mapping | None = None) -> tuple[pd.DataFrame, list[dict], list]:
    """Target values for the accepted molecules from a target provider (NaN + reason).

    Returns ``(values, failures, specs)``; ``values`` is indexed by molecule id.
    """
    prov = get_target_provider(provider, names, options)
    ok = st.ok
    mt = pd.DataFrame({"input_smiles": ok["input_smiles"], "smiles": ok["smiles"],
                       "inchikey": ok["inchikey"]}, index=pd.Index(st.ids, name="id"))
    tv = prov.compute(mt, st.mols)
    vals = tv.values.astype(float)
    have = {(str(f.get("row")), str(f.get("feature"))) for f in tv.failures}
    for r, c in zip(*np.nonzero(np.isnan(vals.to_numpy()))):
        if (str(vals.index[r]), str(vals.columns[c])) not in have:
            raise ValueError(f"provider {provider!r} returned NaN without a reason at "
                             f"({vals.index[r]!r}, {vals.columns[c]!r})")
    return vals, list(tv.failures), list(prov.specs())


# =================================================================== result
def _jsonable(obj: Any) -> Any:
    if isinstance(obj, Mapping):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return f if math.isfinite(f) else None
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    return obj


def _write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False, lineterminator="\n")


@dataclass
class ExperimentResult:
    """Output of a runner: ``summary`` (JSON-able dict) and named tables."""

    name: str
    summary: dict
    tables: dict[str, pd.DataFrame] = field(default_factory=dict)

    def write(self, out_dir: str | Path) -> list[Path]:
        """Write ``<table>.csv`` for every table and ``summary.json`` (NaN -> null).

        CSVs: UTF-8, ``\\n`` line endings, no index, missing values as empty cells (read
        them with ``keep_default_na=False, na_values=[""]``).
        """
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        paths = []
        for name in sorted(self.tables):
            p = out / f"{name}.csv"
            _write_csv(self.tables[name], p)
            paths.append(p)
        p = out / "summary.json"
        p.write_text(json.dumps(_jsonable({"experiment": self.name, **self.summary}),
                                sort_keys=True, indent=1, allow_nan=False) + "\n",
                     encoding="utf-8")
        paths.append(p)
        return paths


# ============================================================ shared helpers
def _ci_of_mean(res: ev.NestedCVResult, alpha: float) -> dict:
    """Mean outer-fold R² with the corrected repeated-CV t interval (NaN if any fold NaN)."""
    f = res.folds
    t = stats.corrected_ttest(f["r2"].to_numpy(float), f["n_train"].to_numpy(float),
                              f["n_test"].to_numpy(float), alpha=alpha)
    return dict(mean=t["mean"], ci_low=t["ci_low"], ci_high=t["ci_high"], se=t["se"],
                df=t["df"])


def q2_f3(y_true, y_pred, y_train_sets: Sequence[np.ndarray]) -> float:
    """Q²_F3 (Consonni, Ballabio & Todeschini 2009) for pooled out-of-sample predictions.

        Q²_F3 = 1 − [Σ_i (y_i − ŷ_i)² / n_out] / [mean_k Σ_{j∈train_k} (y_j − ȳ_train_k)² / n_train_k]

    One external set gives the published formula; for K-fold or leave-one-out CV the
    numerator pools every out-of-fold prediction and the denominator averages the training
    variance over the K training sets. NaN if undefined.
    """
    y, p = np.asarray(y_true, float), np.asarray(y_pred, float)
    if y.size == 0 or not np.all(np.isfinite(p)):
        return math.nan
    den = [float(np.mean((t - t.mean()) ** 2)) for t in y_train_sets if len(t) > 1]
    if not den or np.mean(den) == 0:
        return math.nan
    return float(1.0 - np.mean((y - p) ** 2) / np.mean(den))


def ccc(y_true, y_pred) -> float:
    """Lin's (1989) concordance correlation coefficient (population moments); NaN if undefined."""
    y, p = np.asarray(y_true, float), np.asarray(y_pred, float)
    if y.size < 2 or not np.all(np.isfinite(p)):
        return math.nan
    num = 2.0 * np.mean((y - y.mean()) * (p - p.mean()))
    den = y.var() + p.var() + (y.mean() - p.mean()) ** 2
    return float(num / den) if den > 0 else math.nan


def _oof_secondary(res: ev.NestedCVResult) -> dict:
    """Q²_F3 and CCC per repeat of a nested CV, averaged (NaN if any repeat is NaN)."""
    pr = res.predictions
    q, c = [], []
    for _, g in pr.groupby("repeat", sort=True):
        y = g["y_true"].to_numpy(float)
        pos = g["pos"].to_numpy()
        trains = []
        for _, gk in g.groupby("fold", sort=True):
            mask = ~np.isin(pos, gk["pos"].to_numpy())
            trains.append(y[mask])
        q.append(q2_f3(y, g["y_pred"].to_numpy(float), trains))
        c.append(ccc(y, g["y_pred"].to_numpy(float)))
    return dict(q2_f3=float(np.mean(q)) if q else math.nan,
                ccc=float(np.mean(c)) if c else math.nan)


def _nested_row(res: ev.NestedCVResult, alpha: float) -> dict:
    s = res.summary()
    ci = _ci_of_mean(res, alpha)
    return dict(n=res.meta["n"], p=res.meta["p"], r2_fold_mean=s["r2_fold_mean"],
                r2_fold_ci_low=ci["ci_low"], r2_fold_ci_high=ci["ci_high"],
                r2_pooled_mean=s["r2_pooled_mean"], mae_fold_mean=s["mae_fold_mean"],
                n_failed_folds=s["n_failed_folds"], n_valid_folds=s["n_valid_folds"],
                min_n_test=res.meta["min_n_test"], underpowered=s["underpowered"],
                underpowered_reason=s["underpowered_reason"], reason=s["reason"],
                family_counts=json.dumps({str(k): v for k, v in s["family_counts"].items()},
                                          sort_keys=True))


def _r2(y: np.ndarray, p: np.ndarray) -> float:
    """R² = 1 − SSE/SST (NaN if n < 2, SST = 0 or a prediction is not finite)."""
    y, p = np.asarray(y, float).ravel(), np.asarray(p, float).ravel()
    sst = float(np.sum((y - y.mean()) ** 2)) if y.size else 0.0
    if y.size < 2 or sst == 0 or not np.all(np.isfinite(p)):
        return math.nan
    return 1.0 - float(np.sum((y - p) ** 2)) / sst


def _bh_column(p: pd.Series, alpha: float) -> np.ndarray:
    q, _ = stats.benjamini_hochberg(p.to_numpy(float), q=alpha)
    return q


# ============================================================== Experiment A
_NAIVE_KINDS = {"naive"}


def _cohort_frame(c) -> tuple[pd.DataFrame, str]:
    """(included manifest rows, kind) of a :class:`~topoqspr.cohort.Cohort` or a frame."""
    if hasattr(c, "manifest"):
        man, kind = c.manifest, getattr(c, "kind", "")
    else:
        man, kind = c, ""
    inc = man[man["decision"] == "include"].reset_index(drop=True)
    return inc, kind


def _cohort_candidates(inc: pd.DataFrame, kind: str) -> pd.DataFrame:
    """Candidate ids/SMILES: naive -> one per molecule (salts separate, as in v0.x);
    otherwise one per curated parent (as the pipeline)."""
    idcol = "molecule_chembl_id" if kind in _NAIVE_KINDS else "parent_chembl_id"
    return pd.DataFrame({"id": inc[idcol].astype(str), "smiles": inc["canonical_smiles"].astype(str),
                         "raw_inchikey": inc["standard_inchi_key"].fillna("").astype(str),
                         "parent_chembl_id": inc["parent_chembl_id"].astype(str)})


def _describe(x: np.ndarray) -> dict:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return dict(n=0, mean=math.nan, sd=math.nan, median=math.nan, q25=math.nan,
                    q75=math.nan, min=math.nan, max=math.nan)
    return dict(n=int(x.size), mean=float(x.mean()),
                sd=float(x.std(ddof=1)) if x.size > 1 else math.nan,
                median=float(np.median(x)), q25=float(np.quantile(x, 0.25)),
                q75=float(np.quantile(x, 0.75)), min=float(x.min()), max=float(x.max()))


def _correlations(ti: pd.DataFrame, y: pd.Series, alpha: float) -> pd.DataFrame:
    rows = []
    for c in ti.columns:
        r = stats.spearman_ci(ti[c].to_numpy(float), y.to_numpy(float), alpha=alpha)
        rows.append(dict(index=c, rho=r["rho"], ci_low=r["ci_low"], ci_high=r["ci_high"],
                         p=r["p"], n=r["n"], reason=r["reason"]))
    df = pd.DataFrame(rows)
    df["q"] = _bh_column(df["p"], alpha)
    return df


def _top_k(cor: pd.DataFrame, k: int, alpha: float | None = None) -> list[str]:
    """Top-k indices by |rho| (NaN excluded; ties in registry/column order).
    With ``alpha``, only indices with BH q <= alpha are eligible."""
    d = cor[np.isfinite(cor["rho"].to_numpy(float))].copy()
    if alpha is not None:
        d = d[d["q"] <= alpha]
    d["_abs"] = d["rho"].abs()
    d["_ord"] = np.arange(len(d))
    d = d.sort_values(["_abs", "_ord"], ascending=[False, True], kind="mergesort")
    return list(d["index"].head(k))


def _top_k_ties(cor: pd.DataFrame, k: int, alpha: float | None = None) -> list[str]:
    """Like :func:`_top_k`, but every index whose |rho| ties the k-th largest |rho|
    (relative tolerance 1e-12) is included, so the set does not depend on registry order;
    its size can exceed k."""
    d = cor[np.isfinite(cor["rho"].to_numpy(float))].copy()
    if alpha is not None:
        d = d[d["q"] <= alpha]
    if d.empty or k <= 0:
        return []
    a = d["rho"].abs().to_numpy(float)
    kth = np.sort(a)[::-1][min(k, a.size) - 1]
    sel = (a >= kth) | np.isclose(a, kth, rtol=1e-12, atol=0.0)
    d = d[sel].assign(_abs=a[sel], _ord=np.arange(len(d))[sel])
    d = d.sort_values(["_abs", "_ord"], ascending=[False, True], kind="mergesort")
    return list(d["index"])


def experiment_a(naive, controlled, *, target_provider: str = "rdkit",
                 target_names: Sequence[str] | None = None,
                 target_options: Mapping | None = None,
                 index_names: Sequence[str] | None = None, top_k: int = 10,
                 alpha: float = 0.05, model: bool = False,
                 families: Sequence = FAST_FAMILIES, seeds: Sequence[int] = ev.DEFAULT_SEEDS,
                 outer_splits: int = 5, inner_splits: int = 3, n_jobs: int = 1,
                 min_n: int = 50, min_n_test: int = 10, min_heavy_atoms: int = 3,
                 mesh_ids: Iterable[str] | None = None, counts_only: bool = False,
                 plan_frozen: bool = False) -> ExperimentResult:
    """Experiment A: naive vs controlled cohort of one disease.

    ``naive``/``controlled`` are :class:`~topoqspr.cohort.Cohort` objects (or manifest
    frames; a frame is treated as a parent-level cohort). Both cohorts' included
    structures are standardised with the v1.0 standardiser (so that the comparison
    isolates *which drugs* each definition retrieves; the naive salt forms collapse onto
    their parents as ``duplicate_inchikey``). With ``model=True`` the nested
    ΔR²(TI − size) of each cohort × target is computed (``compare_feature_sets``,
    identical folds, connectivity groups).

    Quarantine (EXPERIMENT_PLAN §3.3 step 2: before the freeze, candidate diseases get
    "counts only; no properties, no models"). Unless ``plan_frozen``, the target-dependent
    outputs (target distributions, correlations, top-k overlap, ``model``) need the
    controlled cohort's MeSH IDs (``mesh_ids`` united with :func:`cohort_mesh_ids`) to be
    exactly the pilot's and the naive cohort's (:func:`cohort_mesh_ids`: the headings its
    free text matched) to lie in :data:`PILOT_NAIVE_MESH_IDS`; an empty set is refused.
    Otherwise :class:`QuarantineError` is raised unless ``counts_only=True``, which
    computes only counts, overlaps, membership and the size distributions.
    """
    ctl_mesh = sorted(set(cohort_mesh_ids(controlled))
                      | {str(m).strip().upper() for m in (mesh_ids or []) if str(m).strip()})
    nav_mesh = cohort_mesh_ids(naive)
    if counts_only and model:
        raise ValueError("model=True contradicts counts_only=True")
    quarantine = ""
    if not plan_frozen:
        try:
            check_modelling_allowed(ctl_mesh, allowed=PILOT_MESH_IDS)
            try:
                check_modelling_allowed(nav_mesh, allowed=PILOT_NAIVE_MESH_IDS)
            except QuarantineError as exc:
                raise QuarantineError(f"naive cohort: {exc}") from None
        except QuarantineError as exc:
            quarantine = str(exc)
    if quarantine and not counts_only:
        raise QuarantineError(f"{quarantine}. Experiment A on this cohort pair: pass "
                              "counts_only=True (CLI --counts-only) before the freeze")
    full = not counts_only
    cohorts = {}
    for label, c, default_kind in (("naive", naive, "naive"), ("controlled", controlled, "controlled")):
        inc, kind = _cohort_frame(c)
        kind = kind or default_kind
        cand = _cohort_candidates(inc, kind)
        st = prepare_structures(cand["id"], cand["smiles"], min_heavy_atoms=min_heavy_atoms)
        cohorts[label] = dict(inc=inc, cand=cand, st=st, kind=kind)

    nv, ct = cohorts["naive"], cohorts["controlled"]
    from .cohort import jaccard

    def okkeys(c):
        return set(c["st"].ok["inchikey"]) - {""}

    j_raw = jaccard(nv["cand"]["raw_inchikey"], ct["cand"]["raw_inchikey"])
    j_parent_id = jaccard(nv["cand"]["parent_chembl_id"], ct["cand"]["parent_chembl_id"])
    j_std = jaccard(okkeys(nv), okkeys(ct))
    counts = {}
    for label, c in cohorts.items():
        t = c["st"].table
        counts[label] = dict(
            manifest_included=int(len(c["inc"])),
            distinct_parent_chembl_ids=int(c["cand"]["parent_chembl_id"].nunique()),
            standardised_ok=int((t["status"] == "ok").sum()),
            rejected={k: int(v) for k, v in sorted(
                t.loc[t["status"] != "ok", "reason"].str.split(":").str[0].value_counts().items())},
            connectivity_groups=c["st"].n_connectivity_groups)
    inter = okkeys(nv) & okkeys(ct)

    # membership table on standardised parent InChIKeys
    mem = {}
    for label, c in cohorts.items():
        ok = c["st"].ok
        for rid, key in zip(ok.index, ok["inchikey"]):
            mem.setdefault(key, {"naive": "", "controlled": ""})[label] = rid
    membership = pd.DataFrame([dict(inchikey=k, naive_id=v["naive"], controlled_id=v["controlled"],
                                    in_naive=bool(v["naive"]), in_controlled=bool(v["controlled"]))
                               for k, v in sorted(mem.items())])

    # features and targets per cohort (targets and indices only outside counts-only mode)
    tgt_names: list[str] = []
    for label, c in cohorts.items():
        st = c["st"]
        c["size"] = size_table(st)
        if full:
            c["ti"], _ = index_table(st, index_names)
            c["y"], c["y_fail"], specs = compute_targets(st, target_provider, target_names,
                                                         target_options)
            tgt_names = [s.name for s in specs]

    # distributions + KS
    variables = list(ft.SIZE_RINGS_COLUMNS) + tgt_names
    drows = []
    for v in variables:
        src = "size" if v in ft.SIZE_RINGS_COLUMNS else "y"
        xn = cohorts["naive"][src][v].to_numpy(float)
        xc = cohorts["controlled"][src][v].to_numpy(float)
        dn, dc = _describe(xn), _describe(xc)
        row = dict(variable=v, kind="size" if src == "size" else "target")
        row.update({f"naive_{k}": val for k, val in dn.items()})
        row.update({f"controlled_{k}": val for k, val in dc.items()})
        if dn["n"] >= 2 and dc["n"] >= 2:
            ks = _st.ks_2samp(xn[np.isfinite(xn)], xc[np.isfinite(xc)], method="auto")
            row.update(ks_stat=float(ks.statistic), ks_p=float(ks.pvalue), reason="")
        else:
            row.update(ks_stat=math.nan, ks_p=math.nan, reason="fewer than 2 values in a cohort")
        drows.append(row)
    dist = pd.DataFrame(drows)
    dist["ks_q"] = _bh_column(dist["ks_p"], alpha)

    # correlations and top-k overlap
    crows, orows = [], []
    for t in tgt_names:
        tops = {}
        for label, c in cohorts.items():
            cor = _correlations(c["ti"], c["y"][t], alpha)
            cor.insert(0, "target", t)
            cor.insert(0, "cohort", label)
            crows.append(cor)
            tops[label] = (_top_k(cor, top_k), _top_k(cor, top_k, alpha),
                           _top_k_ties(cor, top_k), _top_k_ties(cor, top_k, alpha))
        for variant, i in (("all", 0), ("bh_significant", 1)):
            a, b = set(tops["naive"][i]), set(tops["controlled"][i])
            at, bt = set(tops["naive"][i + 2]), set(tops["controlled"][i + 2])
            orows.append(dict(target=t, variant=variant, k=top_k,
                              naive_top=";".join(tops["naive"][i]),
                              controlled_top=";".join(tops["controlled"][i]),
                              n_naive=len(a), n_controlled=len(b), n_overlap=len(a & b),
                              jaccard=jaccard(a, b),
                              naive_top_with_ties=";".join(tops["naive"][i + 2]),
                              controlled_top_with_ties=";".join(tops["controlled"][i + 2]),
                              n_naive_with_ties=len(at), n_controlled_with_ties=len(bt),
                              n_overlap_with_ties=len(at & bt),
                              jaccard_with_ties=jaccard(at, bt)))
    correlations = pd.concat(crows, ignore_index=True) if crows else pd.DataFrame()
    overlap = pd.DataFrame(orows)

    tables = dict(membership=membership, distributions=dist)
    if full:
        tables.update(correlations=correlations, top_overlap=overlap)
    for label, c in cohorts.items():
        tables[f"molecules_{label}"] = c["st"].table

    # optional modelling step
    if model:
        mrows = []
        for label, c in cohorts.items():
            st = c["st"]
            groups_all = cv_groups(st, "random")
            sz = c["size"][list(ft.SIZE_COLUMNS)]
            for t in tgt_names:
                y = c["y"][t]
                keep = np.flatnonzero(np.isfinite(y.to_numpy(float)))
                base = dict(cohort=label, target=t, n_target_missing=int(len(y) - keep.size))
                try:
                    cmp_ = ev.compare_feature_sets(
                        {"size": sz.iloc[keep], "TI": c["ti"].iloc[keep]}, y.iloc[keep],
                        reference="size", families=families, alpha=alpha, seeds=seeds,
                        outer_splits=outer_splits, inner_splits=inner_splits,
                        groups=None if groups_all is None else groups_all[keep],
                        n_jobs=n_jobs, min_n=min_n, min_n_test=min_n_test)
                except ValueError as exc:
                    mrows.append(dict(base, reason=f"not evaluated: {exc}"[:300]))
                    continue
                s = cmp_.summary.set_index("set")
                mrows.append(dict(base, n=int(keep.size),
                                  r2_size=s.loc["size", "r2_fold_mean"],
                                  r2_TI=s.loc["TI", "r2_fold_mean"],
                                  d_r2_mean=s.loc["TI", "d_r2_mean"],
                                  d_r2_ci_low=s.loc["TI", "d_r2_ci_low"],
                                  d_r2_ci_high=s.loc["TI", "d_r2_ci_high"],
                                  d_r2_p=s.loc["TI", "d_r2_p"],
                                  underpowered=bool(s.loc["TI", "underpowered"]),
                                  reason="; ".join(x for x in (s.loc["size", "reason"],
                                                               s.loc["TI", "reason"]) if x)))
        tables["model"] = pd.DataFrame(mrows)

    summary = dict(
        counts=counts,
        jaccard=dict(parent_chembl_id=j_parent_id, structures_used=j_raw,
                     standardised_parent_inchikey=j_std),
        intersection_standardised=len(inter),
        controlled_in_naive_share=(len(inter) / len(okkeys(ct)) if okkeys(ct) else math.nan),
        mode="counts_only" if counts_only else "full", quarantine=quarantine,
        mesh_ids=dict(controlled=ctl_mesh, naive=nav_mesh),
        targets=tgt_names, top_k=top_k, alpha=alpha, model=model,
        ks_note="KS p-values compare overlapping, non-independent samples; descriptive only",
        settings=dict(index_names=list(index_names) if index_names else "all",
                      families=[getattr(f, "name", f) for f in families], seeds=list(seeds),
                      outer_splits=outer_splits, inner_splits=inner_splits) if model else {},
    )
    return ExperimentResult("A", summary, tables)


# ============================================================== Experiment C
def _group_stats(groups: np.ndarray) -> dict:
    _, counts = np.unique(np.asarray(groups, dtype=object).astype(str), return_counts=True)
    return dict(n_groups=int(counts.size), largest_group=int(counts.max()),
                largest_group_share=float(counts.max() / counts.sum()),
                n_singletons=int(np.sum(counts == 1)))


def experiment_c(X: pd.DataFrame, y, *, connectivity: Sequence[str] | None = None,
                 scaffolds: Sequence[str] | None = None, mols: Sequence | None = None,
                 families: Sequence = FAST_FAMILIES, legacy_families: Mapping | None = None,
                 seeds: Sequence[int] = ev.DEFAULT_SEEDS, outer_splits: int = 5,
                 inner_splits: int = 3, alpha: float = 0.05, ad: str | None = "knn",
                 n_jobs: int = 1, min_n: int = 50, min_n_test: int = 10,
                 label: str = "") -> ExperimentResult:
    """Experiment C: the optimism chain on one dataset (EXPERIMENT_PLAN §8 C).

    1. **training R²**: :func:`~topoqspr.evaluate.legacy_protocol`'s best family refitted
       on all rows and scored on them (``r2_train``); also the v1.0 inner-CV selection
       among ``families`` refitted on all rows (``r2_train_selected``, secondary);
    2. **ordinary CV R²**: the v0.x protocol's best-family fold-mean R² (``r2_cv``);
    3. **nested CV R²** (random outer folds; members of a connectivity group, i.e. the
       same InChIKey first block, kept in one fold): ``r2_nested``;
    4. **scaffold-grouped nested CV R²** (Bemis–Murcko groups from ``scaffolds`` or
       ``mols``; connectivity members follow their first member): ``r2_scaffold``.

    B_train = r2_train − r2_nested, B_CV = r2_cv − r2_nested,
    B_scaffold = r2_nested − r2_scaffold (fold-mean R²; pooled variants too).
    CIs: the nested and scaffold means get the corrected repeated-CV t interval
    (Nadeau & Bengio 2003). B_train and B_CV treat r2_train / r2_cv as fixed numbers
    (their CI is the nested interval reflected). B_scaffold's CI combines the two corrected
    standard errors as if the two CV estimates were independent (approximate; they share
    the data, so it is likely conservative). Rows with a missing ``y`` are dropped and
    counted; NaN features are median-imputed inside nested CV but dropped by the v0.x
    protocol (as v0.x did; reported as ``legacy_n_used``). B_train and B_CV always compare
    R² values on the same molecules: when the v0.x protocol drops rows, the nested CV is
    re-run on its complete-case rows (step ``nested_random_legacy_rows``) and that is the
    reference (``B_rows`` says which rows; ``B_train_all_rows``/``B_CV_all_rows`` keep the
    all-rows contrast for comparison). If the v0.x protocol cannot run (< 10 complete-case
    rows) its values are NaN with ``legacy_reason``.
    """
    if not isinstance(X, pd.DataFrame):
        X = pd.DataFrame(np.asarray(X, dtype=float))
    y = pd.Series(np.asarray(y, dtype=float).ravel(), index=X.index)
    if len(y) != len(X):
        raise ValueError("X and y lengths differ")
    keep = np.flatnonzero(np.isfinite(y.to_numpy()))
    n_missing = int(len(y) - keep.size)
    Xk, yk = X.iloc[keep], y.iloc[keep]
    conn = None
    if connectivity is not None:
        conn = np.asarray(connectivity, dtype=object)[keep]
    if scaffolds is None and mols is not None:
        scaffolds = ev.scaffold_groups(list(mols))
    scaf = None
    if scaffolds is not None:
        scaf = np.asarray(scaffolds, dtype=object)[keep]
        if conn is not None:  # connectivity members follow their first member's scaffold
            first = {}
            for c, s in zip(conn, scaf):
                first.setdefault(c, s)
            scaf = np.asarray([first[c] for c in conn], dtype=object)
    random_groups = None
    if conn is not None and len(set(conn)) < len(conn):
        random_groups = conn
    common = dict(seeds=seeds, outer_splits=outer_splits, inner_splits=inner_splits, ad=ad,
                  n_jobs=n_jobs, min_n=min_n, min_n_test=min_n_test)

    try:
        legacy = ev.legacy_protocol(Xk, yk, families=legacy_families)
        legacy_reason = ""
    except ValueError as exc:
        legacy_reason = f"v0.x protocol not evaluated: {exc}"
        cc_rows = np.all(np.isfinite(Xk.to_numpy(float)), axis=1)
        legacy = ev.LegacyResult(pd.DataFrame(columns=["family", "R2_mean", "R2_std",
                                                       "MAE_mean", "params", "reason"]),
                                 None, math.nan, math.nan, int(cc_rows.sum()),
                                 list(Xk.index[~cc_rows]))
    sel = ev.select_and_refit(Xk.to_numpy(float), yk.to_numpy(float), families, inner_splits,
                              random_state=0, groups=random_groups)
    r2_train_sel = math.nan
    if sel.estimator is not None:
        r2_train_sel = _r2(yk.to_numpy(float),
                           np.asarray(sel.estimator.predict(Xk.to_numpy(float)), dtype=float))
    nested = ev.nested_cv(Xk, yk, families, groups=random_groups, **common)
    rows = [dict(step="nested_random", **_nested_row(nested, alpha))]
    ci_n = _ci_of_mean(nested, alpha)
    sn = nested.summary()
    scaffold_reason, scaffold_res, sc_stats = "", None, {}
    if scaf is None:
        scaffold_reason = "no scaffolds or molecules given"
    else:
        sc_stats = _group_stats(scaf)
        try:
            scaffold_res = ev.nested_cv(Xk, yk, families, groups=scaf, **common)
        except ValueError as exc:
            scaffold_reason = f"scaffold CV not possible: {exc}"
    if scaffold_res is not None:
        rows.append(dict(step="nested_scaffold", **_nested_row(scaffold_res, alpha)))
        ci_s = _ci_of_mean(scaffold_res, alpha)
        ss = scaffold_res.summary()
    else:
        rows.append(dict(step="nested_scaffold", reason=scaffold_reason))
        ci_s = dict(mean=math.nan, ci_low=math.nan, ci_high=math.nan, se=math.nan, df=math.nan)
        ss = dict(r2_fold_mean=math.nan, r2_pooled_mean=math.nan, reason=scaffold_reason)

    r2n, r2s = sn["r2_fold_mean"], ss["r2_fold_mean"]
    b_sc = r2n - r2s
    # B_train / B_CV must compare R² values on the same molecules: the v0.x protocol drops
    # rows with any NaN feature, so the nested reference is re-run on those rows if needed.
    cc = np.all(np.isfinite(Xk.to_numpy(float)), axis=1)
    cc_reason, nested_cc = "", nested
    if legacy_reason:
        nested_cc, cc_reason = None, "v0.x protocol not evaluated, no same-rows reference"
    elif not cc.all():
        sub = None
        if random_groups is not None:
            g = random_groups[cc]
            sub = g if len(set(g)) < len(g) else None
        try:
            nested_cc = ev.nested_cv(Xk[cc], yk[cc], families, groups=sub, **common)
        except ValueError as exc:
            nested_cc = None
            cc_reason = f"nested CV on the v0.x rows not possible: {exc}"[:300]
        rows.append(dict(step="nested_random_legacy_rows", reason=cc_reason) if nested_cc is None
                    else dict(step="nested_random_legacy_rows", **_nested_row(nested_cc, alpha)))
    if nested_cc is not None:
        scc, ci_c = nested_cc.summary(), _ci_of_mean(nested_cc, alpha)
    else:
        scc = dict(r2_fold_mean=math.nan, r2_pooled_mean=math.nan, reason=cc_reason)
        ci_c = dict(mean=math.nan, ci_low=math.nan, ci_high=math.nan, se=math.nan, df=math.nan)
    r2c = scc["r2_fold_mean"]
    b_train, b_cv = legacy.train_r2 - r2c, legacy.cv_r2 - r2c
    b_rows = ("all rows" if cc.all() else
              f"v0.x complete-case rows ({int(cc.sum())} of {len(cc)}; NaN-feature rows dropped)")
    if math.isfinite(ci_s["se"]) and math.isfinite(ci_n["se"]):
        se = math.sqrt(ci_n["se"] ** 2 + ci_s["se"] ** 2)
        df = min(ci_n["df"], ci_s["df"])
        crit = float(_st.t.ppf(1 - alpha / 2, df))
        b_sc_ci = (b_sc - crit * se, b_sc + crit * se)
    else:
        b_sc_ci = (math.nan, math.nan)
    chain_rows = [
        dict(step="train", r2=legacy.train_r2, ci_low=math.nan, ci_high=math.nan,
             detail=legacy_reason or f"v0.x best family {legacy.best_family} refitted on all "
             f"{legacy.n_used} complete-case rows"),
        dict(step="train_selected", r2=r2_train_sel, ci_low=math.nan, ci_high=math.nan,
             detail=f"v1.0 inner-CV selection ({sel.family}) refitted on all rows"),
        dict(step="cv_legacy", r2=legacy.cv_r2, ci_low=math.nan, ci_high=math.nan,
             detail=legacy_reason or "v0.x protocol: tuned on all rows, best-of-families "
             "fold-mean R2"),
        dict(step="nested_random", r2=r2n, ci_low=ci_n["ci_low"], ci_high=ci_n["ci_high"],
             detail=sn["reason"]),
    ]
    if not cc.all():
        chain_rows.append(dict(step="nested_random_legacy_rows", r2=r2c, ci_low=ci_c["ci_low"],
                               ci_high=ci_c["ci_high"],
                               detail=scc["reason"] or f"same {int(cc.sum())} rows as the v0.x "
                               "protocol (reference of B_train and B_CV)"))
    chain_rows.append(dict(step="nested_scaffold", r2=r2s, ci_low=ci_s["ci_low"],
                           ci_high=ci_s["ci_high"], detail=ss["reason"]))
    chain = pd.DataFrame(chain_rows)
    optimism = pd.DataFrame([
        dict(term="B_train", value=b_train, ci_low=legacy.train_r2 - ci_c["ci_high"],
             ci_high=legacy.train_r2 - ci_c["ci_low"], pooled=legacy.train_r2 - scc["r2_pooled_mean"],
             ci_method="r2_train fixed; nested corrected-t CI reflected", rows=b_rows),
        dict(term="B_CV", value=b_cv, ci_low=legacy.cv_r2 - ci_c["ci_high"],
             ci_high=legacy.cv_r2 - ci_c["ci_low"], pooled=legacy.cv_r2 - scc["r2_pooled_mean"],
             ci_method="r2_cv fixed; nested corrected-t CI reflected", rows=b_rows),
        dict(term="B_scaffold", value=b_sc, ci_low=b_sc_ci[0], ci_high=b_sc_ci[1],
             pooled=sn["r2_pooled_mean"] - ss["r2_pooled_mean"],
             ci_method="sqrt(se_random^2 + se_scaffold^2), df = min; approximate",
             rows="all rows"),
    ])
    legacy_tab = legacy.table.copy()
    summary = dict(
        label=label, n=int(len(yk)), n_target_missing=n_missing, p=int(X.shape[1]),
        r2_train=legacy.train_r2, r2_train_selected=r2_train_sel, r2_cv=legacy.cv_r2,
        r2_nested=r2n, r2_nested_pooled=sn["r2_pooled_mean"], r2_scaffold=r2s,
        r2_scaffold_pooled=ss["r2_pooled_mean"], B_train=b_train, B_CV=b_cv, B_scaffold=b_sc,
        B_rows=b_rows, r2_nested_legacy_rows=r2c,
        B_train_all_rows=legacy.train_r2 - r2n, B_CV_all_rows=legacy.cv_r2 - r2n,
        legacy_reason=legacy_reason, legacy_rows_reason=cc_reason,
        legacy_best_family=legacy.best_family, legacy_n_used=legacy.n_used,
        legacy_dropped=[str(d) for d in legacy.dropped],
        connectivity_groups=0 if random_groups is None else int(
            np.sum(np.unique(random_groups.astype(str), return_counts=True)[1] > 1)),
        random_split="balanced_group_kfold (connectivity groups)" if random_groups is not None
        else "KFold(shuffle=True)",
        scaffold=sc_stats, scaffold_reason=scaffold_reason,
        underpowered=bool(nested.meta["underpowered"]),
        settings=dict(families=[getattr(f, "name", f) for f in families],
                      legacy_families=sorted(legacy_families) if legacy_families else "all v0.x",
                      seeds=list(seeds), outer_splits=outer_splits, inner_splits=inner_splits,
                      ad=ad, alpha=alpha),
    )
    tables = dict(chain=chain, optimism=optimism, nested=pd.DataFrame(rows), legacy=legacy_tab,
                  folds_random=nested.folds)
    if scaffold_res is not None:
        tables["folds_scaffold"] = scaffold_res.folds
    return ExperimentResult("C", summary, tables)


# ============================================================== Experiment D
_FORMS = {"linear": 1, "quadratic": 2, "cubic": 3}


def _norm(s: str) -> str:
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode()
    s = s.lower().replace("index", "").replace("indices", "")
    return re.sub(r"[^a-z0-9]", "", s)


#: Literature names -> registry names (see docs/INDEX_VERIFICATION.md).
INDEX_ALIASES: dict[str, str] = {
    "first zagreb": "M1", "second zagreb": "M2", "forgotten": "F", "hyper zagreb": "HM",
    "randic": "R", "connectivity": "R", "product connectivity": "R",
    "reciprocal randic": "RR", "sum connectivity": "SC", "harmonic": "H",
    "geometric arithmetic": "GA", "atom bond connectivity": "ABC",
    "augmented zagreb": "AZI", "sombor": "SO", "inverse sum indeg": "ISI",
    "redefined second zagreb": "ISI", "ReZ2": "ISI", "symmetric division deg": "SDD",
    "symmetric division": "SDD", "sigma": "sigma", "atom bond sum connectivity": "ABS",
    "redefined third zagreb": "ReZ3", "arithmetic geometric": "AG",
    "first gourava": None, "wiener": "W", "balaban": "J", "balaban j": "J", "szeged": "Sz",
    "edge szeged": "Sz_e", "padmakar ivan": "PI", "pi": "PI", "vertex pi": "PI_v",
    "mostar": "Mo", "edge mostar": "Mo_e", "edge wiener": "W_e", "degree distance": "DD",
    "gutman": "Gut", "hosoya": "Z", "hosoya z": "Z", "degree entropy": "I_deg",
    "M1": "M1", "M2": "M2", "HM": "HM", "AZI": "AZI", "HZ": "HM",
}


def _literature_names() -> dict[str, str]:
    out = {}
    for s in REGISTRY:
        m = re.match(r"literature name (\S+)", s.description or "")
        if m:
            out[m.group(1)] = s.name
    return out


def map_index_name(label: str) -> tuple[str | None, str]:
    """Map a study's index label to a registry name: ``(name or None, how)``.

    Order: exact registry name; literature name recorded in the registry (``ABC_4``,
    ``GA_5``, ``ND_1``...); :data:`INDEX_ALIASES` (normalised: case, spaces, accents and
    the word "index" ignored); a unique case-insensitive registry match. Anything else is
    unmappable (``None`` with the reason) and is reported, never guessed.
    """
    lab = str(label).strip()
    if lab in INDEX_NAMES:
        return lab, "registry name"
    lit = _literature_names()
    if lab in lit:
        return lit[lab], "literature name in registry"
    aliases = {_norm(k): v for k, v in INDEX_ALIASES.items()}
    n = _norm(lab)
    if n in aliases:
        if aliases[n] is None:
            return None, "known index not in the v1.0 registry"
        return aliases[n], "alias"
    ci = [s for s in INDEX_NAMES if _norm(s) == n]
    if len(ci) == 1:
        return ci[0], "case-insensitive registry name"
    if len(ci) > 1:
        return None, f"ambiguous: matches {ci}"
    return None, "no registry index with this name"


def load_study_spec(path: str | Path) -> dict:
    """Read a study spec (JSON, or YAML if PyYAML is installed). See docs/EXPERIMENTS.md."""
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() in (".yaml", ".yml"):
        try:
            import yaml  # type: ignore
        except ImportError as exc:  # pragma: no cover - depends on the environment
            raise ValueError("YAML study specs need PyYAML; use JSON") from exc
        spec = yaml.safe_load(text)
    else:
        spec = json.loads(text)
    if not isinstance(spec, dict):
        raise ValueError(f"{path}: the study spec must be a mapping")
    return spec


_SPEC_KEYS = {"name", "reference", "notes", "compounds", "properties", "indices",
              "model_forms", "reported_fits", "reported_index_values"}


def _name(v: Any) -> str:
    """A compound name as used everywhere in Experiment D: ``str(v).strip()`` (papers often
    number their compounds, so JSON ints are accepted); ``None``/bool/blank are refused."""
    if v is None or isinstance(v, bool) or not str(v).strip():
        raise ValueError(f"every compound needs a non-blank 'name' (got {v!r})")
    return str(v).strip()


def _form(v: Any) -> str:
    f = str(v).strip().lower()
    if f not in _FORMS:
        raise ValueError(f"unknown model form {v!r}; allowed: {sorted(_FORMS)}")
    return f


def _check_spec(spec: Mapping) -> dict:
    """Validate a study spec and return a normalised copy: compound names (and the keys
    of ``reported_index_values``) as ``str(name).strip()``, properties as strings, model
    forms lower-case, every reported fit with a known form (one whose property is not
    analysed is listed as not compared, not refused)."""
    unknown = sorted(set(spec) - _SPEC_KEYS)
    if unknown:
        raise ValueError(f"unknown study-spec keys {unknown}; allowed: {sorted(_SPEC_KEYS)}")
    for k in ("compounds", "properties", "indices"):
        if not spec.get(k):
            raise ValueError(f"study spec needs a non-empty {k!r}")
    out = dict(spec)
    out["model_forms"] = list(dict.fromkeys(_form(f) for f in spec.get("model_forms", ["linear"])))
    out["properties"] = [str(p).strip() for p in spec["properties"]]
    comps = []
    for c in spec["compounds"]:
        if not isinstance(c, Mapping):
            raise ValueError(f"each compound must be a mapping, got {c!r}")
        comps.append(dict(c, name=_name(c.get("name"))))
    names = [c["name"] for c in comps]
    if len(set(names)) != len(names):
        dup = sorted({n for n in names if names.count(n) > 1})
        raise ValueError(f"every compound needs a unique 'name' (after str/strip): {dup}")
    out["compounds"] = comps
    riv = spec.get("reported_index_values") or {}
    out["reported_index_values"] = {_name(k): v for k, v in riv.items()}
    fits = []
    for r in spec.get("reported_fits") or []:
        if not isinstance(r, Mapping) or "property" not in r or "index" not in r:
            raise ValueError(f"a reported fit needs 'property' and 'index': {r!r}")
        prop = str(r["property"]).strip()
        fits.append(dict(r, property=prop, form=_form(r.get("form", "linear"))))
    out["reported_fits"] = fits
    return out


class SingleIndexPolynomial(RegressorMixin, BaseEstimator):
    """y = Σ_{j≤degree} b_j x^j on one column of X (the genre's model form).

    ``column`` selects the descriptor; NaN in that column is replaced by the training
    median of the column (fitted on training rows only); fitting needs at least
    ``degree + 2`` rows with ``degree + 1`` distinct values (else ``ValueError``, so the
    grid point is not selectable). x is centred and scaled with training moments before
    the least-squares fit (numerical stability; the fitted polynomial is the same).
    """

    def __init__(self, column: int = 0, degree: int = 1):
        self.column = column
        self.degree = degree

    def fit(self, X, y):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).ravel()
        x = X[:, self.column]
        fin = np.isfinite(x)
        if not fin.any():
            raise ValueError("column is all NaN in the training rows")
        self.median_ = float(np.median(x[fin]))
        x = np.where(fin, x, self.median_)
        if len(y) < self.degree + 2 or np.unique(x).size < self.degree + 1:
            raise ValueError(f"degree {self.degree} needs >= {self.degree + 2} rows and "
                             f">= {self.degree + 1} distinct values")
        self.mu_, self.sd_ = float(x.mean()), float(x.std())
        if self.sd_ == 0:
            raise ValueError("constant column")
        V = np.vander((x - self.mu_) / self.sd_, self.degree + 1, increasing=True)
        self.coef_, *_ = np.linalg.lstsq(V, y, rcond=None)
        self.n_features_in_ = X.shape[1]
        return self

    def predict(self, X):
        x = np.asarray(X, dtype=float)[:, self.column]
        x = np.where(np.isfinite(x), x, self.median_)
        V = np.vander((x - self.mu_) / self.sd_, self.degree + 1, increasing=True)
        return V @ self.coef_


def single_index_family(p: int, degrees: Sequence[int] = (1,), name: str = "single_index") -> ModelFamily:
    """Family choosing one column and one polynomial degree by inner CV (the genre's
    "best index, best form" selection, nested)."""
    degrees = sorted({int(d) for d in degrees})
    return ModelFamily(name, lambda rs: SingleIndexPolynomial(),
                       {"column": list(range(max(p, 1))), "degree": degrees},
                       constraint=lambda g, n, pp: g["column"] < pp and n >= g["degree"] + 2)


def _poly_fit(x: np.ndarray, y: np.ndarray, degree: int) -> dict:
    """In-sample least-squares polynomial fit: R², adjusted R², coefficients (ascending
    powers of the original x). NaN + reason if undefined."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = int(x.size)
    out = dict(n=n, r2=math.nan, adj_r2=math.nan, coefficients="", reason="")
    if n < degree + 2:
        out["reason"] = f"n={n} < degree+2 (fit saturated or undefined)"
        return out
    if np.unique(x).size < degree + 1:
        out["reason"] = f"fewer than {degree + 1} distinct index values"
        return out
    sst = float(np.sum((y - y.mean()) ** 2))
    if sst == 0:
        out["reason"] = "constant property"
        return out
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        poly = np.polynomial.Polynomial.fit(x, y, degree)
    sse = float(np.sum((y - poly(x)) ** 2))
    r2 = 1.0 - sse / sst
    out.update(r2=r2, adj_r2=1.0 - (1.0 - r2) * (n - 1) / (n - degree - 1),
               coefficients=json.dumps([float(c) for c in poly.convert().coef]))
    return out


def _loo(x: np.ndarray, y: np.ndarray, degree: int) -> dict:
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = int(x.size)
    out = dict(loo_q2=math.nan, loo_rmse=math.nan, loo_q2_f3=math.nan, loo_ccc=math.nan,
               loo_reason="")
    if n < degree + 3:
        out["loo_reason"] = f"n={n} < degree+3"
        return out
    pred = np.full(n, np.nan)
    trains = []
    for i in range(n):
        m = np.arange(n) != i
        try:
            est = SingleIndexPolynomial(0, degree).fit(x[m, None], y[m])
        except ValueError as exc:
            out["loo_reason"] = f"fold {i}: {exc}"
            return out
        pred[i] = est.predict(x[i:i + 1, None])[0]
        trains.append(y[m])
    sst = float(np.sum((y - y.mean()) ** 2))
    press = float(np.sum((y - pred) ** 2))
    out.update(loo_q2=1.0 - press / sst if sst > 0 else math.nan,
               loo_rmse=math.sqrt(press / n), loo_q2_f3=q2_f3(y, pred, trains),
               loo_ccc=ccc(y, pred))
    return out


def _r2_projection(Q: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """R² of least squares on the column space of Q (orthonormal, includes the
    intercept) for every column of Y."""
    fitted = Q @ (Q.T @ Y)
    sse = np.sum((Y - fitted) ** 2, axis=0)
    sst = np.sum((Y - Y.mean(axis=0)) ** 2, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(sst > 0, 1.0 - sse / sst, np.nan)


def _basis(x: np.ndarray, degree: int) -> np.ndarray | None:
    if np.unique(x).size < degree + 1 or x.size < degree + 2:
        return None
    z = (x - x.mean()) / x.std()
    Q, _ = np.linalg.qr(np.vander(z, degree + 1, increasing=True))
    return Q


def _best_r2_arg(Xc: np.ndarray, Y: np.ndarray, degrees: Sequence[int]
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """max over columns of Xc and degrees of in-sample polynomial R², per column of Y.

    Each column is fitted on the rows where it is finite (as :func:`_poly_fit` does), so
    an index that is undefined for some compounds is not dropped. Returns
    ``(best, column, degree)``; ``best`` is NaN and column/degree ``-1`` where no fit is
    defined. Ties keep the first (column, degree) in input order.
    """
    best = np.full(Y.shape[1], -np.inf)
    arg_j = np.full(Y.shape[1], -1)
    arg_d = np.full(Y.shape[1], -1)
    for j in range(Xc.shape[1]):
        m = np.isfinite(Xc[:, j])
        x, Ym = Xc[m, j], Y[m]
        for d in degrees:
            Q = _basis(x, d)
            if Q is None:
                continue
            r = _r2_projection(Q, Ym)
            upd = np.isfinite(r) & (r > best)
            best[upd], arg_j[upd], arg_d[upd] = r[upd], j, d
    return np.where(np.isfinite(best), best, np.nan), arg_j, arg_d


def _best_r2(Xc: np.ndarray, Y: np.ndarray, degrees: Sequence[int]) -> np.ndarray:
    """:func:`_best_r2_arg` without the arg-max."""
    return _best_r2_arg(Xc, Y, degrees)[0]


def _perm_p(null: np.ndarray, obs: float) -> float:
    v = null[np.isfinite(null)]
    if not math.isfinite(obs) or v.size == 0:
        return math.nan
    return (1 + int(np.sum(v >= obs))) / (1 + v.size)


def experiment_d(spec: Mapping, *, families: Sequence = ("ridge", "pls", "linear"),
                 seeds: Sequence[int] = ev.DEFAULT_SEEDS, outer_splits: int = 5,
                 inner_splits: int = 3, alpha: float = 0.05, n_random: int = 1000,
                 n_permutations: int = 1000, random_seed: int = 0, n_jobs: int = 1,
                 min_heavy_atoms: int = 3, scaffold: bool = True, min_n: int = 50,
                 min_n_test: int = 10) -> ExperimentResult:
    """Experiment D: replication harness for one published single-index TI-QSPR study.

    Steps per property (docs/EXPERIMENTS.md §D):
    (i) recompute the study's indices (labels mapped to the registry; unmappable labels
    and compounds without a usable SMILES are reported); compare with the study's printed
    index values if given; (ii) in-sample polynomial fit of the property on each index
    for each model form, compared with the reported R²; (iii) LOOCV of each fit;
    (iv) nested CV of the genre's procedure (choose index and form by inner CV) on the
    study's indices and, as baselines, on size and composition, plus the ``families``
    ladder (null, size, composition, study TIs) with ΔR² against size; scaffold-grouped
    nested CV of the genre procedure; (v) controls for the best in-sample R²:
    ``n_random`` draws of as many i.i.d. N(0,1) pseudo-descriptors, ``n_random`` draws of
    as many registry indices, and ``n_permutations`` y-permutations (Rücker et al. 2007);
    (vi) Q²_F3 and CCC of LOOCV and nested predictions.

    ``min_n``/``min_n_test`` only set the ``underpowered`` flag of the nested rows
    (``evaluate.nested_cv``); the defaults (50, 10) are the main-panel thresholds, so a
    genre-sized study is flagged at them by design. Compound names are used as
    ``str(name).strip()`` everywhere (numbered compounds are fine).
    """
    spec = _check_spec(spec)
    forms = list(spec["model_forms"])
    degrees = [_FORMS[f] for f in forms]
    form_of = {d: f for f, d in _FORMS.items()}
    props = [str(p) for p in spec["properties"]]

    # compounds
    crow, ids, smis = [], [], []
    for c in spec["compounds"]:
        smi = (c.get("smiles") or "").strip()
        if not smi:
            crow.append(dict(name=c["name"], status="unresolved",
                             reason="no SMILES in the spec (name resolution is not automated; "
                                    "record a SMILES and its source)"))
            continue
        ids.append(str(c["name"]))
        smis.append(smi)
    st = prepare_structures(ids, smis, min_heavy_atoms=min_heavy_atoms)
    stt = st.table.set_index("id")
    for c in spec["compounds"]:
        if c["name"] in stt.index:
            r = stt.loc[c["name"]]
            crow.append(dict(name=c["name"], status=r["status"], reason=r["reason"],
                             input_smiles=r["input_smiles"], smiles=r["smiles"],
                             inchikey=r["inchikey"], connectivity_group=r["connectivity_group"]))
    order = {c["name"]: i for i, c in enumerate(spec["compounds"])}
    compounds = pd.DataFrame(crow).sort_values("name", key=lambda s: s.map(order),
                                                kind="mergesort").reset_index(drop=True)

    # index mapping
    labels = [i if isinstance(i, str) else i.get("label") or i.get("name") for i in spec["indices"]]
    mrows, mapped = [], {}
    for lab in labels:
        name, how = map_index_name(lab)
        mrows.append(dict(study_label=lab, registry_name=name or "", status="mapped" if name else
                          "unmappable", how=how))
        if name and name not in mapped.values():
            mapped[lab] = name
    index_map = pd.DataFrame(mrows)
    reg_names = list(dict.fromkeys(mapped.values()))
    all_ti, ti_fail = index_table(st)
    ti = all_ti[reg_names] if reg_names else pd.DataFrame(index=all_ti.index)
    values = ti.copy()
    values.insert(0, "compound", values.index)
    values = values.reset_index(drop=True)

    # reported index values
    chk = []
    for comp, vals in (spec.get("reported_index_values") or {}).items():
        for lab, v in vals.items():
            name = mapped.get(lab) or map_index_name(lab)[0]
            got = math.nan
            if name and comp in all_ti.index and name in all_ti.columns:
                got = float(all_ti.loc[comp, name])
            rep = float(v)
            diff = got - rep
            rel = abs(diff) / abs(rep) if rep != 0 else abs(diff)
            chk.append(dict(compound=comp, study_label=lab, registry_name=name or "",
                            reported=rep, recomputed=got, diff=diff, rel_diff=rel,
                            agree=bool(math.isfinite(got) and (abs(diff) <= 0.01 or rel <= 1e-3)),
                            reason="" if math.isfinite(got) else "not recomputed (unmapped, "
                            "unresolved compound or undefined index)"))
    index_check = pd.DataFrame(chk, columns=["compound", "study_label", "registry_name",
                                             "reported", "recomputed", "diff", "rel_diff",
                                             "agree", "reason"])

    # property values (NaN + reason if missing)
    yrows, yfail = {}, []
    for c in spec["compounds"]:
        if c["name"] not in st.ids:
            continue
        pv_ = c.get("properties") or {}
        row = {}
        for p in props:
            v = pv_.get(p)
            try:
                row[p] = float(v) if v is not None and v != "" else math.nan
            except (TypeError, ValueError):
                row[p] = math.nan
            if not math.isfinite(row[p]):
                row[p] = math.nan
                yfail.append(dict(row=c["name"], feature=p, reason=f"missing or non-numeric: {v!r}"))
        yrows[c["name"]] = row
    Y = pd.DataFrame.from_dict(yrows, orient="index").reindex(st.ids)[props]

    # reported fits lookup (every reported fit is listed in reported_fits.csv; one that
    # matches no fitted (property, index, form) says why)
    rep_fit, rep_rows = {}, []
    for r in spec["reported_fits"]:
        lab = r["index"]
        nm = mapped.get(lab) or map_index_name(lab)[0]
        r2 = r.get("r2")
        if r2 is None and r.get("r") is not None:
            r2 = float(r["r"]) ** 2
        r2 = float(r2) if r2 is not None else math.nan
        key = (r["property"], nm, r["form"])
        why = ""
        if r["property"] not in props:
            why = f"property {r['property']!r} not among the analysed properties {props}"
        elif nm is None:
            why = f"index label unmappable ({map_index_name(lab)[1]})"
        elif nm not in reg_names:
            why = f"index {nm} is not among the study's indices"
        elif r["form"] not in forms:
            why = f"form {r['form']!r} not in model_forms {forms} (add it to compare)"
        elif not math.isfinite(r2):
            why = "no reported r2 or r"
        elif key in rep_fit:
            why = "duplicate of an earlier reported fit (the first is compared)"
        if not why:
            rep_fit[key] = (r2, r)
        rep_rows.append(dict(property=r["property"], study_label=str(lab),
                             registry_name=nm or "", form=r["form"], reported_r2=r2,
                             compared=not why, reason=why))

    size = size_table(st)
    comp = ft.composition_features(st.mols, pd.Index(st.ids, name="id"))
    conn = st.connectivity_labels()
    groups_r = cv_groups(st, "random")
    groups_s = cv_groups(st, "scaffold") if scaffold else None
    rng = np.random.default_rng(random_seed)
    k = len(reg_names)

    fits, nested_rows, ctrl_rows = [], [], []
    for p in props:
        yv = Y[p].to_numpy(float)
        keep = np.flatnonzero(np.isfinite(yv))
        yk = yv[keep]
        # (ii)+(iii) in-sample and LOO per index x form
        for lab, name in mapped.items():
            x = ti[name].to_numpy(float)[keep]
            for d in degrees:
                f = _poly_fit(x, yk, d)
                lo = _loo(x, yk, d)
                rr, rraw = rep_fit.get((p, name, form_of[d]), (math.nan, None))
                fits.append(dict(property=p, study_label=lab, index=name, form=form_of[d],
                                 **f, reported_r2=rr, r2_minus_reported=f["r2"] - rr, **lo))
        # (iv) nested CV
        sets = {"null": pd.DataFrame(index=ti.index[keep]),
                "size": size[list(ft.SIZE_COLUMNS)].iloc[keep],
                "composition": comp.iloc[keep], "TI_study": ti.iloc[keep]}
        gk = None if groups_r is None else groups_r[keep]
        common = dict(seeds=seeds, outer_splits=outer_splits, inner_splits=inner_splits,
                      groups=gk, ad=None, n_jobs=n_jobs, min_n=min_n, min_n_test=min_n_test)
        learners = {"genre": None, "families": tuple(families)}
        for learner in learners:
            for setname, Xs in sets.items():
                fams = learners[learner]
                if learner == "genre":
                    fams = (single_index_family(Xs.shape[1], degrees),)
                base = dict(property=p, learner=learner, set=setname, split="random")
                try:
                    res = ev.nested_cv(Xs, pd.Series(yk, index=Xs.index), fams, **common)
                except ValueError as exc:
                    nested_rows.append(dict(base, reason=f"not evaluated: {exc}"[:300]))
                    continue
                nested_rows.append(dict(base, **_nested_row(res, alpha), **_oof_secondary(res),
                                        _res=res))
        if scaffold and k:
            base = dict(property=p, learner="genre", set="TI_study", split="scaffold")
            try:
                res = ev.nested_cv(sets["TI_study"], pd.Series(yk, index=sets["TI_study"].index),
                                   (single_index_family(k, degrees),),
                                   **dict(common, groups=groups_s[keep]))
                nested_rows.append(dict(base, **_nested_row(res, alpha), **_oof_secondary(res),
                                        **{f"scaffold_{a}": b for a, b in
                                           _group_stats(groups_s[keep]).items()}))
            except ValueError as exc:
                nested_rows.append(dict(base, reason=f"scaffold CV not possible: {exc}"[:300]))
        # (v) controls: observed statistic and every null from the same computation, each
        # index on the rows where it is finite (as in fits.csv)
        Xs = ti.to_numpy(float)[keep] if k else np.zeros((keep.size, 0))
        fin = np.isfinite(Xs)
        row = dict(property=p, n=int(keep.size), k=k,
                   rows_per_index=json.dumps({reg_names[j]: int(fin[:, j].sum())
                                              for j in range(k)}))
        obs, best_lab, best_n = math.nan, "", math.nan
        if k:
            b, bj, bd = _best_r2_arg(Xs, yk[:, None], degrees)
            obs = float(b[0])
            if bj[0] >= 0:
                best_lab = f"{reg_names[bj[0]]}:{form_of[int(bd[0])]}"
                best_n = int(fin[:, bj[0]].sum())
            excl = {}
            for j, nm in enumerate(reg_names):
                if not math.isfinite(float(_best_r2(Xs[:, [j]], yk[:, None], degrees)[0])):
                    excl[nm] = (f"no polynomial fit defined on its {int(fin[:, j].sum())} "
                                "finite rows (too few rows or distinct values)")
            row["excluded_indices"] = json.dumps(excl, sort_keys=True)
        fit_best = [f["r2"] for f in fits if f["property"] == p and math.isfinite(f["r2"])]
        fb = max(fit_best) if fit_best else math.nan
        agree = (math.isnan(fb) and math.isnan(obs)) or (
            math.isfinite(fb) and math.isfinite(obs) and abs(fb - obs) <= 1e-9)
        row.update(observed_best_r2=obs, observed_best=best_lab, observed_best_n=best_n,
                   observed_matches_fits=bool(agree),
                   observed_reason="" if agree else
                   f"best R2 in fits.csv {fb!r} differs from the controls' {obs!r}")
        if k:
            Z = rng.standard_normal((n_random, keep.size, k))
            Z[:, ~fin] = np.nan  # pseudo-descriptor j is missing where study index j is
            pseudo = np.array([float(_best_r2(Z[b_], yk[:, None], degrees)[0])
                               for b_ in range(n_random)])
            row.update(pseudo_null_mean=float(np.nanmean(pseudo)),
                       pseudo_null_p95=float(np.nanquantile(pseudo, 0.95)),
                       pseudo_p=_perm_p(pseudo, obs))
            Xall = all_ti.to_numpy(float)[keep]
            pool_r2 = {}
            for j, c in enumerate(all_ti.columns):
                v = float(_best_r2(Xall[:, [j]], yk[:, None], degrees)[0])
                if math.isfinite(v):
                    pool_r2[c] = v
            pool = list(pool_r2)
            n_partial = int(sum(1 for j, c in enumerate(all_ti.columns)
                                if c in pool_r2 and not np.all(np.isfinite(Xall[:, j]))))
            row.update(registry_pool=len(pool), registry_pool_partial_nan=n_partial)
            if len(pool) >= k:
                arr = np.array([pool_r2[c] for c in pool])
                reg = np.array([np.nanmax(arr[rng.choice(len(pool), k, replace=False)])
                                for _ in range(n_random)])
                row.update(registry_null_mean=float(np.nanmean(reg)),
                           registry_null_p95=float(np.nanquantile(reg, 0.95)),
                           registry_p=_perm_p(reg, obs))
            else:
                row["registry_reason"] = f"pool of {len(pool)} defined indices < k={k}"
            perms = np.stack([yk[rng.permutation(yk.size)] for _ in range(n_permutations)],
                             axis=1)
            permnull = _best_r2(Xs, perms, degrees)
            row.update(perm_null_mean=float(np.nanmean(permnull)),
                       perm_null_p95=float(np.nanquantile(permnull, 0.95)),
                       perm_p=_perm_p(permnull, obs))
        ctrl_rows.append(row)

    # ΔR² against size within each learner (identical folds; corrected t)
    for p in props:
        for learner in ("genre", "families"):
            ref = next((r for r in nested_rows if r.get("property") == p and r.get("learner") == learner
                        and r.get("set") == "size" and r.get("split") == "random" and "_res" in r), None)
            if ref is None:
                continue
            rf = ref["_res"].folds.set_index(["repeat", "fold"])
            pv = []
            for r in nested_rows:
                if (r.get("property"), r.get("learner"), r.get("split")) != (p, learner, "random") \
                        or "_res" not in r or r["set"] == "size":
                    continue
                f = r["_res"].folds.set_index(["repeat", "fold"])
                t = stats.corrected_ttest((f["r2"] - rf["r2"]).to_numpy(float),
                                          f["n_train"].to_numpy(float), f["n_test"].to_numpy(float),
                                          alpha=alpha)
                r.update(d_r2_vs_size=t["mean"], d_r2_ci_low=t["ci_low"],
                         d_r2_ci_high=t["ci_high"], d_r2_p=t["p"])
    for r in nested_rows:
        r.pop("_res", None)
    nested = pd.DataFrame(nested_rows)

    unresolved = int((compounds["status"] != "ok").sum())
    summary = dict(
        study=spec.get("name", ""), reference=spec.get("reference", ""),
        n_compounds_spec=len(spec["compounds"]), n_compounds_used=len(st.ids),
        n_compounds_not_used=unresolved, properties=props, model_forms=forms,
        n_index_labels=len(labels), n_mapped=len(reg_names),
        unmappable=[r["study_label"] for r in mrows if r["status"] != "mapped"],
        index_values_checked=int(len(index_check)),
        index_values_agree=int(index_check["agree"].sum()) if len(index_check) else 0,
        missing_property_values=yfail,
        reported_fits=len(rep_rows),
        reported_fits_unmatched=[dict(property=r["property"], study_label=r["study_label"],
                                      form=r["form"], reason=r["reason"])
                                 for r in rep_rows if not r["compared"]],
        settings=dict(families=[getattr(f, "name", f) for f in families], seeds=list(seeds),
                      outer_splits=outer_splits, inner_splits=inner_splits, alpha=alpha,
                      n_random=n_random, n_permutations=n_permutations,
                      random_seed=random_seed, scaffold=scaffold, min_n=min_n,
                      min_n_test=min_n_test),
    )
    tables = dict(compounds=compounds, index_map=index_map, index_values=values,
                  index_check=index_check, fits=pd.DataFrame(fits), nested=nested,
                  controls=pd.DataFrame(ctrl_rows),
                  reported_fits=pd.DataFrame(rep_rows, columns=[
                      "property", "study_label", "registry_name", "form", "reported_r2",
                      "compared", "reason"]))
    return ExperimentResult("D", summary, tables)


# ========================================================= size dependence
def _poly2_factory(rs):
    return Pipeline([("impute", SimpleImputer(strategy="median")),
                     ("clean", DropConstantAndDuplicate()),
                     ("poly", PolynomialFeatures(2, include_bias=False)),
                     ("scale", StandardScaler()), ("model", Ridge())])


#: Ridge on all degree-2 monomials of the (scaled) size features: captures the
#: quadratic growth of distance-based indices with size.
POLY2_FAMILY = ModelFamily("poly2_ridge", _poly2_factory,
                           {"model__alpha": [1e-4, 1e-2, 1.0, 100.0]})

#: Default learners for :func:`size_dependence` (chosen in the inner loop).
SIZE_DEPENDENCE_FAMILIES: tuple = ("linear", POLY2_FAMILY, "knn")


def size_dependence(descriptors: pd.DataFrame, size: pd.DataFrame, *,
                    family_of: Mapping[str, str] | None = None,
                    families: Sequence = SIZE_DEPENDENCE_FAMILIES, groups=None,
                    seeds: Sequence[int] = ev.DEFAULT_SEEDS, outer_splits: int = 5,
                    inner_splits: int = 3, alpha: float = 0.05, n_jobs: int = 1,
                    min_n: int = 50, min_n_test: int = 10) -> ExperimentResult:
    """How much of each descriptor *is* molecular size (EXPERIMENT_PLAN §7).

    For every column of ``descriptors`` and for the predictor sets ``size``
    (``n_heavy``, ``n_bonds``) and ``size_rings`` (+ ``n_rings``): repeated nested CV
    R² of the descriptor regressed on the set (``families`` selected in the inner loop;
    default OLS, degree-2 ridge, kNN). Rows where the descriptor is NaN are dropped
    for that descriptor (``n_missing``). ``family_of`` maps descriptor -> index family
    (default: the registry's family; other columns get ``"other"``). ``groups`` are CV
    groups aligned with the rows (connectivity groups: identical graphs have identical
    descriptors and sizes). Returns tables ``per_descriptor`` and ``per_family``.
    """
    if not descriptors.index.equals(size.index):
        raise ValueError("descriptors and size must share one index")
    fam = {s.name: s.family for s in REGISTRY}
    if family_of is not None:
        fam.update(family_of)
    sets = {"size": list(ft.SIZE_COLUMNS), "size_rings": list(ft.SIZE_RINGS_COLUMNS)}
    rows = []
    for col in descriptors.columns:
        y = descriptors[col].to_numpy(float)
        keep = np.flatnonzero(np.isfinite(y))
        for sname, cols in sets.items():
            base = dict(descriptor=col, family=fam.get(col, "other"), set=sname,
                        n_missing=int(len(y) - keep.size))
            if keep.size < 2 * outer_splits or np.ptp(y[keep]) == 0:
                why = "constant descriptor" if keep.size and np.ptp(y[keep]) == 0 else \
                    f"only {keep.size} defined values"
                rows.append(dict(base, n=int(keep.size), r2_fold_mean=math.nan, reason=why))
                continue
            Xs = size[cols].iloc[keep]
            try:
                res = ev.nested_cv(Xs, pd.Series(y[keep], index=Xs.index), families,
                                   seeds=seeds, outer_splits=outer_splits,
                                   inner_splits=inner_splits, ad=None, n_jobs=n_jobs,
                                   groups=None if groups is None else np.asarray(groups)[keep],
                                   min_n=min_n, min_n_test=min_n_test)
            except ValueError as exc:
                rows.append(dict(base, n=int(keep.size), r2_fold_mean=math.nan,
                                 reason=f"not evaluated: {exc}"[:300]))
                continue
            rows.append(dict(base, **_nested_row(res, alpha)))
    per = pd.DataFrame(rows)
    frows = []
    fam_order = list(INDEX_FAMILIES) + sorted(set(per["family"]) - set(INDEX_FAMILIES))
    for f in fam_order:
        for sname in sets:
            d = per[(per["family"] == f) & (per["set"] == sname)]
            if d.empty:
                continue
            r = d["r2_fold_mean"].to_numpy(float)
            fin = r[np.isfinite(r)]
            frows.append(dict(family=f, set=sname, n_descriptors=int(len(d)),
                              n_defined=int(fin.size),
                              r2_median=float(np.median(fin)) if fin.size else math.nan,
                              r2_mean=float(fin.mean()) if fin.size else math.nan,
                              r2_min=float(fin.min()) if fin.size else math.nan,
                              r2_max=float(fin.max()) if fin.size else math.nan,
                              n_r2_ge_0_9=int(np.sum(fin >= 0.9)),
                              n_r2_ge_0_99=int(np.sum(fin >= 0.99))))
    summary = dict(n=int(len(descriptors)), n_descriptors=int(descriptors.shape[1]),
                   settings=dict(families=[getattr(f, "name", f) for f in families],
                                 seeds=list(seeds), outer_splits=outer_splits,
                                 inner_splits=inner_splits, grouped=groups is not None))
    return ExperimentResult("size_dependence", summary,
                            dict(per_descriptor=per, per_family=pd.DataFrame(frows)))
