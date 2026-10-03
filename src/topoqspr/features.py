"""Baseline feature sets for the evaluation ladder.

Every feature set is computed from an RDKit ``Mol`` with hydrogen-suppressed
graph semantics: the molecular graph G = (V, E) has one vertex per heavy atom
(atomic number > 1) and one edge per bond between two heavy atoms. Explicit
hydrogen atoms in the input are ignored for the graph and counted only in the
composition set.

Feature sets (columns, in this order):

``size``
    ``n_heavy`` = |V|, ``n_bonds`` = |E|.
``size_rings``
    ``size`` + ``n_rings`` = cyclomatic number mu(G) = |E| - |V| + c(G), where
    c(G) is the number of connected components. mu(G) equals the size of the
    smallest set of smallest rings (SSSR, ``len(Chem.GetSSSR(mol))``); note
    that ``rdMolDescriptors.CalcNumRings`` counts the *symmetrised* SSSR and
    can be larger (cubane: 6 vs mu = 5). Because mu is an exact linear
    function of |V|, |E| and c, ``size_rings`` adds information beyond ``size``
    only through c(G) (multi-fragment inputs) for a linear model; tree models
    may still use it differently. Reported as a sensitivity check.
``composition``
    counts of C, N, O, S, P, F, Cl, Br, I, all other heavy atoms
    (``n_other``) and all hydrogens (implicit + explicit, ``n_H``).
``rdkit2d``
    RDKit ``Descriptors.descList`` minus a per-target exclusion list
    (:data:`TARGET_EXCLUSIONS`). A descriptor is excluded when it *is* the
    target, is computed from the same per-atom contributions as the target,
    or is a function that takes the target as an input (e.g. QED uses MW,
    Crippen logP, TPSA, HBD, HBA and rotatable bonds; Bickerton et al. 2012).
    Alternative estimators of the same physical quantity (e.g. Crippen
    ``MolLogP`` when the target is PubChem XLogP3 or an experimental logP)
    are *kept*: a practitioner's conventional baseline would use them.

Topological-index features are produced by a separate module and are passed
to the evaluation as a DataFrame; this module does not compute them.

Failures are never imputed. A value that cannot be computed, or that is not
finite, is stored as NaN and a record ``{"row", "feature", "reason"}`` is
appended to ``DataFrame.attrs["failures"]``.

References
----------
Wildman SA, Crippen GM (1999) J Chem Inf Comput Sci 39:868-873
(atomic logP and MR contributions). doi:10.1021/ci990307l
Ertl P, Rohde B, Selzer P (2000) J Med Chem 43:3714-3717 (TPSA).
doi:10.1021/jm000942e
Bickerton GR et al. (2012) Nat Chem 4:90-98 (QED). doi:10.1038/nchem.1243
"""

from __future__ import annotations

import math
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

__all__ = [
    "SIZE_COLUMNS",
    "SIZE_RINGS_COLUMNS",
    "COMPOSITION_COLUMNS",
    "TARGET_EXCLUSIONS",
    "rdkit2d_names",
    "size_features",
    "size_rings_features",
    "composition_features",
    "rdkit2d_features",
    "build_feature_sets",
]

SIZE_COLUMNS = ("n_heavy", "n_bonds")
SIZE_RINGS_COLUMNS = SIZE_COLUMNS + ("n_rings",)
_ELEMENTS = ("C", "N", "O", "S", "P", "F", "Cl", "Br", "I")
COMPOSITION_COLUMNS = tuple(f"n_{e}" for e in _ELEMENTS) + ("n_other", "n_H")

# --- descriptor families sharing per-atom contributions ---------------------
_CRIPPEN_LOGP = ("MolLogP", "BCUT2D_LOGPHI", "BCUT2D_LOGPLOW") + tuple(
    f"SlogP_VSA{i}" for i in range(1, 13)
)
_CRIPPEN_MR = ("MolMR", "BCUT2D_MRHI", "BCUT2D_MRLOW") + tuple(
    f"SMR_VSA{i}" for i in range(1, 11)
)
# Both Crippen families use the same Wildman-Crippen atom typing, so either
# Crippen target excludes both families.
_CRIPPEN = _CRIPPEN_LOGP + _CRIPPEN_MR
_MASS = ("MolWt", "ExactMolWt", "HeavyAtomMolWt", "BCUT2D_MWHI", "BCUT2D_MWLOW")
_QED = ("qed",)

#: target name -> descriptors removed from ``rdkit2d`` for that target.
#: ``None`` (no exclusion) must be requested explicitly, as ``target=None``.
TARGET_EXCLUSIONS: dict[str, tuple[str, ...]] = {
    # mass (RDKit and PubChem spellings are the same quantity)
    "MolWt": _MASS + _QED,
    "PC_MolecularWeight": _MASS + _QED,
    "PC_ExactMass": _MASS + _QED,
    "PC_MonoisotopicMass": _MASS + _QED,
    # Crippen logP / MR (Wildman & Crippen 1999)
    "LogP": _CRIPPEN + _QED,
    "MolLogP": _CRIPPEN + _QED,
    "MolMR": _CRIPPEN,
    # Ertl TPSA. LabuteASA is kept: it is total (not polar) approximate
    # surface area from Labute's per-atom contributions, not Ertl's polar
    # fragment contributions, and no other descList entry uses TPSA terms.
    "TPSA": ("TPSA",) + _QED,
    "PC_TPSA": ("TPSA",) + _QED,
    # Lipinski-style counts
    "HBD": ("NumHDonors", "NHOHCount") + _QED,
    "PC_HBondDonorCount": ("NumHDonors", "NHOHCount") + _QED,
    "HBA": ("NumHAcceptors", "NOCount") + _QED,
    "PC_HBondAcceptorCount": ("NumHAcceptors", "NOCount") + _QED,
    "RotBonds": ("NumRotatableBonds",) + _QED,
    "PC_RotatableBondCount": ("NumRotatableBonds",) + _QED,
    # PubChem complexity is the Bertz/Hendrickson/Ihlenfeldt formula.
    "PC_Complexity": ("BertzCT",),
    # Alternative estimators / measurements of logP: Crippen is a separately
    # fitted model, kept as a legitimate conventional baseline.
    "PC_XLogP": (),
    "XLogP3": (),
}


def _heavy_graph(mol: Chem.Mol) -> tuple[int, int, int]:
    """Return (|V|, |E|, c(G)) of the hydrogen-suppressed graph."""
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    heavy_set = set(heavy)
    edges = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in heavy_set and b.GetEndAtomIdx() in heavy_set
    ]
    # union-find for connected components
    parent = {v: v for v in heavy}

    def find(v: int) -> int:
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
    n_comp = len({find(v) for v in heavy})
    return len(heavy), len(edges), n_comp


def _index(mols: Sequence, index: Iterable | None) -> pd.Index:
    if index is None:
        return pd.RangeIndex(len(mols))
    idx = pd.Index(index)
    if len(idx) != len(mols):
        raise ValueError(f"index has {len(idx)} labels for {len(mols)} molecules")
    return idx


def _frame(rows: list[dict], columns: Sequence[str], index: pd.Index,
           failures: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=list(columns), index=index, dtype=float)
    df.attrs["failures"] = failures
    return df


def _none_row(columns: Sequence[str], label, failures: list[dict], reason: str) -> dict:
    for c in columns:
        failures.append({"row": label, "feature": c, "reason": reason})
    return {c: math.nan for c in columns}


def size_features(mols: Sequence[Chem.Mol | None], index=None) -> pd.DataFrame:
    """``n_heavy`` = |V| and ``n_bonds`` = |E| of the hydrogen-suppressed graph."""
    idx = _index(mols, index)
    rows, failures = [], []
    for label, mol in zip(idx, mols):
        if mol is None:
            rows.append(_none_row(SIZE_COLUMNS, label, failures, "mol is None"))
            continue
        n, m, _ = _heavy_graph(mol)
        rows.append({"n_heavy": n, "n_bonds": m})
    return _frame(rows, SIZE_COLUMNS, idx, failures)


def size_rings_features(mols: Sequence[Chem.Mol | None], index=None) -> pd.DataFrame:
    """``size`` plus ``n_rings`` = |E| - |V| + c(G) (cyclomatic number)."""
    idx = _index(mols, index)
    rows, failures = [], []
    for label, mol in zip(idx, mols):
        if mol is None:
            rows.append(_none_row(SIZE_RINGS_COLUMNS, label, failures, "mol is None"))
            continue
        n, m, c = _heavy_graph(mol)
        rows.append({"n_heavy": n, "n_bonds": m, "n_rings": m - n + c})
    return _frame(rows, SIZE_RINGS_COLUMNS, idx, failures)


def composition_features(mols: Sequence[Chem.Mol | None], index=None) -> pd.DataFrame:
    """Element counts of heavy atoms plus the total hydrogen count ``n_H``.

    ``n_H`` = number of H atoms after ``Chem.AddHs`` (implicit and explicit
    hydrogens; isotopes such as D count as H).
    """
    idx = _index(mols, index)
    rows, failures = [], []
    for label, mol in zip(idx, mols):
        if mol is None:
            rows.append(_none_row(COMPOSITION_COLUMNS, label, failures, "mol is None"))
            continue
        row = dict.fromkeys(COMPOSITION_COLUMNS, 0)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() <= 1:
                continue
            key = f"n_{atom.GetSymbol()}"
            row[key if key in row and key != "n_H" else "n_other"] += 1
        try:
            withh = Chem.AddHs(mol)
            row["n_H"] = sum(1 for a in withh.GetAtoms() if a.GetAtomicNum() == 1)
        except Exception as exc:  # noqa: BLE001 - recorded, not swallowed
            row["n_H"] = math.nan
            failures.append({"row": label, "feature": "n_H", "reason": repr(exc)})
        rows.append(row)
    return _frame(rows, COMPOSITION_COLUMNS, idx, failures)


def rdkit2d_names(target: str | None, extra_exclude: Iterable[str] = ()) -> list[str]:
    """Descriptor names kept for ``target`` (``descList`` order).

    ``target=None`` keeps every descriptor; an unknown target name raises
    ``KeyError`` so that each new target gets an explicit exclusion decision.
    """
    if target is None:
        excluded = set()
    else:
        if target not in TARGET_EXCLUSIONS:
            raise KeyError(
                f"no rdkit2d exclusion rule for target {target!r}; add it to "
                f"TARGET_EXCLUSIONS or pass target=None explicitly"
            )
        excluded = set(TARGET_EXCLUSIONS[target])
    excluded |= set(extra_exclude)
    return [name for name, _ in Descriptors.descList if name not in excluded]


def rdkit2d_features(mols: Sequence[Chem.Mol | None], target: str | None,
                     index=None, extra_exclude: Iterable[str] = ()) -> pd.DataFrame:
    """RDKit 2D descriptors minus the exclusions for ``target``.

    Descriptors are computed on ``Chem.RemoveHs(mol)``. Exceptions and
    non-finite values become NaN with a failure record.
    """
    idx = _index(mols, index)
    names = rdkit2d_names(target, extra_exclude)
    funcs = dict(Descriptors.descList)
    rows, failures = [], []
    for label, mol in zip(idx, mols):
        if mol is None:
            rows.append(_none_row(names, label, failures, "mol is None"))
            continue
        try:
            hmol = Chem.RemoveHs(mol)
        except Exception as exc:  # noqa: BLE001
            rows.append(_none_row(names, label, failures, f"RemoveHs: {exc!r}"))
            continue
        row = {}
        for name in names:
            try:
                val = float(funcs[name](hmol))
            except Exception as exc:  # noqa: BLE001
                row[name] = math.nan
                failures.append({"row": label, "feature": name, "reason": repr(exc)})
                continue
            if not math.isfinite(val):
                failures.append({"row": label, "feature": name,
                                 "reason": f"non-finite value {val!r}"})
                val = math.nan
            row[name] = val
        rows.append(row)
    return _frame(rows, names, idx, failures)


def build_feature_sets(mols: Sequence[Chem.Mol | None], target: str | None,
                       index=None,
                       ti: pd.DataFrame | None = None,
                       extra_sets: Mapping[str, pd.DataFrame] | None = None,
                       ) -> dict[str, pd.DataFrame]:
    """Return the baseline ladder as ``{name: DataFrame}`` with a shared index.

    Keys: ``null`` (zero columns; evaluated with the training-fold mean),
    ``size``, ``size_rings``, ``composition``, ``rdkit2d`` and, when given,
    ``TI`` (``ti`` reindexed to the molecule index; rows missing from ``ti``
    become NaN with a failure record) plus any ``extra_sets``.
    """
    idx = _index(mols, index)
    sets: dict[str, pd.DataFrame] = {
        "null": pd.DataFrame(index=idx),
        "size": size_features(mols, idx),
        "size_rings": size_rings_features(mols, idx),
        "composition": composition_features(mols, idx),
        "rdkit2d": rdkit2d_features(mols, target, idx),
    }
    if ti is not None:
        missing = idx.difference(ti.index)
        t = ti.reindex(idx).astype(float)
        failures = [{"row": r, "feature": "*", "reason": "row absent from TI frame"}
                    for r in missing]
        inf_r, inf_c = np.nonzero(np.isinf(t.to_numpy()))
        failures += [{"row": idx[i], "feature": t.columns[j],
                      "reason": "non-finite value"} for i, j in zip(inf_r, inf_c)]
        t = t.replace([np.inf, -np.inf], np.nan)
        t.attrs["failures"] = list(ti.attrs.get("failures", [])) + failures
        sets["TI"] = t
    for name, df in (extra_sets or {}).items():
        if name in sets:
            raise ValueError(f"feature set name {name!r} already used")
        if not df.index.equals(idx):
            raise ValueError(f"feature set {name!r} index does not match molecules")
        sets[name] = df
    return sets
