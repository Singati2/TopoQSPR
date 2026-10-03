"""Target providers: where the modelled property values come from.

The pipeline asks each configured provider for

* :meth:`TargetProvider.specs` - one :class:`TargetSpec` per target column
  (name, EXPERIMENT_PLAN §4 class, the ``rdkit2d`` exclusion rule, units,
  source), and
* :meth:`TargetProvider.compute` - the values for a table of standardised
  molecules, as a :class:`TargetValues` (NaN + a failure record wherever a
  value is missing; never 0), and
* :meth:`TargetProvider.provenance` - what the run record must say about the
  source (version, licence, file hashes, ...).

The molecule table passed to ``compute`` is indexed by molecule ID (ChEMBL
parent ID or the ``id`` column of a SMILES file) with columns
``input_smiles`` (as supplied), ``smiles`` (canonical SMILES of the
standardised parent, :func:`topoqspr.graph.standardize`) and ``inchikey``
(full standard InChIKey of that parent). Experimental providers match on
``inchikey`` and must document their matching rule (full key vs first
block) in their ``provenance()``. The RDKit ``Mol`` objects of the
standardised parents are passed alongside, in the same order.

Shipped provider: ``rdkit`` (:class:`RDKitTargetProvider`) - computed
targets Crippen logP, Crippen molar refractivity, Ertl TPSA (class 2,
"computed physicochemical": genre replication only, no headline claims) and
molecular weight (class 1, sanity check only). Experimental (class 3)
providers plug in through :func:`register_target_provider` or by naming
``"package.module:factory"`` as the provider in the study config; the
factory is called as ``factory(names=<list or None>, **options)``.

References
----------
Wildman SA, Crippen GM (1999) J Chem Inf Comput Sci 39:868-873. doi:10.1021/ci990307l
Ertl P, Rohde B, Selzer P (2000) J Med Chem 43:3714-3717. doi:10.1021/jm000942e
"""

from __future__ import annotations

import importlib
import math
from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Protocol, Sequence, runtime_checkable

import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import Crippen, Descriptors

from .features import TARGET_EXCLUSIONS

__all__ = [
    "TARGET_CLASSES",
    "TargetSpec",
    "TargetValues",
    "TargetProvider",
    "RDKitTargetProvider",
    "register_target_provider",
    "get_target_provider",
    "PROVIDERS",
]

#: EXPERIMENT_PLAN §4 target classes.
TARGET_CLASSES: dict[int, str] = {
    1: "structural / near-deterministic (sanity check only)",
    2: "computed physicochemical (genre replication; no headline claims)",
    3: "external experimental (headline analysis)",
}


@dataclass(frozen=True)
class TargetSpec:
    """One target column.

    ``rdkit2d_exclusion`` is the key into
    :data:`topoqspr.features.TARGET_EXCLUSIONS` that decides which RDKit
    descriptors are removed from the ``rdkit2d`` baseline for this target;
    ``None`` means "exclude nothing" and has to be written on purpose (there
    is no default), so that every new target gets an explicit leakage
    decision.
    """

    name: str
    target_class: int
    rdkit2d_exclusion: str | None
    units: str
    source: str
    description: str = ""

    def __post_init__(self) -> None:
        if self.target_class not in TARGET_CLASSES:
            raise ValueError(f"target {self.name!r}: class must be one of {sorted(TARGET_CLASSES)}")
        if self.rdkit2d_exclusion is not None and self.rdkit2d_exclusion not in TARGET_EXCLUSIONS:
            raise KeyError(f"target {self.name!r}: no rdkit2d exclusion rule "
                           f"{self.rdkit2d_exclusion!r} in features.TARGET_EXCLUSIONS")

    def as_dict(self) -> dict:
        return dict(name=self.name, target_class=self.target_class,
                    target_class_label=TARGET_CLASSES[self.target_class],
                    rdkit2d_exclusion=self.rdkit2d_exclusion, units=self.units,
                    source=self.source, description=self.description)


@dataclass
class TargetValues:
    """Values of a provider's targets: ``values`` has the molecule index and
    one float column per spec (NaN = missing); ``failures`` holds
    ``{"row", "feature", "reason"}`` records for every NaN."""

    values: pd.DataFrame
    failures: list[dict] = field(default_factory=list)


@runtime_checkable
class TargetProvider(Protocol):
    """Interface every target source implements (see the module docstring)."""

    name: str

    def specs(self) -> list[TargetSpec]:
        ...

    def compute(self, molecules: pd.DataFrame, mols: Sequence[Chem.Mol]) -> TargetValues:
        ...

    def provenance(self) -> dict:
        ...


# ----------------------------------------------------------------- RDKit
_RDKIT_TARGETS: dict[str, tuple[int, str, str, Callable[[Chem.Mol], float], str]] = {
    # name: (class, exclusion key, units, function, description)
    "MolWt": (1, "MolWt", "g/mol", Descriptors.MolWt,
              "average molecular weight (RDKit Descriptors.MolWt)"),
    "MolLogP": (2, "MolLogP", "log10 units", Crippen.MolLogP,
                "Wildman-Crippen logP (RDKit Crippen.MolLogP)"),
    "TPSA": (2, "TPSA", "angstrom^2", Descriptors.TPSA,
             "Ertl topological polar surface area, N and O only (RDKit Descriptors.TPSA)"),
    "MolMR": (2, "MolMR", "cm^3/mol", Crippen.MolMR,
              "Wildman-Crippen molar refractivity (RDKit Crippen.MolMR)"),
}


class RDKitTargetProvider:
    """Computed targets from RDKit, on ``Chem.RemoveHs`` of the standardised parent
    (the same molecule and H handling as :func:`topoqspr.features.rdkit2d_features`)."""

    name = "rdkit"
    AVAILABLE: tuple[str, ...] = tuple(_RDKIT_TARGETS)

    def __init__(self, names: Sequence[str] | None = None) -> None:
        names = list(self.AVAILABLE) if names is None else [str(n) for n in names]
        unknown = [n for n in names if n not in _RDKIT_TARGETS]
        if unknown:
            raise KeyError(f"unknown RDKit targets {unknown}; available: {list(self.AVAILABLE)}")
        if len(set(names)) != len(names):
            raise ValueError(f"duplicate target names in {names}")
        self.names = names

    def specs(self) -> list[TargetSpec]:
        out = []
        for n in self.names:
            cls, excl, units, _, desc = _RDKIT_TARGETS[n]
            out.append(TargetSpec(n, cls, excl, units, "rdkit", desc))
        return out

    def compute(self, molecules: pd.DataFrame, mols: Sequence[Chem.Mol]) -> TargetValues:
        if len(mols) != len(molecules):
            raise ValueError(f"{len(mols)} molecules for {len(molecules)} table rows")
        rows, failures = [], []
        for label, mol in zip(molecules.index, mols):
            row = {}
            try:
                with rdBase.BlockLogs():
                    hmol = Chem.RemoveHs(mol)
            except Exception as exc:  # noqa: BLE001 - recorded, never imputed
                for n in self.names:
                    row[n] = math.nan
                    failures.append({"row": label, "feature": n, "reason": f"RemoveHs: {exc!r}"})
                rows.append(row)
                continue
            for n in self.names:
                fn = _RDKIT_TARGETS[n][3]
                try:
                    v = float(fn(hmol))
                except Exception as exc:  # noqa: BLE001
                    row[n] = math.nan
                    failures.append({"row": label, "feature": n, "reason": repr(exc)})
                    continue
                if not math.isfinite(v):
                    failures.append({"row": label, "feature": n, "reason": f"non-finite value {v!r}"})
                    v = math.nan
                row[n] = v
            rows.append(row)
        df = pd.DataFrame(rows, index=molecules.index, columns=self.names, dtype=float)
        return TargetValues(df, failures)

    def provenance(self) -> dict:
        return dict(provider=self.name, rdkit_version=rdBase.rdkitVersion,
                    computed_on="Chem.RemoveHs(standardised parent)",
                    targets=list(self.names))


# --------------------------------------------------------------- registry
#: provider name -> factory(names=None | list, **options) -> TargetProvider
PROVIDERS: dict[str, Callable[..., TargetProvider]] = {
    "rdkit": lambda names=None, **opts: RDKitTargetProvider(names, **opts),
}


def register_target_provider(name: str, factory: Callable[..., TargetProvider]) -> None:
    """Make ``factory`` available as ``{"provider": name}`` in study configs."""
    if name in PROVIDERS:
        raise ValueError(f"target provider {name!r} already registered")
    PROVIDERS[name] = factory


def get_target_provider(name: str, names: Sequence[str] | None = None,
                        options: Mapping[str, Any] | None = None) -> TargetProvider:
    """Instantiate a registered provider, or ``"module.path:attr"`` (imported)."""
    options = dict(options or {})
    if name in PROVIDERS:
        factory = PROVIDERS[name]
    elif ":" in name:
        mod, _, attr = name.partition(":")
        factory = getattr(importlib.import_module(mod), attr)
    else:
        raise KeyError(f"unknown target provider {name!r}; registered: {sorted(PROVIDERS)} "
                       "(or give 'package.module:factory')")
    prov = factory(names=None if names is None else list(names), **options)
    if not isinstance(prov, TargetProvider):
        raise TypeError(f"target provider {name!r} does not implement TargetProvider")
    return prov
