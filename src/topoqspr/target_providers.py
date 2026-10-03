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

Shipped providers:

* ``rdkit`` (:class:`RDKitTargetProvider`) - computed targets Crippen logP,
  Crippen molar refractivity, Ertl TPSA (class 2, "computed
  physicochemical": genre replication only, no headline claims) and
  molecular weight (class 1, sanity check only);
* ``expdata`` (:class:`ExpdataTargetProvider`) - the licence-checked
  measured targets of :mod:`topoqspr.expdata` (class 3), joined on the full
  standard InChIKey of the standardised parent.

Other providers plug in through :func:`register_target_provider` or by
naming ``"package.module:factory"`` as the provider in the study config; the
factory is called as ``factory(names=<list or None>, **options)``. A provider
may also define ``set_runtime(**runtime)``: the pipeline calls it before
``compute`` with machine-local settings that are not part of the hashed
config (``cache_dir``, ``chembl_cache_dir``, ``offline``, ``client``).

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
    "ExpdataTargetProvider",
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


# ------------------------------------------------------------- expdata
def _json_safe(x):
    """NaN/inf -> None, numpy scalars -> Python (the run record forbids NaN)."""
    if isinstance(x, Mapping):
        return {str(k): _json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_json_safe(v) for v in x]
    if hasattr(x, "item") and not isinstance(x, (str, bytes)):
        x = x.item()
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


class ExpdataTargetProvider:
    """Measured (class 3) targets from :func:`topoqspr.expdata.get_target`.

    Matching rule: a molecule gets the value of the ``get_target`` row whose
    ``inchikey`` equals the **full** standard InChIKey of its standardised
    parent (``molecules["inchikey"]``); there is no first-block fallback, so a
    stereoisomer never inherits its partner's measurement. Every molecule
    without a value is NaN with one of these reasons:

    * ``no_inchikey`` - the standardised parent has no InChIKey;
    * ``no_measurement`` - no curated record of this target has the full key;
    * ``no_measurement_stereo_mismatch`` - no full-key match, but a record of
      another stereoisomer (same first block) exists; it is not used;
    * ``replicates_disagree`` - the parent's replicate SD exceeds the target's
      ``disagreement_sd`` (the ``get_target`` rule, applied here so that the
      reason is recorded; ``drop_disagreeing=False`` keeps such values).

    Options (JSON, part of the hashed config): ``drop_disagreeing`` (default
    True), ``keep_salt_forms`` (default False; sensitivity analysis),
    ``chembl_release`` (release the AstraZeneca ChEMBL records are pinned to;
    default: the release of the client passed by the pipeline, else
    ``"ChEMBL_37"``). Machine-local settings come through :meth:`set_runtime`.

    ``loader`` (Python API only) replaces :func:`topoqspr.expdata.get_target`,
    called as ``loader(name, cache_dir, offline=..., drop_disagreeing=False,
    keep_salt_forms=..., client=...)``; the offline tests use it.
    """

    name = "expdata"

    def __init__(self, names: Sequence[str] | None = None, *, drop_disagreeing: bool = True,
                 keep_salt_forms: bool = False, chembl_release: str | None = None,
                 loader: Callable[..., pd.DataFrame] | None = None) -> None:
        from . import expdata

        self._ed = expdata
        available = tuple(expdata.TARGETS)
        names = list(available) if names is None else [str(n) for n in names]
        unknown = [n for n in names if n not in expdata.TARGETS]
        if unknown:
            raise KeyError(f"unknown expdata targets {unknown}; available: {list(available)}")
        if len(set(names)) != len(names):
            raise ValueError(f"duplicate target names in {names}")
        self.names = names
        self.drop_disagreeing = bool(drop_disagreeing)
        self.keep_salt_forms = bool(keep_salt_forms)
        self.chembl_release = chembl_release
        self._loader = loader
        self._runtime: dict[str, Any] = {}
        self._prov: dict[str, dict] = {}

    def set_runtime(self, *, cache_dir=None, chembl_cache_dir=None, offline: bool = False,
                    client=None, **_ignored) -> None:
        """Machine-local settings (not hashed): the expdata cache directory
        (default ``<topoqspr cache>/expdata``), the ChEMBL cache and offline
        mode for the AstraZeneca (ChEMBL-hosted) targets, or a ready client."""
        self._runtime = dict(cache_dir=cache_dir, chembl_cache_dir=chembl_cache_dir,
                             offline=bool(offline), client=client)

    def specs(self) -> list[TargetSpec]:
        out = []
        for n in self.names:
            t = self._ed.TARGETS[n]
            src = self._ed.SOURCES[t.source]
            # Measured values are not computed from any RDKit descriptor; alternative
            # estimators of the same quantity (e.g. Crippen logP for measured logP) are
            # kept by the leakage rule (features.TARGET_EXCLUSIONS docstring).
            out.append(TargetSpec(n, 3, None, t.unit, f"expdata:{t.source}",
                                  f"{t.endpoint}; {t.description}; licence {src.licence}"))
        return out

    def _client(self):
        rt = self._runtime
        if rt.get("client") is not None:
            if self.chembl_release and rt["client"].release != self.chembl_release:
                raise RuntimeError(f"expdata provider pins {self.chembl_release} but the client "
                                   f"serves {rt['client'].release}")
            return rt["client"]
        from .chembl import ChemblClient

        return ChemblClient(rt.get("chembl_cache_dir"), offline=bool(rt.get("offline")),
                            release=self.chembl_release or "ChEMBL_37")

    def _load(self, name: str) -> pd.DataFrame:
        t = self._ed.TARGETS[name]
        src = self._ed.SOURCES[t.source]
        client = self._client() if src.chembl_assay else None
        loader = self._loader or self._ed.get_target
        return loader(name, self._runtime.get("cache_dir"), offline=bool(self._runtime.get("offline")),
                      drop_disagreeing=False, keep_salt_forms=self.keep_salt_forms, client=client)

    def compute(self, molecules: pd.DataFrame, mols: Sequence[Chem.Mol]) -> TargetValues:
        if len(mols) != len(molecules):
            raise ValueError(f"{len(mols)} molecules for {len(molecules)} table rows")
        keys = molecules["inchikey"].fillna("").astype(str).tolist()
        cols, failures = {}, []
        self._prov = {}
        for n in self.names:
            t = self._ed.TARGETS[n]
            src = self._ed.SOURCES[t.source]
            df = self._load(n)
            missing = [c for c in ("inchikey", "value", "replicates_disagree") if c not in df.columns]
            if missing:
                raise ValueError(f"expdata target {n!r}: loader returned no {missing} column(s)")
            if df["inchikey"].duplicated().any():
                raise ValueError(f"expdata target {n!r}: duplicate InChIKeys in get_target output")
            sds = df["sd"] if "sd" in df.columns else pd.Series(math.nan, index=df.index)
            by_key = {str(k): (float(v), bool(d), float(s)) for k, v, d, s in zip(
                df["inchikey"], df["value"], df["replicates_disagree"], sds)}
            blocks = {k.split("-")[0] for k in by_key}
            vals, n_match, n_dis = [], 0, 0
            for label, key in zip(molecules.index, keys):
                v, why = math.nan, ""
                if not key:
                    why = "no_inchikey: standardised parent has no InChIKey"
                elif key in by_key:
                    val, dis, sd = by_key[key]
                    n_match += 1
                    if dis and self.drop_disagreeing:
                        n_dis += 1
                        why = (f"replicates_disagree: replicate SD {sd:.3g} > "
                               f"{t.disagreement_sd:g} {t.unit}")
                    elif not math.isfinite(val):
                        why = f"non-finite measured value {val!r}"
                    else:
                        v = val
                elif key.split("-")[0] in blocks:
                    why = ("no_measurement_stereo_mismatch: only another stereoisomer or "
                           "isotopologue (same InChIKey first block) has a value; not used "
                           "(full-key matching)")
                else:
                    why = f"no_measurement: no curated {n} record with this full InChIKey"
                if why:
                    failures.append({"row": label, "feature": n, "reason": why})
                vals.append(v)
            cols[n] = vals
            self._prov[n] = dict(
                target=n, endpoint=t.endpoint, unit=t.unit, target_class=3, source=src.key,
                title=src.title, citation=src.citation, doi=src.doi, version=src.version,
                licence=src.licence, redistribution=src.redistribution,
                url=src.url or None, chembl_assay=src.chembl_assay or None,
                raw_sha256=df.attrs.get("raw_sha256"),
                processing_fingerprint=df.attrs.get("processing_fingerprint"),
                n_parents_in_source=int(len(df)), n_matched_full_key=int(n_match),
                n_dropped_disagreeing=int(n_dis),
                n_values=int(sum(math.isfinite(x) for x in vals)))
        values = pd.DataFrame(cols, index=molecules.index, columns=self.names, dtype=float)
        return TargetValues(values, failures)

    def provenance(self) -> dict:
        return _json_safe(dict(
            provider=self.name,
            matching="full standard InChIKey of the standardised parent "
                     "(topoqspr.graph.standardize); no first-block fallback",
            drop_disagreeing=self.drop_disagreeing, keep_salt_forms=self.keep_salt_forms,
            chembl_release=self.chembl_release,
            processing_version=self._ed.PROCESSING_VERSION,
            targets=[self._prov.get(n, dict(target=n, computed=False)) for n in self.names]))


# --------------------------------------------------------------- registry
#: provider name -> factory(names=None | list, **options) -> TargetProvider
PROVIDERS: dict[str, Callable[..., TargetProvider]] = {
    "rdkit": lambda names=None, **opts: RDKitTargetProvider(names, **opts),
    "expdata": lambda names=None, **opts: ExpdataTargetProvider(names, **opts),
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
