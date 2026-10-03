"""Structure standardisation and the hydrogen-suppressed molecular graph.

A SMILES string is turned into a simple, undirected, connected graph
``G = (V, E)`` whose vertices are the heavy (non-hydrogen) atoms of the
standardised parent molecule and whose edges are the bonds between them.
Bond orders, charges, aromaticity and element types are deliberately
discarded: every index in :mod:`topoqspr.indices` is a graph invariant of
this hydrogen-suppressed skeleton, as in the chemical-graph-theory papers
the indices come from.

Standardisation (in this order)
-------------------------------
1. Reject input containing whitespace after trimming (``unparsable``):
   RDKit reads text after the first space as the molecule name, so
   ``"CCO junk"`` would otherwise be accepted as ethanol. CXSMILES
   extensions (``"CCO |$...$|"``) are therefore not accepted either.
2. Parse with RDKit (``Chem.MolFromSmiles``, default sanitisation).
3. ChEMBL Structure Pipeline: ``standardize_mol`` then ``get_parent_mol``
   (Bento et al., J. Cheminform. 12 (2020) 51). This removes counter-ions,
   salts and solvents listed in the pipeline's dictionaries, neutralises
   charges where possible, and strips isotope labels (``[2H]`` becomes H).
4. Parent fallback, only if more than one fragment survives step 3. The
   pipeline keeps *every* fragment when all of them are on its salt list
   (sodium benzoate, choline chloride, ammonium stearate), and its list
   lacks a few pharmaceutical counter-ions. So, in order:

   a. drop fragments without heavy atoms (``[H+]``);
   b. drop monoatomic inorganic ions and ammonia/ammonium (one heavy atom
      from :data:`SIMPLE_COUNTERION_ELEMENTS`);
   c. drop fragments in :data:`SUPPLEMENTARY_COUNTERIONS` (matched on the
      neutralised canonical SMILES);
   d. merge identical fragments.

   Steps a-c are each skipped if they would remove every fragment. If
   exactly one fragment is left and it is organic (some carbon bonded to H
   or C, so lithium carbonate and sodium cyanide are left alone and later
   rejected as ``disallowed_element``), it is neutralised with the pipeline's
   ``uncharge_mol`` and becomes the parent; what was dropped is recorded in
   ``StandardizationResult.parent_note``. Otherwise the input is rejected
   (``parent_ambiguous`` when every remaining fragment is itself on the
   ChEMBL salt list, e.g. choline salicylate; ``mixture`` otherwise).
5. Reject the result if

   * it contains an element outside ``allowed_elements``
     (default :data:`ORGANIC_ELEMENTS`)  -> ``disallowed_element``
   * it contains no carbon                -> ``no_carbon``
   * more than one fragment remains        -> ``mixture`` or
     ``parent_ambiguous`` (step 4)
   * it has fewer than ``min_heavy_atoms`` heavy atoms (default 3)
                                           -> ``too_few_heavy_atoms``
   * the heavy-atom graph is disconnected (possible only when hydrogens
     bridge heavy atoms)                   -> ``disconnected_graph``

Nothing is silently repaired: every rejection carries a machine-readable
``code`` and a human-readable ``detail``.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from functools import cached_property, lru_cache
from typing import Iterable, Optional, Sequence

from chembl_structure_pipeline import standardizer as _csp
from rdkit import Chem, rdBase

__all__ = [
    "ORGANIC_ELEMENTS",
    "SIMPLE_COUNTERION_ELEMENTS",
    "SUPPLEMENTARY_COUNTERIONS",
    "REJECTION_CODES",
    "MolGraph",
    "Rejection",
    "StandardizationResult",
    "standardize",
    "graph_from_mol",
]

#: Elements accepted by default (hydrogen-suppressed organic chemistry).
ORGANIC_ELEMENTS: frozenset[str] = frozenset(
    {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I"}
)

#: Elements whose single-heavy-atom fragments are dropped as counter-ions in
#: the parent fallback: alkali and alkaline-earth metals, halogens, and N
#: (ammonia / ammonium).
SIMPLE_COUNTERION_ELEMENTS: frozenset[str] = frozenset(
    {"Li", "Na", "K", "Rb", "Cs", "Mg", "Ca", "Sr", "Ba", "F", "Cl", "Br", "I", "N"}
)

#: Pharmaceutical counter-ions missing from the ChEMBL 1.2.4 salt list, as
#: neutral RDKit canonical SMILES -> name. Of 24 common counter-ions tested
#: as fluoxetine salts (besylate, camsylate, edisylate, esylate, isethionate,
#: napadisylate, closylate, hippurate, lactobionate, gluceptate, lauryl
#: sulfate, orotate, embonate, oxalate, adipate, benzoate, salicylate,
#: olamine, diolamine, trolamine, benzathine, cyclamate, stearate,
#: xinafoate) only xinafoate was not stripped.
SUPPLEMENTARY_COUNTERIONS: dict[str, str] = {
    "O=C(O)c1cc2ccccc2cc1O": "xinafoic acid (1-hydroxy-2-naphthoic acid)",
}

#: Every rejection code :func:`standardize` can emit.
REJECTION_CODES: tuple[str, ...] = (
    "unparsable",
    "standardization_failed",
    "disallowed_element",
    "no_carbon",
    "mixture",
    "parent_ambiguous",
    "too_few_heavy_atoms",
    "disconnected_graph",
)


# ---------------------------------------------------------------------------
# Graph
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class MolGraph:
    """A simple undirected graph on vertices ``0 .. n-1``.

    ``edges`` is normalised to sorted ``(i, j)`` tuples with ``i < j`` and the
    tuple itself is sorted, so two graphs built from the same edge set compare
    equal regardless of input order. ``degrees[v]`` is the number of
    neighbours of ``v``.
    """

    n: int
    edges: tuple[tuple[int, int], ...]
    degrees: tuple[int, ...] = field(compare=False)

    def __post_init__(self) -> None:
        deg = [0] * self.n
        for i, j in self.edges:
            deg[i] += 1
            deg[j] += 1
        if tuple(deg) != tuple(self.degrees):
            raise ValueError("degrees do not match edges; build graphs with MolGraph.from_edges")

    @classmethod
    def from_edges(cls, n: int, edges: Iterable[Sequence[int]]) -> "MolGraph":
        """Build a validated graph (no loops, no multi-edges, indices in range)."""
        if n < 0:
            raise ValueError("n must be non-negative")
        norm = set()
        for e in edges:
            i, j = int(e[0]), int(e[1])
            if i == j:
                raise ValueError(f"self-loop at vertex {i}")
            if not (0 <= i < n and 0 <= j < n):
                raise ValueError(f"edge ({i}, {j}) out of range for n={n}")
            key = (i, j) if i < j else (j, i)
            if key in norm:
                raise ValueError(f"duplicate edge {key}")
            norm.add(key)
        deg = [0] * n
        for i, j in norm:
            deg[i] += 1
            deg[j] += 1
        return cls(n=n, edges=tuple(sorted(norm)), degrees=tuple(deg))

    @property
    def m(self) -> int:
        """Number of edges."""
        return len(self.edges)

    @cached_property
    def adjacency(self) -> tuple[tuple[int, ...], ...]:
        """Sorted neighbour lists."""
        adj: list[list[int]] = [[] for _ in range(self.n)]
        for i, j in self.edges:
            adj[i].append(j)
            adj[j].append(i)
        return tuple(tuple(sorted(a)) for a in adj)

    @cached_property
    def max_degree(self) -> int:
        """Delta(G), the maximum vertex degree (0 for the empty graph)."""
        return max(self.degrees, default=0)

    @cached_property
    def n_components(self) -> int:
        """Number of connected components (0 for the empty graph)."""
        seen = [False] * self.n
        comps = 0
        for s in range(self.n):
            if seen[s]:
                continue
            comps += 1
            seen[s] = True
            queue = deque([s])
            while queue:
                v = queue.popleft()
                for u in self.adjacency[v]:
                    if not seen[u]:
                        seen[u] = True
                        queue.append(u)
        return comps

    @property
    def is_connected(self) -> bool:
        return self.n > 0 and self.n_components == 1

    def relabel(self, perm: Sequence[int]) -> "MolGraph":
        """Return the isomorphic graph with vertex ``v`` renamed ``perm[v]``."""
        if sorted(perm) != list(range(self.n)):
            raise ValueError("perm must be a permutation of range(n)")
        return MolGraph.from_edges(self.n, ((perm[i], perm[j]) for i, j in self.edges))


def graph_from_mol(mol: Chem.Mol) -> MolGraph:
    """Hydrogen-suppressed graph of an RDKit molecule.

    Vertices are the atoms with atomic number > 1, numbered in RDKit atom
    order; every bond between two heavy atoms is an edge. Hydrogens
    (including isotopic and explicit ones) are never vertices.
    """
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    index = {old: new for new, old in enumerate(heavy)}
    edges = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in index and j in index:
            edges.append((index[i], index[j]))
    return MolGraph.from_edges(len(heavy), edges)


# ---------------------------------------------------------------------------
# Standardisation
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Rejection:
    """Why a structure was not turned into a graph."""

    code: str
    detail: str

    def __post_init__(self) -> None:
        if self.code not in REJECTION_CODES:
            raise ValueError(f"unknown rejection code {self.code!r}")


@dataclass(frozen=True)
class StandardizationResult:
    """Outcome of :func:`standardize`.

    Exactly one of ``graph`` and ``rejection`` is set. ``canonical_smiles`` and
    ``inchikey`` describe the standardised parent whenever it could be
    computed (also for most rejections, which helps auditing). ``inchikey`` is
    the full 27-character key, or ``None`` if RDKit's InChI layer failed.
    ``parent_note`` is set when the parent fallback (module docstring,
    step 4) removed fragments that the ChEMBL pipeline kept.
    """

    input_smiles: str
    mol: Optional[Chem.Mol] = field(default=None, compare=False, repr=False)
    graph: Optional[MolGraph] = None
    canonical_smiles: Optional[str] = None
    inchikey: Optional[str] = None
    chembl_exclude_flag: Optional[bool] = None
    parent_note: Optional[str] = None
    rejection: Optional[Rejection] = None

    @property
    def ok(self) -> bool:
        return self.rejection is None


def _inchikey(mol: Chem.Mol) -> Optional[str]:
    try:
        key = Chem.MolToInchiKey(mol)
    except Exception:  # pragma: no cover - RDKit built without InChI
        return None
    return key or None


@lru_cache(maxsize=1)
def _chembl_salts() -> tuple[Chem.Mol, ...]:
    """Salt patterns of the ChEMBL pipeline (its ``data/salts.smi``)."""
    pats = []
    with open(_csp._salts_file) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.strip().split("\t")
            if len(parts) == 2:
                pat = Chem.MolFromSmarts(parts[1])
                if pat is not None:
                    pats.append(pat)
    return tuple(pats)


def _on_chembl_salt_list(frag: Chem.Mol) -> bool:
    """The test of ``get_fragment_parent_mol``: equal size and a substructure match."""
    f = Chem.RemoveHs(frag, sanitize=False)
    f.UpdatePropertyCache(strict=False)
    Chem.SetAromaticity(f)
    return any(
        f.GetNumAtoms() == p.GetNumAtoms()
        and f.GetNumBonds() == p.GetNumBonds()
        and f.HasSubstructMatch(p)
        for p in _chembl_salts()
    )


def _neutral_smiles(frag: Chem.Mol) -> str:
    return Chem.MolToSmiles(_csp.uncharge_mol(Chem.Mol(frag)))


def _n_heavy(frag: Chem.Mol) -> int:
    return sum(1 for a in frag.GetAtoms() if a.GetAtomicNum() > 1)


def _is_simple_counterion(frag: Chem.Mol) -> bool:
    heavy = [a for a in frag.GetAtoms() if a.GetAtomicNum() > 1]
    return len(heavy) == 1 and heavy[0].GetSymbol() in SIMPLE_COUNTERION_ELEMENTS


def _is_organic(frag: Chem.Mol) -> bool:
    """Some carbon bonded to hydrogen or to another carbon (excludes CO3, CN-)."""
    for a in frag.GetAtoms():
        if a.GetAtomicNum() == 6 and (
            a.GetTotalNumHs() > 0 or any(nb.GetAtomicNum() == 6 for nb in a.GetNeighbors())
        ):
            return True
    return False


def _resolve_parent(parent: Chem.Mol) -> tuple[Chem.Mol, Optional[str], Optional[str]]:
    """Parent fallback (module docstring, step 4).

    Returns ``(mol, note, ambiguity)``. ``note`` lists the dropped fragments
    (``None`` if the parent is unchanged); ``ambiguity`` is
    ``"parent_ambiguous"`` or ``"mixture"`` when more than one distinct
    fragment remains (the parent is then returned unchanged), else ``None``.
    """
    frags = list(Chem.GetMolFrags(parent, asMols=True, sanitizeFrags=False))
    if len(frags) <= 1:
        return parent, None, None
    dropped: list[str] = []
    rules = (
        (lambda f: _n_heavy(f) == 0, "no heavy atom"),
        (_is_simple_counterion, "inorganic counter-ion"),
        (lambda f: _neutral_smiles(f) in SUPPLEMENTARY_COUNTERIONS, "supplementary counter-ion"),
    )
    for pred, label in rules:
        hit = [pred(f) for f in frags]
        if any(hit) and not all(hit):
            dropped.extend(f"{Chem.MolToSmiles(f)} ({label})" for f, h in zip(frags, hit) if h)
            frags = [f for f, h in zip(frags, hit) if not h]
    unique: dict[str, Chem.Mol] = {}
    for f in frags:
        unique.setdefault(_neutral_smiles(f), f)
    if len(unique) > 1:
        salts_only = all(_on_chembl_salt_list(f) for f in unique.values())
        return parent, None, "parent_ambiguous" if salts_only else "mixture"
    if not _is_organic(frags[0]):
        # e.g. lithium carbonate or sodium cyanide: removing the cation would
        # turn an inorganic salt into carbonic acid / HCN, so leave it alone
        return parent, None, None
    if len(frags) > 1:
        dropped.append(f"{len(frags) - 1} duplicate cop{'y' if len(frags) == 2 else 'ies'}")
    mol = _csp.uncharge_mol(Chem.Mol(next(iter(unique.values()))))
    return mol, "parent fallback dropped: " + "; ".join(sorted(dropped)), None


def standardize(
    smiles: str,
    *,
    allowed_elements: Iterable[str] = ORGANIC_ELEMENTS,
    min_heavy_atoms: int = 3,
) -> StandardizationResult:
    """Standardise ``smiles`` and build its hydrogen-suppressed graph.

    See the module docstring for the exact pipeline and rejection rules.
    """
    allowed = frozenset(allowed_elements)
    smi = "" if smiles is None else str(smiles).strip()

    def reject(code: str, detail: str, **kw) -> StandardizationResult:
        return StandardizationResult(input_smiles=smi, rejection=Rejection(code, detail), **kw)

    if any(c.isspace() for c in smi):
        return reject(
            "unparsable",
            "SMILES contains whitespace; RDKit would read the text after it as a "
            "molecule name (CXSMILES extensions are not supported)",
        )

    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None or mol.GetNumAtoms() == 0:
            return reject("unparsable", "RDKit could not parse the SMILES" if smi else "empty SMILES")
        try:
            std = _csp.standardize_mol(mol)
            parent, exclude = _csp.get_parent_mol(std)
            parent = Chem.Mol(parent)
            Chem.SanitizeMol(parent)
            parent, note, ambiguity = _resolve_parent(parent)
            if note is not None:
                parent = Chem.Mol(parent)
                Chem.SanitizeMol(parent)
        except Exception as exc:  # the pipeline raises heterogeneous errors
            return reject("standardization_failed", f"{type(exc).__name__}: {exc}")

        can = Chem.MolToSmiles(parent)
        key = _inchikey(parent)
    prov = dict(
        mol=parent,
        canonical_smiles=can,
        inchikey=key,
        chembl_exclude_flag=bool(exclude),
        parent_note=note,
    )

    symbols = sorted({a.GetSymbol() for a in parent.GetAtoms()})
    bad = [s for s in symbols if s not in allowed]
    if bad:
        return reject("disallowed_element", "elements outside the allowed set: " + ",".join(bad), **prov)
    if "C" not in symbols:
        return reject("no_carbon", "parent structure contains no carbon atom", **prov)
    n_frag = len(Chem.GetMolFrags(parent))
    if n_frag > 1:
        if ambiguity == "parent_ambiguous":
            return reject(
                "parent_ambiguous",
                f"{n_frag} fragments remain and each is on the ChEMBL salt list, "
                "so no parent can be chosen",
                **prov,
            )
        return reject("mixture", f"{n_frag} fragments remain after parent extraction", **prov)
    graph = graph_from_mol(parent)
    if graph.n < min_heavy_atoms:
        return reject(
            "too_few_heavy_atoms", f"{graph.n} heavy atoms (< {min_heavy_atoms})", **prov
        )
    if not graph.is_connected:
        return reject(
            "disconnected_graph",
            f"hydrogen-suppressed graph has {graph.n_components} components",
            **prov,
        )
    return StandardizationResult(input_smiles=smi, graph=graph, **prov)
