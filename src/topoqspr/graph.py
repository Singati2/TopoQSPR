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
1. Parse with RDKit (``Chem.MolFromSmiles``, default sanitisation).
2. ChEMBL Structure Pipeline: ``standardize_mol`` then ``get_parent_mol``
   (Bento et al., J. Cheminform. 12 (2020) 51). This removes counter-ions,
   salts and solvents listed in the pipeline's dictionaries, neutralises
   charges where possible, and strips isotope labels (``[2H]`` becomes H).
3. Reject the result if

   * it contains an element outside ``allowed_elements``
     (default :data:`ORGANIC_ELEMENTS`)  -> ``disallowed_element``
   * it contains no carbon                -> ``no_carbon``
   * more than one fragment remains        -> ``mixture``
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
from functools import cached_property
from typing import Iterable, Optional, Sequence

from chembl_structure_pipeline import standardizer as _csp
from rdkit import Chem, rdBase

__all__ = [
    "ORGANIC_ELEMENTS",
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

#: Every rejection code :func:`standardize` can emit.
REJECTION_CODES: tuple[str, ...] = (
    "unparsable",
    "standardization_failed",
    "disallowed_element",
    "no_carbon",
    "mixture",
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
    """

    input_smiles: str
    mol: Optional[Chem.Mol] = field(default=None, compare=False, repr=False)
    graph: Optional[MolGraph] = None
    canonical_smiles: Optional[str] = None
    inchikey: Optional[str] = None
    chembl_exclude_flag: Optional[bool] = None
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

    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None or mol.GetNumAtoms() == 0:
            return reject("unparsable", "RDKit could not parse the SMILES" if smi else "empty SMILES")
        try:
            std = _csp.standardize_mol(mol)
            parent, exclude = _csp.get_parent_mol(std)
            parent = Chem.Mol(parent)
            Chem.SanitizeMol(parent)
        except Exception as exc:  # the pipeline raises heterogeneous errors
            return reject("standardization_failed", f"{type(exc).__name__}: {exc}")

        can = Chem.MolToSmiles(parent)
        key = _inchikey(parent)
    prov = dict(mol=parent, canonical_smiles=can, inchikey=key, chembl_exclude_flag=bool(exclude))

    symbols = sorted({a.GetSymbol() for a in parent.GetAtoms()})
    bad = [s for s in symbols if s not in allowed]
    if bad:
        return reject("disallowed_element", "elements outside the allowed set: " + ",".join(bad), **prov)
    if "C" not in symbols:
        return reject("no_carbon", "parent structure contains no carbon atom", **prov)
    n_frag = len(Chem.GetMolFrags(parent))
    if n_frag > 1:
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
