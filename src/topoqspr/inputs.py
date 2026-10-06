"""Structure and graph inputs beyond SMILES: MOL/SDF, PDB and adjacency lists/matrices.

Molecular formats (SMILES, MOL/SDF, PDB) are converted to SMILES with RDKit and then go
through :func:`topoqspr.graph.standardize`, exactly like the study (parent structure,
hydrogen-suppressed graph, explicit rejection reasons). Adjacency input describes an
abstract graph (no elements, no bond orders): it is turned into a :class:`MolGraph`
directly, so graph-theoretic indices are available for any simple connected graph, as in
ChemGraphX, while chemistry-specific steps do not apply.

Every record yields an :class:`InputRecord` with either a graph or a rejection reason.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from typing import Iterable, Optional

from rdkit import Chem, RDLogger

from .graph import MolGraph, standardize

__all__ = ["InputRecord", "parse_input", "FORMATS", "MAX_RECORDS", "MAX_VERTICES"]

FORMATS = ("smiles", "sdf", "mol", "pdb", "adjacency_matrix", "adjacency_list")
MAX_RECORDS = 300
MAX_VERTICES = 500

RDLogger.DisableLog("rdApp.*")


@dataclass
class InputRecord:
    label: str
    kind: str  # "molecule" or "graph"
    graph: Optional[MolGraph] = None
    smiles: Optional[str] = None
    reason: Optional[str] = None
    inchikey: Optional[str] = None
    input_text: Optional[str] = None


def _molecule(label: str, mol) -> InputRecord:
    if mol is None:
        return InputRecord(label, "molecule", reason="unparsable structure")
    smi = Chem.MolToSmiles(mol)
    res = standardize(smi)
    graph = getattr(res, "graph", None)
    if graph is None:
        rej = getattr(res, "rejection", None)
        return InputRecord(label, "molecule", smiles=smi,
                           reason=str(getattr(rej, "code", None) or getattr(res, "code", None) or res))
    return InputRecord(label, "molecule", graph=graph, smiles=getattr(res, "canonical_smiles", smi),
                       inchikey=getattr(res, "inchikey", None))


def _graph(label: str, n: int, edges: Iterable[tuple[int, int]]) -> InputRecord:
    if n > MAX_VERTICES:
        return InputRecord(label, "graph", reason=f"more than {MAX_VERTICES} vertices")
    try:
        g = MolGraph.from_edges(n, edges)
    except ValueError as exc:
        return InputRecord(label, "graph", reason=f"invalid graph: {exc}")
    if n == 0 or not g.edges:
        return InputRecord(label, "graph", reason="graph has no edges")
    if not _connected(g):
        return InputRecord(label, "graph", reason="graph is disconnected")
    return InputRecord(label, "graph", graph=g)


def _connected(g: MolGraph) -> bool:
    adj = {i: set() for i in range(g.n)}
    for u, v in g.edges:
        adj[u].add(v)
        adj[v].add(u)
    seen, stack = {0}, [0]
    while stack:
        for w in adj[stack.pop()] - seen:
            seen.add(w)
            stack.append(w)
    return len(seen) == g.n


def _parse_matrix(text: str) -> list[list[int]]:
    """JSON ``[[0,1],[1,0]]`` (optionally after ``name =``) or whitespace/comma rows."""
    body = text.split("=", 1)[1] if "=" in text and "[" in text.split("=", 1)[1] else text
    body = body.strip()
    if body.startswith("["):
        rows = json.loads(body)
    else:
        rows = [[int(x) for x in re.split(r"[\s,;]+", line.strip()) if x != ""]
                for line in body.splitlines() if line.strip()]
    return [[int(x) for x in r] for r in rows]


def _matrix_graph(label: str, text: str) -> InputRecord:
    try:
        m = _parse_matrix(text)
    except (ValueError, json.JSONDecodeError) as exc:
        return InputRecord(label, "graph", reason=f"unparsable adjacency matrix: {exc}")
    n = len(m)
    if any(len(r) != n for r in m):
        return InputRecord(label, "graph", reason="adjacency matrix is not square")
    if any(m[i][j] not in (0, 1) for i in range(n) for j in range(n)):
        return InputRecord(label, "graph", reason="adjacency matrix must contain only 0 and 1")
    if any(m[i][i] for i in range(n)):
        return InputRecord(label, "graph", reason="adjacency matrix has self-loops")
    if any(m[i][j] != m[j][i] for i in range(n) for j in range(n)):
        return InputRecord(label, "graph", reason="adjacency matrix is not symmetric")
    return _graph(label, n, [(i, j) for i in range(n) for j in range(i + 1, n) if m[i][j]])


def _list_graph(label: str, text: str) -> InputRecord:
    """JSON ``{"1": [2, 3], ...}`` (ChemGraphX's 1-based adjacency list), lines ``1: 2 3``
    (a vertex and its neighbours) or lines ``1 2`` (one edge per line)."""
    try:
        body = text.strip()
        if body.startswith("{"):
            adj = {int(k): [int(x) for x in v] for k, v in json.loads(body).items()}
        else:
            adj = {}
            for line in body.splitlines():
                if not line.strip():
                    continue
                if ":" in line:
                    head, _, tail = line.partition(":")
                    adj.setdefault(int(head), []).extend(int(x) for x in re.split(r"[\s,;]+", tail.strip()) if x)
                else:
                    pair = [int(x) for x in re.split(r"[\s,;]+", line.strip()) if x]
                    if len(pair) != 2:
                        raise ValueError(f"line {line.strip()!r} is neither 'u: v w ...' nor an edge 'u v'")
                    adj.setdefault(pair[0], []).append(pair[1])
    except (ValueError, json.JSONDecodeError) as exc:
        return InputRecord(label, "graph", reason="unparsable adjacency list (accepted: JSON {\"1\": [2, 3]}, "
                           f"lines '1: 2 3', or one edge '1 2' per line): {exc}")
    verts = sorted(set(adj) | {w for vs in adj.values() for w in vs})
    index = {v: i for i, v in enumerate(verts)}
    edges = set()
    for u, vs in adj.items():
        for v in vs:
            if u == v:
                return InputRecord(label, "graph", reason="adjacency list has self-loops")
            edges.add(tuple(sorted((index[u], index[v]))))
    return _graph(label, len(verts), sorted(edges))


def parse_input(text: str, fmt: str, label: str = "input") -> list[InputRecord]:
    """Parse ``text`` in format ``fmt`` (one of :data:`FORMATS`) into records."""
    if fmt not in FORMATS:
        raise ValueError(f"unknown format {fmt!r}; allowed: {FORMATS}")
    if fmt == "smiles":
        out = []
        for i, line in enumerate(l for l in text.splitlines() if l.strip()):
            parts = line.split()
            rec = _molecule(parts[1] if len(parts) > 1 else f"{label}_{i + 1}",
                            Chem.MolFromSmiles(parts[0]))
            rec.input_text = parts[0]
            out.append(rec)
        records = out
    elif fmt in ("sdf", "mol"):
        blocks = [b for b in re.split(r"^\$\$\$\$\s*$", text, flags=re.M) if b.strip()]
        records = []
        for i, block in enumerate(blocks):
            mol = Chem.MolFromMolBlock(block, removeHs=True)
            name = (mol.GetProp("_Name").strip() if mol is not None and mol.HasProp("_Name") else "")
            records.append(_molecule(name or f"{label}_{i + 1}", mol))
    elif fmt == "pdb":
        records = [_molecule(label, Chem.MolFromPDBBlock(text, removeHs=True))]
    elif fmt == "adjacency_matrix":
        records = [_matrix_graph(label, text)]
    else:
        records = [_list_graph(label, text)]
    if len(records) > MAX_RECORDS:
        raise ValueError(f"more than {MAX_RECORDS} records")
    return records
