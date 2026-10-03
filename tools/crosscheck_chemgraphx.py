"""Cross-check TopoQSPR indices against ChemGraphX.

ChemGraphX (Jacob, Clement, Arockiaraj, Jindal, J. Cheminform. 18:39, 2026;
MIT; https://github.com/aditya-jindal/chemgraphx) is an independent
implementation of degree, neighbourhood-degree-sum and distance indices. This
script feeds the *same* hydrogen-suppressed graph to both tools and compares
every index they share.

    # needs a ChemGraphX clone and python-igraph (ChemGraphX's distance code
    # uses it); matplotlib is imported by ChemGraphX's graph.py but never used,
    # so a stub is installed when it is missing.
    PYTHONPATH=src python tools/crosscheck_chemgraphx.py \
        --chemgraphx /path/to/chemgraphx [--write-fixture]

Inputs given to ChemGraphX (each graph through three routes; ChemGraphX must
see a graph isomorphic to TopoQSPR's on every route, and the routes' values
are recorded with their largest relative gap, which comes only from summing
the edges in a different order):

* ``.mol`` V2000 and ``.mol`` V3000, read by ``server/graph.py``
  (``read_mol_v2`` / ``read_mol_v3``), exactly as an uploaded file. For
  molecules the molblock is RDKit's of the hydrogen-free standardised parent
  (heteroatoms and bond orders kept; ChemGraphX ignores both); for the
  closed-form graphs it is an all-carbon single-bond molblock.
* the adjacency-list JSON of the text-area route (``create_graph``), with
  vertices numbered 1..n.

The ChemGraphX classes are called as ``main.py`` calls them
(``ColumnF(file, textarea)`` for degree and distance indices,
``ColumnCD(file, textarea)`` for the degree-sum ones). Values are read from
the instances' ``degree_indices`` / ``distance_indices`` dicts, which the
Flask endpoints only re-label (``capitalize_dictionary``).

``--write-fixture`` writes ``tests/fixtures/chemgraphx/values.json`` (graph
edge lists, the ChemGraphX values and the minimal-example probes), which the
offline test ``tests/test_chemgraphx_crosscheck.py`` compares against
TopoQSPR. Output is deterministic (sorted keys, no timestamps).
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib
import importlib.metadata as md
import io
import json
import math
import sys
import types
import warnings
from pathlib import Path
from typing import Callable, Optional

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "chemgraphx" / "values.json"
DRUGS_JSON = ROOT / "tests" / "fixtures" / "indices" / "drugs.json"

#: ChemGraphX commit the fixture was produced with (clone of 2026-10-02).
CHEMGRAPHX_COMMIT = "d8fea04ffc90ec4acb93748a79392c45c93adc1d"
CHEMGRAPHX_URL = "https://github.com/aditya-jindal/chemgraphx"
#: ChemGraphX files whose code is exercised; their SHA-256 goes in the fixture.
CHEMGRAPHX_FILES = (
    "server/graph.py",
    "server/basic_values.py",
    "server/four_columns.py",
    "server/column_ab.py",
    "server/column_cd.py",
    "server/column_e.py",
    "server/column_f.py",
)

#: Relative tolerance for "agree" (same as the Mordred cross-check).
REL_TOL = 1e-10

# ---------------------------------------------------------------------------
# Index mapping
# ---------------------------------------------------------------------------
#: (ChemGraphX block, ChemGraphX key, TopoQSPR name, relation)
#: block: "degree" (ColumnAB/ColumnF.degree_indices, vertex weight d_u),
#:        "degree_sum" (ColumnCD.degree_indices, vertex weight S_u),
#:        "distance" (ColumnF.distance_indices).
#: relation: "equal", or "W_e - m(m-1)/2" (ChemGraphX edge Wiener uses the
#:           minimum endpoint distance, 0 for adjacent edges; TopoQSPR uses
#:           W(L(G)), whose distances are 1 larger for every pair e != f).
MAPPING: tuple[tuple[str, str, str, str], ...] = (
    ("degree", "first_zagreb", "M1", "equal"),
    ("degree", "second_zagreb", "M2", "equal"),
    ("degree", "randić", "R", "equal"),
    ("degree", "atom_bond_connectivity", "ABC", "equal"),
    ("degree", "harmonic", "H", "equal"),
    ("degree", "sum_connectivity", "SC", "equal"),
    ("degree", "hyper_zagreb", "HM", "equal"),
    ("degree", "geometric_arithmetic", "GA", "equal"),
    ("degree", "sigma", "sigma", "equal"),
    ("degree", "somber", "SO", "equal"),
    ("degree", "forgotten", "F", "equal"),
    ("degree", "symmetric_division_degree", "SDD", "equal"),
    ("degree", "augmented_zagreb", "AZI", "equal"),
    ("degree", "bi_zagreb", "BM", "equal"),
    ("degree", "tri_zagreb", "TM", "equal"),
    ("degree", "geometric_harmonic", "GH", "equal"),
    ("degree", "geometric_bi_zagreb", "GBM", "equal"),
    ("degree", "geometric_tri_zagreb", "GTM", "equal"),
    ("degree", "harmonic_geometric", "HG", "equal"),
    ("degree", "harmonic_bi_zagreb", "HBM", "equal"),
    ("degree", "harmonic_tri_zagreb", "HTM", "equal"),
    ("degree", "bi_zagreb_geometric", "BMG", "equal"),
    ("degree", "bi_zagreb_harmonic", "BMH", "equal"),
    ("degree", "tri_zagreb_geometric", "TMG", "equal"),
    ("degree", "tri_zagreb_harmonic", "TMH", "equal"),
    ("degree_sum", "randić", "S_R", "equal"),
    ("degree_sum", "atom_bond_connectivity", "S_ABC", "equal"),
    ("degree_sum", "harmonic", "S_H", "equal"),
    ("degree_sum", "sum_connectivity", "S_SC", "equal"),
    ("degree_sum", "geometric_arithmetic", "S_GA", "equal"),
    ("degree_sum", "symmetric_division_degree", "S_SDD", "equal"),
    ("degree_sum", "geometric_bi_zagreb", "S_GBM", "equal"),
    ("degree_sum", "tri_zagreb_harmonic", "S_TMH", "equal"),
    ("distance", "wiener_index", "W", "equal"),
    ("distance", "szeged", "Sz", "equal"),
    ("distance", "edge_szeged", "Sz_e", "equal"),
    ("distance", "PI", "PI_v", "equal"),  # ChemGraphX "PI" = n_u + n_v (vertex PI)
    ("distance", "edge_PI", "PI", "equal"),  # ChemGraphX "edge PI" = m_u + m_v (Khadikar PI)
    ("distance", "mostar", "Mo", "equal"),
    ("distance", "edge_mostar", "Mo_e", "equal"),
    ("distance", "edge_wiener", "W_e", "W_e - m(m-1)/2"),
)

RELATIONS: dict[str, Callable[[float, int], float]] = {
    "equal": lambda t, m: t,
    "W_e - m(m-1)/2": lambda t, m: t - m * (m - 1) // 2,
}

# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------
#: Benzenoid PAHs (SMILES).
PAHS: tuple[tuple[str, str], ...] = (
    ("benzene", "c1ccccc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
    ("phenanthrene", "c1ccc2c(c1)ccc1ccccc12"),
    ("tetracene", "c1ccc2cc3cc4ccccc4cc3cc2c1"),
    ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34"),
    ("chrysene", "c1ccc2c(c1)ccc1c3ccccc3ccc21"),
    ("triphenylene", "c1ccc2c(c1)c1ccccc1c1ccccc21"),
    ("benzo[a]pyrene", "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34"),
    ("coronene", "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61"),
)


def _path(n: int) -> list[tuple[int, int]]:
    return [(i, i + 1) for i in range(n - 1)]


def _cycle(n: int) -> list[tuple[int, int]]:
    return [(i, (i + 1) % n) for i in range(n)]


#: Closed-form graphs: name -> (n, edges).
CLOSED_FORM: dict[str, tuple[int, list[tuple[int, int]]]] = {
    "P3": (3, _path(3)),
    "P5": (5, _path(5)),
    "P8": (8, _path(8)),
    "C5": (5, _cycle(5)),
    "C6": (6, _cycle(6)),
    "C9": (9, _cycle(9)),
    "K1,4": (5, [(0, i) for i in range(1, 5)]),
    "K4": (4, [(i, j) for i in range(4) for j in range(i + 1, 4)]),
    "K2,3": (5, [(i, j) for i in range(2) for j in range(2, 5)]),
}


def panel() -> list[dict]:
    """The fixed panel: closed-form graphs, PAHs, the 24 drugs of drugs.json.

    Each entry: name, kind, smiles (None for closed-form graphs), graph
    (a :class:`topoqspr.graph.MolGraph`) and mol (RDKit hydrogen-free parent,
    None for closed-form graphs).
    """
    from rdkit import Chem

    from topoqspr.graph import MolGraph, graph_from_mol, standardize

    out = []
    for name, (n, edges) in CLOSED_FORM.items():
        out.append(dict(name=name, kind="closed_form", smiles=None,
                        graph=MolGraph.from_edges(n, edges), mol=None))
    drugs = json.loads(DRUGS_JSON.read_text())["drugs"]
    mols = [(n, s, "pah") for n, s in PAHS] + [(d["name"], d["smiles"], "drug") for d in drugs]
    for name, smi, kind in mols:
        res = standardize(smi)
        if not res.ok:
            raise RuntimeError(f"{name}: {res.rejection}")
        mol = Chem.RemoveHs(res.mol)
        g = graph_from_mol(mol)
        if g != res.graph:  # RemoveHs must not change the heavy-atom graph
            raise RuntimeError(f"{name}: RemoveHs changed the graph")
        out.append(dict(name=name, kind=kind, smiles=smi, graph=res.graph, mol=mol))
    return out


# ---------------------------------------------------------------------------
# ChemGraphX driver
# ---------------------------------------------------------------------------
class _Upload(io.BytesIO):
    """Stands in for werkzeug's FileStorage: bytes plus a ``filename``."""

    def __init__(self, data: bytes, filename: str):
        super().__init__(data)
        self.filename = filename


def clone_commit(clone: Path) -> Optional[str]:
    """HEAD commit of a git clone, read from ``.git`` without running git."""
    gitdir = clone / ".git"
    try:
        head = (gitdir / "HEAD").read_text().strip()
    except OSError:
        return None
    if not head.startswith("ref: "):
        return head
    ref = head[5:]
    ref_file = gitdir / ref
    if ref_file.exists():
        return ref_file.read_text().strip()
    packed = gitdir / "packed-refs"
    if packed.exists():
        for line in packed.read_text().splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1] == ref:
                return parts[0]
    return None


def file_hashes(clone: Path) -> dict[str, str]:
    return {f: hashlib.sha256((clone / f).read_bytes()).hexdigest() for f in CHEMGRAPHX_FILES}


def load_chemgraphx(clone: Path) -> types.SimpleNamespace:
    """Import ChemGraphX's server modules from a clone."""
    server = clone / "server"
    if not (server / "column_f.py").exists():
        raise FileNotFoundError(f"{server} is not a ChemGraphX server directory")
    try:
        importlib.import_module("matplotlib.pyplot")
    except ImportError:  # graph.py imports pyplot only for an unused plot method
        mpl = types.ModuleType("matplotlib")
        plt = types.ModuleType("matplotlib.pyplot")
        mpl.pyplot = plt  # type: ignore[attr-defined]
        sys.modules.setdefault("matplotlib", mpl)
        sys.modules.setdefault("matplotlib.pyplot", plt)
    if str(server) not in sys.path:
        sys.path.insert(0, str(server))
    column_f = importlib.import_module("column_f")
    column_cd = importlib.import_module("column_cd")
    return types.SimpleNamespace(ColumnF=column_f.ColumnF, ColumnCD=column_cd.ColumnCD)


def adjacency_json(n: int, edges) -> str:
    adj: dict[str, list[int]] = {str(i + 1): [] for i in range(n)}
    for u, v in edges:
        adj[str(u + 1)].append(v + 1)
        adj[str(v + 1)].append(u + 1)
    return json.dumps(adj)


def molblock(entry: dict, v3000: bool) -> bytes:
    from rdkit import Chem

    mol = entry["mol"]
    if mol is None:  # closed-form graph: all-carbon, single bonds
        g = entry["graph"]
        rw = Chem.RWMol()
        for _ in range(g.n):
            rw.AddAtom(Chem.Atom(6))
        for u, v in g.edges:
            rw.AddBond(u, v, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        return Chem.MolToMolBlock(mol, kekulize=False, forceV3000=v3000).encode()
    return Chem.MolToMolBlock(mol, forceV3000=v3000).encode()


def _num(x) -> float:
    return float(x)


def run_chemgraphx(cgx, source) -> tuple[dict[str, dict[str, float]], object]:
    """Run ChemGraphX on one input; ``source`` is ``("mol", bytes, name)`` or
    ``("adjacency", json_str)``. Returns (values by block, ChemGraphX nx graph)."""

    def args():
        if source[0] == "mol":
            return _Upload(source[1], source[2]), None
        return "", source[1]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = cgx.ColumnF(*args())
        cd = cgx.ColumnCD(*args())
    vals = {
        "degree": {k: _num(v) for k, v in f.degree_indices.items()},
        "degree_sum": {k: _num(v) for k, v in cd.degree_indices.items()},
        "distance": {k: _num(v) for k, v in f.distance_indices.items()},
    }
    return vals, f.graph


def chemgraphx_values(cgx, entry: dict) -> dict:
    """ChemGraphX values for one panel entry via all three input routes."""
    import networkx as nx

    g = entry["graph"]
    ours = nx.Graph()
    ours.add_nodes_from(range(g.n))
    ours.add_edges_from(g.edges)
    routes = {
        "mol_v2000": ("mol", molblock(entry, False), "panel.mol"),
        "mol_v3000": ("mol", molblock(entry, True), "panel.mol"),
        "adjacency": ("adjacency", adjacency_json(g.n, g.edges)),
    }
    results = {}
    for route, src in routes.items():
        vals, G = run_chemgraphx(cgx, src)
        if not nx.vf2pp_is_isomorphic(G, ours):
            raise RuntimeError(f"{entry['name']}: ChemGraphX graph via {route} is not our graph")
        results[route] = vals
    ref = results["mol_v2000"]
    # The routes give ChemGraphX the edges in different orders, so floating
    # sums may differ in the last bit; report the largest relative gap.
    gap = max(
        rel_diff(results[r][b][k], ref[b][k]) for r in results for b in ref for k in ref[b]
    )
    return {"values": ref, "route_max_rel_diff": gap, "routes": sorted(routes)}


# ---------------------------------------------------------------------------
# Minimal examples for every discrepancy class
# ---------------------------------------------------------------------------
def probes(cgx) -> dict:
    """ChemGraphX output on the minimal reproducing inputs of the report."""
    from rdkit import Chem

    def adj(n, edges, keys):
        vals, _ = run_chemgraphx(cgx, ("adjacency", adjacency_json(n, edges)))
        flat = {f"{b}.{k}": v for b, d in vals.items() for k, v in d.items()}
        return {k: flat[k] for k in keys}

    def mol(smiles, keys, add_h=False, on_atoms=None):
        m = Chem.MolFromSmiles(smiles)
        if add_h:
            m = Chem.AddHs(m, onlyOnAtoms=on_atoms) if on_atoms else Chem.AddHs(m)
        vals, G = run_chemgraphx(cgx, ("mol", Chem.MolToMolBlock(m).encode(), "probe.mol"))
        flat = {f"{b}.{k}": v for b, d in vals.items() for k, v in d.items()}
        out = {k: flat[k] for k in keys}
        out["n_vertices"] = G.number_of_nodes()
        out["n_edges"] = G.number_of_edges()
        return out

    dist = ["distance.wiener_index", "distance.szeged", "distance.PI", "distance.edge_PI",
            "distance.mostar", "distance.edge_wiener", "distance.vertex_edge_wiener"]
    out = {
        "P3": adj(3, _path(3), dist),
        "P4": adj(4, _path(4), dist),
        "K2": adj(2, _path(2), ["degree.augmented_zagreb", "distance.wiener_index"]),
        "2K2": adj(4, [(0, 1), (2, 3)], dist),
        "P3+K2": adj(5, [(0, 1), (1, 2), (3, 4)], dist),
        "acetic_acid.mol": mol("CC(=O)O", ["degree.first_zagreb"] + dist),
        "sodium_acetate_Na_last.mol": mol("CC(=O)[O-].[Na+]", ["degree.first_zagreb"] + dist),
        "sodium_acetate_Na_first.mol": mol("[Na+].CC(=O)[O-]", ["degree.first_zagreb"] + dist),
        "toluene.mol": mol("Cc1ccccc1", ["degree.first_zagreb"] + dist),
        "toluene_explicit_H.mol": mol("Cc1ccccc1", ["degree.first_zagreb"] + dist, add_h=True),
        "alanine.mol": mol("C[C@H](N)C(=O)O", ["degree.first_zagreb", "distance.wiener_index"]),
        "alanine_one_stereo_H.mol": mol("C[C@H](N)C(=O)O",
                                        ["degree.first_zagreb", "distance.wiener_index"],
                                        add_h=True, on_atoms=(1,)),
    }
    for n in (300, 600, 1000):
        out[f"P{n}"] = adj(n, _path(n), ["distance.wiener_index", "distance.szeged",
                                         "distance.edge_wiener"])
    return out


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------
def rel_diff(a: float, b: float) -> float:
    """|a - b| / |b| (absolute difference when b == 0); NaN if either is non-finite."""
    if not (math.isfinite(a) and math.isfinite(b)):
        return float("nan")
    return abs(a - b) / abs(b) if b != 0 else abs(a - b)


def compare(fixture: dict) -> dict:
    """Compare recorded ChemGraphX values with TopoQSPR, per mapped index.

    Returns {topoqspr_name: {"n": molecules compared, "max_rel": float,
    "worst": molecule name, "failures": [names above REL_TOL]}}.
    """
    from topoqspr.graph import MolGraph
    from topoqspr.indices import compute_indices

    names = sorted({t for _, _, t, _ in MAPPING})
    table = {t: {"n": 0, "max_rel": 0.0, "worst": None, "failures": []} for t in names}
    for mol in fixture["panel"]:
        g = MolGraph.from_edges(mol["n"], mol["edges"])
        topo, _ = compute_indices(g, names)
        for block, key, tname, rel in MAPPING:
            c = mol["chemgraphx"][block][key]
            t = RELATIONS[rel](topo[tname], g.m)
            d = rel_diff(float(c), float(t))
            row = table[tname]
            row["n"] += 1
            if not d <= REL_TOL:
                row["failures"].append(mol["name"])
            if d > row["max_rel"] or row["worst"] is None:
                row["max_rel"], row["worst"] = d, mol["name"]
    return table


def build_fixture(clone: Path) -> dict:
    cgx = load_chemgraphx(clone)
    commit = clone_commit(clone)
    entries = []
    for e in panel():
        res = chemgraphx_values(cgx, e)
        g = e["graph"]
        entries.append({
            "name": e["name"],
            "kind": e["kind"],
            "smiles": e["smiles"],
            "n": g.n,
            "edges": [list(x) for x in g.edges],
            "chemgraphx": res["values"],
            "chemgraphx_routes": res["routes"],
            "chemgraphx_route_max_rel_diff": res["route_max_rel_diff"],
        })
    versions = {}
    for pkg in ("igraph", "networkx", "numpy", "rdkit", "chembl_structure_pipeline"):
        try:
            versions[pkg] = md.version(pkg)
        except md.PackageNotFoundError:
            versions[pkg] = None
    return {
        "chemgraphx": {
            "url": CHEMGRAPHX_URL,
            "commit": commit,
            "file_sha256": file_hashes(clone),
        },
        "generated_with": versions,
        "rel_tol": REL_TOL,
        "mapping": [list(m) for m in MAPPING],
        "panel": entries,
        "probes": probes(cgx),
    }


def _dumps(obj) -> str:
    def fix(x):
        if isinstance(x, float) and not math.isfinite(x):
            return {"nonfinite": repr(x)}  # JSON has no inf/nan
        if isinstance(x, dict):
            return {k: fix(v) for k, v in x.items()}
        if isinstance(x, list):
            return [fix(v) for v in x]
        return x

    return json.dumps(fix(obj), indent=1, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"


def report(fixture: dict) -> str:
    table = compare(fixture)
    by_name = {t: (b, k, r) for b, k, t, r in MAPPING}
    lines = ["| ChemGraphX (block.key) | TopoQSPR | relation | graphs | max rel. diff | status |",
             "|---|---|---|---|---|---|"]
    for _, _, t, _ in MAPPING:
        b, k, r = by_name[t]
        row = table[t]
        status = "agree" if not row["failures"] else f"DIFFER ({len(row['failures'])})"
        lines.append(f"| `{b}.{k}` | `{t}` | {r} | {row['n']} | {row['max_rel']:.1e} | {status} |")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--chemgraphx", type=Path, required=True, help="path to a ChemGraphX clone")
    ap.add_argument("--write-fixture", action="store_true", help=f"write {FIXTURE.relative_to(ROOT)}")
    args = ap.parse_args(argv)
    commit = clone_commit(args.chemgraphx)
    if commit != CHEMGRAPHX_COMMIT:
        print(f"warning: clone is at {commit}, fixture pin is {CHEMGRAPHX_COMMIT}", file=sys.stderr)
    fx = build_fixture(args.chemgraphx)
    gap = max(e["chemgraphx_route_max_rel_diff"] for e in fx["panel"])
    print(f"panel: {len(fx['panel'])} graphs; largest relative gap between the three "
          f"ChemGraphX input routes: {gap:.1e}")
    print(report(fx))
    print(json.dumps(fx["probes"], indent=1, sort_keys=True, default=str))
    if args.write_fixture:
        FIXTURE.parent.mkdir(parents=True, exist_ok=True)
        FIXTURE.write_text(_dumps(fx))
        print(f"wrote {FIXTURE}")
    return 0


if __name__ == "__main__":
    with contextlib.suppress(BrokenPipeError):
        sys.exit(main())
