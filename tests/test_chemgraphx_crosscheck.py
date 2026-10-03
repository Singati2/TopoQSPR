"""Cross-check of TopoQSPR indices against ChemGraphX (docs/CHEMGRAPHX_CROSSCHECK.md).

Offline tests compare TopoQSPR with the ChemGraphX values recorded in
``tests/fixtures/chemgraphx/values.json`` by ``tools/crosscheck_chemgraphx.py``
and check, with TopoQSPR alone, the minimal examples behind every documented
difference. ``test_fixture_is_reproducible_with_chemgraphx`` re-runs
ChemGraphX; it needs a clone (``CHEMGRAPHX_DIR``) and python-igraph and is
skipped otherwise.
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import math
import os
from pathlib import Path

import networkx as nx
import pytest

from topoqspr.graph import MolGraph, standardize
from topoqspr.indices import INDEX_NAMES, compute_indices

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "chemgraphx" / "values.json"


def _load_tool():
    spec = importlib.util.spec_from_file_location(
        "crosscheck_chemgraphx", ROOT / "tools" / "crosscheck_chemgraphx.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


TOOL = _load_tool()


def _decode(x):
    """Undo the tool's JSON encoding of non-finite floats."""
    if isinstance(x, dict):
        if set(x) == {"nonfinite"}:
            return float(x["nonfinite"])
        return {k: _decode(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_decode(v) for v in x]
    return x


FX = _decode(json.loads(FIXTURE.read_text()))
PROBES = FX["probes"]
PANEL = {e["name"]: e for e in FX["panel"]}
TOPO_NAMES = [t for _, _, t, _ in TOOL.MAPPING]


def _graph(entry) -> MolGraph:
    return MolGraph.from_edges(entry["n"], entry["edges"])


def _path(n: int) -> MolGraph:
    return MolGraph.from_edges(n, [(i, i + 1) for i in range(n - 1)])


# ---------------------------------------------------------------------------
# fixture integrity
# ---------------------------------------------------------------------------
def test_fixture_pins_commit_files_and_mapping():
    assert FX["chemgraphx"]["commit"] == TOOL.CHEMGRAPHX_COMMIT
    assert sorted(FX["chemgraphx"]["file_sha256"]) == sorted(TOOL.CHEMGRAPHX_FILES)
    assert [tuple(m) for m in FX["mapping"]] == list(TOOL.MAPPING)
    assert FX["rel_tol"] == TOOL.REL_TOL


def test_mapping_targets_exist_and_are_unique():
    assert len(TOPO_NAMES) == len(set(TOPO_NAMES)) == 41
    assert set(TOPO_NAMES) <= set(INDEX_NAMES)
    assert set(TOOL.RELATIONS) == {r for *_, r in TOOL.MAPPING}
    for e in FX["panel"]:
        for block, key, _, _ in TOOL.MAPPING:
            assert key in e["chemgraphx"][block], (e["name"], block, key)


def test_panel_composition():
    kinds = [e["kind"] for e in FX["panel"]]
    assert len(kinds) >= 30
    assert kinds.count("closed_form") == len(TOOL.CLOSED_FORM)
    assert kinds.count("pah") == len(TOOL.PAHS)
    drugs = json.loads(TOOL.DRUGS_JSON.read_text())["drugs"]
    assert kinds.count("drug") == len(drugs) == 24
    assert [e["name"] for e in FX["panel"]] == (
        list(TOOL.CLOSED_FORM) + [n for n, _ in TOOL.PAHS] + [d["name"] for d in drugs]
    )


def test_fixture_graphs_are_the_topoqspr_graphs():
    """The recorded edge lists are what TopoQSPR builds today."""
    for e in FX["panel"]:
        if e["kind"] == "closed_form":
            n, edges = TOOL.CLOSED_FORM[e["name"]]
            expected = MolGraph.from_edges(n, edges)
        else:
            res = standardize(e["smiles"])
            assert res.ok, (e["name"], res.rejection)
            expected = res.graph
        assert _graph(e) == expected, e["name"]


def test_chemgraphx_input_routes_agree():
    """.mol V2000, .mol V3000 and the adjacency JSON gave the same values."""
    for e in FX["panel"]:
        assert e["chemgraphx_routes"] == ["adjacency", "mol_v2000", "mol_v3000"]
        assert e["chemgraphx_route_max_rel_diff"] <= 1e-12, e["name"]


# ---------------------------------------------------------------------------
# the cross-check itself
# ---------------------------------------------------------------------------
_TABLE = TOOL.compare(FX)


@pytest.mark.parametrize("name", TOPO_NAMES)
def test_overlapping_index_agrees_with_chemgraphx(name):
    row = _TABLE[name]
    assert row["n"] == len(FX["panel"])
    assert row["failures"] == [], f"{name} differs from ChemGraphX on {row['failures']}"
    assert row["max_rel"] <= TOOL.REL_TOL


def test_relative_difference_helper():
    assert TOOL.rel_diff(1.0, 1.0) == 0.0
    assert TOOL.rel_diff(0.0, 0.0) == 0.0
    assert TOOL.rel_diff(2.0, 1.0) == 1.0
    assert math.isnan(TOOL.rel_diff(float("inf"), 1.0))
    assert math.isnan(TOOL.rel_diff(1.0, float("nan")))


# ---------------------------------------------------------------------------
# minimal examples of every documented difference (TopoQSPR side computed)
# ---------------------------------------------------------------------------
def test_pi_naming_is_swapped_on_p3():
    """ChemGraphX 'PI' is the vertex PI (TopoQSPR PI_v); its 'edge PI' is
    Khadikar's PI (TopoQSPR PI)."""
    v, _ = compute_indices(_path(3), ["PI", "PI_v"])
    assert (v["PI"], v["PI_v"]) == (2, 6)
    assert PROBES["P3"]["distance.edge_PI"] == v["PI"]
    assert PROBES["P3"]["distance.PI"] == v["PI_v"]


@pytest.mark.parametrize("n", [3, 4])
def test_edge_wiener_convention_on_paths(n):
    """ChemGraphX edge Wiener = W(L(G)) - C(m, 2) (minimum endpoint distance)."""
    g = _path(n)
    v, _ = compute_indices(g, ["W_e"])
    assert v["W_e"] == (n - 2) * (n - 1) * n // 6  # L(P_n) = P_{n-1}
    assert PROBES[f"P{n}"]["distance.edge_wiener"] == v["W_e"] - math.comb(g.m, 2)


def test_edge_wiener_offset_identity_on_the_panel():
    """W_e - C(m,2) equals the sum over edge pairs of the minimum endpoint
    distance, computed here independently with networkx."""
    for e in FX["panel"]:
        g = _graph(e)
        G = nx.Graph(g.edges)
        D = dict(nx.all_pairs_shortest_path_length(G))
        d0 = sum(
            min(D[a][x], D[a][y], D[b][x], D[b][y])
            for (a, b), (x, y) in itertools.combinations(g.edges, 2)
        )
        v, _ = compute_indices(g, ["W_e"])
        assert v["W_e"] - math.comb(g.m, 2) == d0, e["name"]


def test_k2_azi_is_nan_here_and_infinite_in_chemgraphx():
    v, r = compute_indices(_path(2), ["AZI", "W"])
    assert math.isnan(v["AZI"]) and "division by zero" in r["AZI"]
    assert v["W"] == 1
    assert PROBES["K2"]["degree.augmented_zagreb"] == math.inf


@pytest.mark.parametrize(
    "name, n, edges", [("2K2", 4, [(0, 1), (2, 3)]), ("P3+K2", 5, [(0, 1), (1, 2), (3, 4)])]
)
def test_disconnected_graph_distance_indices(name, n, edges):
    """TopoQSPR: NaN with a reason. ChemGraphX: W infinite, Sz/PI/Mo finite."""
    v, r = compute_indices(MolGraph.from_edges(n, edges), ["W", "Sz", "PI_v", "Mo", "W_e"])
    for k in v:
        assert math.isnan(v[k]) and "components" in r[k]
    p = PROBES[name]
    assert p["distance.wiener_index"] == math.inf
    assert math.isfinite(p["distance.szeged"]) and math.isfinite(p["distance.PI"])


def test_salt_counter_ion_position_changes_chemgraphx_wiener():
    """An isolated Na in a .mol file: ChemGraphX W is finite when Na is the
    last atom and infinite when it is the first; TopoQSPR uses the parent."""
    res = standardize("CC(=O)[O-].[Na+]")
    assert res.ok and res.graph.n == 4
    v, _ = compute_indices(res.graph, ["W", "M1"])
    last, first = PROBES["sodium_acetate_Na_last.mol"], PROBES["sodium_acetate_Na_first.mol"]
    assert last["n_vertices"] == first["n_vertices"] == 5
    assert last["distance.wiener_index"] == v["W"] == 9
    assert first["distance.wiener_index"] == math.inf
    assert last["degree.first_zagreb"] == first["degree.first_zagreb"] == v["M1"]
    assert PROBES["acetic_acid.mol"]["distance.wiener_index"] == 9


def test_explicit_hydrogens_become_chemgraphx_vertices():
    """ChemGraphX keeps every atom of the .mol file; TopoQSPR never makes H a vertex."""
    for smi in ("Cc1ccccc1", "[H]C([H])([H])c1ccccc1"):
        res = standardize(smi)
        assert res.ok and res.graph.n == 7
    g = standardize("Cc1ccccc1").graph
    v, _ = compute_indices(g, ["M1", "W"])
    assert PROBES["toluene.mol"]["n_vertices"] == 7
    assert PROBES["toluene.mol"]["degree.first_zagreb"] == v["M1"] == 30
    assert PROBES["toluene.mol"]["distance.wiener_index"] == v["W"] == 42
    assert PROBES["toluene_explicit_H.mol"]["n_vertices"] == 15
    assert PROBES["toluene_explicit_H.mol"]["degree.first_zagreb"] == 78
    # a single stereo hydrogen drawn explicitly changes the ChemGraphX graph
    a = standardize("C[C@H](N)C(=O)O").graph
    va, _ = compute_indices(a, ["M1", "W"])
    assert PROBES["alanine.mol"]["degree.first_zagreb"] == va["M1"] == 22
    assert PROBES["alanine_one_stereo_H.mol"]["n_vertices"] == 7
    assert PROBES["alanine_one_stereo_H.mol"]["degree.first_zagreb"] == 30


def test_float32_accumulation_on_long_paths():
    """ChemGraphX sums distances in float32: exact for P300 (W < 2**24), not
    for P600 / P1000. TopoQSPR is exact (W(P_n) = C(n+1, 3))."""
    for n in (300, 600, 1000):
        exact = math.comb(n + 1, 3)
        v, _ = compute_indices(_path(n), ["W"])
        assert v["W"] == exact
        rec = PROBES[f"P{n}"]["distance.wiener_index"]
        if exact < 2**24:
            assert rec == exact
        else:
            assert rec != exact and abs(rec - exact) / exact < 1e-5


# ---------------------------------------------------------------------------
# optional: re-run ChemGraphX
# ---------------------------------------------------------------------------
@pytest.mark.network
def test_fixture_is_reproducible_with_chemgraphx():
    clone = os.environ.get("CHEMGRAPHX_DIR")
    if not clone or not (Path(clone) / "server" / "column_f.py").exists():
        pytest.skip("set CHEMGRAPHX_DIR to a ChemGraphX clone to re-run it")
    pytest.importorskip("igraph")
    clone = Path(clone)
    if TOOL.clone_commit(clone) != TOOL.CHEMGRAPHX_COMMIT:
        pytest.skip(f"clone is not at the pinned commit {TOOL.CHEMGRAPHX_COMMIT}")
    live = _decode(json.loads(TOOL._dumps(TOOL.build_fixture(clone))))
    assert live["chemgraphx"] == FX["chemgraphx"]
    assert [e["name"] for e in live["panel"]] == list(PANEL)
    for e in live["panel"]:
        rec = PANEL[e["name"]]
        assert e["edges"] == rec["edges"]
        for block, vals in e["chemgraphx"].items():
            for k, x in vals.items():
                assert TOOL.rel_diff(x, rec["chemgraphx"][block][k]) <= 1e-12, (e["name"], k)
    assert live["probes"] == FX["probes"]
