"""Tests for topoqspr.indices.

Layers of evidence, from strongest to weakest:

1. closed forms on paths, cycles, stars and complete graphs;
2. independent oracles: brute-force matching enumeration (Hosoya Z),
   networkx shortest paths / line graphs (distance indices), and Mordred
   (fixture of 24 approved drugs);
3. values printed in the source papers (hydrochlorothiazide, gemcitabine);
4. algebraic identities between indices, invariance under atom renumbering,
   and explicit NaN-with-reason behaviour.
"""

from __future__ import annotations

import itertools
import json
import math
import random
import re
import time
from pathlib import Path

import networkx as nx
import pytest
from rdkit import Chem

from topoqspr.graph import MolGraph, standardize
from topoqspr.indices import (
    FAMILIES,
    INDEX_NAMES,
    PHI,
    REGISTRY,
    RK_VALUES,
    UndefinedIndex,
    compute_indices,
    HOSOYA_MAX_FRONTIER,
    hosoya_z,
    modified_reverse_degree,
)

ROOT = Path(__file__).resolve().parents[1]
FIX = Path(__file__).resolve().parent / "fixtures" / "indices"
DRUGS = json.loads((FIX / "drugs.json").read_text())["drugs"]
MORDRED = json.loads((FIX / "mordred_reference.json").read_text())


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def path(n: int) -> MolGraph:
    return MolGraph.from_edges(n, [(i, i + 1) for i in range(n - 1)])


def cycle(n: int) -> MolGraph:
    return MolGraph.from_edges(n, [(i, (i + 1) % n) for i in range(n)])


def star(n: int) -> MolGraph:
    """K_{1,n-1}: centre 0 and n-1 leaves."""
    return MolGraph.from_edges(n, [(0, i) for i in range(1, n)])


def complete(n: int) -> MolGraph:
    return MolGraph.from_edges(n, itertools.combinations(range(n), 2))


def to_nx(g: MolGraph) -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(g.n))
    G.add_edges_from(g.edges)
    return G


def from_nx(G: nx.Graph) -> MolGraph:
    mapping = {v: i for i, v in enumerate(sorted(G.nodes))}
    return MolGraph.from_edges(len(mapping), [(mapping[u], mapping[v]) for u, v in G.edges])


def fib(k: int) -> int:
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def lucas(k: int) -> int:
    a, b = 2, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def values(g: MolGraph, names=None):
    """compute_indices, asserting the only undefined indices are R_k with k > Delta."""
    v, why = compute_indices(g, names)
    delta = g.max_degree
    expected = {n for n in v for k in RK_VALUES if k > delta and n.startswith(f"R{k}_")}
    assert set(why) == expected
    return v, why


@pytest.fixture(scope="module")
def drug_graphs():
    out = {}
    for d in DRUGS:
        r = standardize(d["smiles"])
        assert r.ok, (d["name"], r.rejection)
        out[d["name"]] = r
    return out


@pytest.fixture(scope="module")
def drug_values(drug_graphs):
    return {k: compute_indices(r.graph) for k, r in drug_graphs.items()}


# ---------------------------------------------------------------------------
# registry integrity
# ---------------------------------------------------------------------------
def test_registry_is_well_formed():
    assert len(INDEX_NAMES) == len(set(INDEX_NAMES)) == len(REGISTRY)
    bib = (ROOT / "docs" / "index_refs.bib").read_text()
    bib_keys = set(re.findall(r"@\w+\{([^,]+),", bib))
    for spec in REGISTRY:
        assert spec.family in FAMILIES, spec.name
        assert spec.latex.strip(), spec.name
        assert spec.ref in bib_keys, f"{spec.name}: bib key {spec.ref} missing"
        assert callable(spec.func)


def test_every_index_is_documented():
    doc = (ROOT / "docs" / "INDEX_VERIFICATION.md").read_text()
    for name in INDEX_NAMES:
        assert f"`{name}`" in doc, f"{name} missing from INDEX_VERIFICATION.md"


def test_every_legacy_name_is_mapped():
    """Each v0.x index name has a row in the legacy mapping table."""
    cfg = ROOT / "config.py"
    if not cfg.exists():
        pytest.skip("legacy config.py not present")
    text = cfg.read_text()
    block = text[text.index("TOPO_INDICES = [") : text.index("]", text.index("TOPO_INDICES = ["))]
    legacy = re.findall(r'"([^"]+)"', block)
    assert len(legacy) == 135
    doc = (ROOT / "docs" / "INDEX_VERIFICATION.md").read_text()
    table = doc[doc.index("## Legacy name mapping") :]
    for name in legacy:
        assert re.search(rf"^\| `{re.escape(name)}` \|", table, re.M), f"legacy {name} unmapped"


def test_values_are_in_registry_order_and_deterministic(drug_graphs):
    g = drug_graphs["imatinib"].graph
    v1, r1 = compute_indices(g)
    v2, r2 = compute_indices(g)
    assert list(v1) == list(INDEX_NAMES)
    assert v1 == v2 and r1 == r2


def test_unknown_index_name_raises():
    with pytest.raises(KeyError):
        compute_indices(path(4), ["M1", "nope"])


def test_subset_selection():
    v, why = compute_indices(path(5), ["W", "M1"])
    assert list(v) == ["M1", "W"]  # registry order
    assert v == {"M1": 14, "W": 20} and why == {}


# ---------------------------------------------------------------------------
# closed forms
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("n", range(3, 31))
def test_path_closed_forms(n):
    v, why = values(path(n))
    assert v["W"] == (n**3 - n) // 6
    assert v["Z"] == fib(n + 1)
    assert v["M1"] == 4 * n - 6
    assert v["M2"] == 4 * n - 8
    assert v["F"] == 8 * n - 14
    assert v["R"] == pytest.approx((n - 3) / 2 + 2 / math.sqrt(2))
    assert v["Sz"] == v["W"]  # Szeged = Wiener for trees
    assert v["PI_v"] == n * (n - 1)  # bipartite: n * m
    assert v["Mo"] == sum(abs(n - 2 * i) for i in range(1, n))
    assert v["DD"] == 4 * v["W"] - n * (n - 1)  # trees (Gutman 1994)
    assert v["Gut"] == 4 * v["W"] - (2 * n - 1) * (n - 1)  # trees (Gutman 1994)
    assert v["W_e"] == ((n - 1) ** 3 - (n - 1)) // 6  # L(P_n) = P_{n-1}
    assert v["J"] == pytest.approx(
        (n - 1) * sum(1 / math.sqrt(s(i, n) * s(i + 1, n)) for i in range(n - 1))
    )


def s(i: int, n: int) -> int:
    """Distance sum of vertex i in P_n."""
    return i * (i + 1) // 2 + (n - 1 - i) * (n - i) // 2


@pytest.mark.parametrize("n", range(3, 31))
def test_cycle_closed_forms(n):
    v, why = values(cycle(n))
    assert v["W"] == (n**3 // 8 if n % 2 == 0 else (n**3 - n) // 8)
    assert v["Z"] == lucas(n)
    assert v["M1"] == v["M2"] == 4 * n
    assert v["Sz"] == (n**3 // 4 if n % 2 == 0 else n * (n - 1) ** 2 // 4)
    assert v["PI"] == (n * (n - 2) if n % 2 == 0 else n * (n - 1))
    assert v["Mo"] == v["Mo_e"] == 0  # distance-balanced
    assert v["J"] == pytest.approx(n**3 / (4 * v["W"]))  # m/(mu+1) = n/2, every D_u = 2W/n
    assert v["I_deg"] == pytest.approx(math.log2(n))  # regular graph
    assert v["sigma"] == 0


@pytest.mark.parametrize("n", range(3, 21))
def test_star_closed_forms(n):
    v, why = values(star(n))
    k = n - 1
    assert v["M1"] == k + k * k
    assert v["M2"] == k * k
    assert v["F"] == k * (1 + k * k)
    assert v["Z"] == n
    assert v["W"] == k * k
    assert v["R"] == pytest.approx(k / math.sqrt(k))
    assert v["sigma"] == k * (k - 1) ** 2
    # Kulli reverse degree c = Delta - d + 1: centre 1, leaves Delta = k
    assert v["R1_M1"] == k * (1 + k)
    assert v["R1_M2"] == k * k
    assert v["HCM2"] == k * k * k


@pytest.mark.parametrize("n", range(3, 16))
def test_cycle_closed_forms_for_every_degree_index(n):
    """C_n is 2-regular: every edge has end degrees (2, 2). Hand-derived values."""
    v, _ = values(cycle(n))
    r2 = math.sqrt(2)
    expected = {
        "M1": 4 * n, "M2": 4 * n, "F": 8 * n, "HM": 16 * n, "R": n / 2, "RR": 2 * n,
        "SC": n / 2, "H": n / 2, "GA": n, "ABC": n / r2, "AZI": 8 * n, "SO": 2 * r2 * n,
        "ISI": n, "SDD": 2 * n, "sigma": 0, "ABS": n / r2, "ReZ3": 16 * n, "AG": n,
        "HG": n / 4, "GBM": n / 4, "BMG": 4 * n, "BMH": 16 * n, "TMH": 24 * n,
    }
    degree_names = {s.name for s in REGISTRY if s.family == "degree"}
    assert set(expected) == degree_names
    for k, x in expected.items():
        assert v[k] == pytest.approx(x, rel=1e-12), k
    # neighbourhood sums are all 4, so S_X(C_n) = X evaluated at (4, 4) per edge
    assert v["S_R"] == pytest.approx(n / 4) and v["S_RR"] == pytest.approx(4 * n)
    assert v["S_SC"] == pytest.approx(n / math.sqrt(8)) and v["S_H"] == pytest.approx(n / 4)
    assert v["S_GA"] == pytest.approx(n) and v["S_SDD"] == pytest.approx(2 * n)
    assert v["S_ABC"] == pytest.approx(n * math.sqrt(6) / 4) and v["S_ReZ3"] == 128 * n
    assert v["S_GBM"] == pytest.approx(n / 6) and v["S_TMH"] == pytest.approx(192 * n)
    # Delta = 2: M_1R(2) = 1 and M_2R(2) = 2 for every vertex
    assert v["R1_M1"] == 2 * n and v["R1_M2"] == n and v["R1_ABC"] == 0.0
    assert v["R2_M1"] == 4 * n and v["R2_M2"] == 4 * n and v["R2_ISI"] == pytest.approx(n)
    assert v["HCM2"] == n


@pytest.mark.parametrize("n", range(3, 16))
def test_star_closed_forms_for_irregular_indices(n):
    k = n - 1  # every edge joins degree 1 to degree k
    v, _ = values(star(n))
    assert v["H"] == pytest.approx(2 * k / (k + 1))
    assert v["GA"] == pytest.approx(2 * k * math.sqrt(k) / (k + 1))
    assert v["SC"] == pytest.approx(k / math.sqrt(k + 1))
    assert v["ABC"] == pytest.approx(k * math.sqrt((k - 1) / k))
    assert v["ABS"] == pytest.approx(k * math.sqrt((k - 1) / (k + 1)))
    assert v["SO"] == pytest.approx(k * math.sqrt(1 + k * k))
    assert v["ISI"] == pytest.approx(k * k / (k + 1))
    assert v["RR"] == pytest.approx(k * math.sqrt(k))
    assert v["I_deg"] == pytest.approx(0.5 + 0.5 * math.log2(2 * k))
    # every S equals k (centre: k leaves of degree 1; leaf: the centre of degree k)
    assert v["S_RR"] == pytest.approx(k * k) and v["S_GA"] == pytest.approx(k)
    if k > 1:
        assert v["AZI"] == pytest.approx(k * (k / (k - 1)) ** 3)


def involutions(n: int) -> int:
    t = [1, 1]
    for i in range(2, n + 1):
        t.append(t[-1] + (i - 1) * t[-2])
    return t[n]


@pytest.mark.parametrize("n", range(3, 12))
def test_complete_graph_closed_forms(n):
    v, why = values(complete(n))
    assert v["W"] == n * (n - 1) // 2
    assert v["M1"] == n * (n - 1) ** 2
    assert v["M2"] == n * (n - 1) ** 3 // 2
    assert v["Z"] == involutions(n)  # matchings of K_n = involutions of n elements
    assert v["Sz"] == n * (n - 1) // 2
    assert v["Mo"] == 0
    assert v["I_deg"] == pytest.approx(math.log2(n))
    # every reverse degree of a regular graph is Delta + k - Delta = k
    for k in RK_VALUES:
        if k <= n - 1:
            assert v[f"R{k}_M2"] == (n * (n - 1) // 2) * k * k


def test_benzene_and_naphthalene_hosoya():
    assert hosoya_z(standardize("c1ccccc1").graph) == 18
    assert hosoya_z(standardize("c1ccc2ccccc2c1").graph) == 148


# ---------------------------------------------------------------------------
# Hosoya Z: brute force and speed
# ---------------------------------------------------------------------------
def brute_force_matchings(g: MolGraph) -> int:
    count = 0
    for r in range(g.m + 1):
        for subset in itertools.combinations(g.edges, r):
            ends = [x for e in subset for x in e]
            if len(ends) == len(set(ends)):
                count += 1
    return count


@pytest.mark.parametrize("seed", range(40))
def test_hosoya_matches_brute_force_on_random_graphs(seed):
    rng = random.Random(seed)
    n = rng.randint(1, 10)
    m = rng.randint(0, min(14, n * (n - 1) // 2))
    g = from_nx(nx.gnm_random_graph(n, m, seed=seed))  # may be disconnected
    assert hosoya_z(g) == brute_force_matchings(g)


def test_hosoya_is_exact_beyond_float_precision():
    z = hosoya_z(path(100))
    assert isinstance(z, int) and z == fib(101) and z > 2**53


def ladder(rungs: int) -> MolGraph:
    L = rungs
    e = [(i, i + 1) for i in range(L - 1)] + [(L + i, L + i + 1) for i in range(L - 1)]
    e += [(i, L + i) for i in range(0, L, 2)]
    return MolGraph.from_edges(2 * L, e)


def test_hosoya_speed_on_drug_size_graphs(drug_graphs):
    big = [
        ladder(75),  # 150 vertices, 37 fused rings
        standardize("N" + "C(Cc1ccccc1)C(=O)N" * 13 + "CC(=O)O").graph,  # 148 heavy atoms
        drug_graphs["vancomycin"].graph,
        drug_graphs["cyclosporine"].graph,
    ]
    for g in big:
        t0 = time.perf_counter()
        hosoya_z(g)
        assert time.perf_counter() - t0 < 0.5, g.n


def test_hosoya_reference_values_for_drugs(drug_graphs):
    # Values recomputed independently by the expert audit (scratch ref.py,
    # memoised edge recursion); the v0.x code returned 452 for imatinib.
    expected = {
        "aspirin": 335,
        "donepezil": 965_380,
        "sitagliptin": 394_051,
        "imatinib": 79_015_504,
        "atorvastatin": 267_895_296,
    }
    for name, z in expected.items():
        assert hosoya_z(drug_graphs[name].graph) == z


# ---------------------------------------------------------------------------
# networkx oracle for distance indices
# ---------------------------------------------------------------------------
def nx_distance_oracle(g: MolGraph) -> dict:
    G = to_nx(g)
    D = dict(nx.all_pairs_shortest_path_length(G))
    V, E = list(G.nodes), [tuple(e) for e in G.edges]
    deg = dict(G.degree)
    pairs = list(itertools.combinations(V, 2))

    def d_edge_vertex(f, x):
        return min(D[f[0]][x], D[f[1]][x])

    sz = sze = pi = piv = mo = moe = 0
    for u, v in E:
        nu = sum(1 for x in V if D[u][x] < D[v][x])
        nv = sum(1 for x in V if D[v][x] < D[u][x])
        mu = sum(1 for f in E if d_edge_vertex(f, u) < d_edge_vertex(f, v))
        mv = sum(1 for f in E if d_edge_vertex(f, v) < d_edge_vertex(f, u))
        sz += nu * nv
        sze += mu * mv
        pi += mu + mv
        piv += nu + nv
        mo += abs(nu - nv)
        moe += abs(mu - mv)
    tr = {x: sum(D[x].values()) for x in V}
    m, n = G.number_of_edges(), G.number_of_nodes()
    J = m / (m - n + 2) * sum((tr[u] * tr[v]) ** -0.5 for u, v in E)
    return {
        "W": nx.wiener_index(G),
        "W_e": nx.wiener_index(nx.line_graph(G)),
        "J": J,
        "Sz": sz,
        "Sz_e": sze,
        "PI": pi,
        "PI_v": piv,
        "Mo": mo,
        "Mo_e": moe,
        "DD": sum((deg[a] + deg[b]) * D[a][b] for a, b in pairs),
        "Gut": sum(deg[a] * deg[b] * D[a][b] for a, b in pairs),
    }


def _check_against_oracle(g: MolGraph):
    v, why = compute_indices(g)
    ref = nx_distance_oracle(g)
    for k, x in ref.items():
        assert v[k] == pytest.approx(x, rel=1e-12), k
    assert not (set(ref) & set(why))


@pytest.mark.parametrize("seed", range(25))
def test_distance_indices_match_networkx_on_random_connected_graphs(seed):
    rng = random.Random(seed)
    n = rng.randint(3, 16)
    while True:
        G = nx.gnm_random_graph(n, rng.randint(n - 1, min(n * (n - 1) // 2, 2 * n)), seed=rng.randrange(10**9))
        if nx.is_connected(G):
            break
    _check_against_oracle(from_nx(G))


def test_distance_indices_match_networkx_on_drugs(drug_graphs):
    for name in ("aspirin", "caffeine", "morphine", "imatinib", "memantine", "sildenafil"):
        _check_against_oracle(drug_graphs[name].graph)


# ---------------------------------------------------------------------------
# Mordred oracle (fixture recorded by tests/fixtures/indices/make_fixtures.py)
# ---------------------------------------------------------------------------
def test_mordred_fixture_covers_at_least_15_drugs():
    assert len(MORDRED["values"]) >= 15
    assert set(MORDRED["values"]) == {d["name"] for d in DRUGS}


@pytest.mark.parametrize("drug", [d["name"] for d in DRUGS])
def test_overlapping_indices_match_mordred(drug, drug_graphs):
    ref = MORDRED["values"][drug]
    r = drug_graphs[drug]
    assert r.inchikey == ref["inchikey"]
    v, _ = compute_indices(r.graph, MORDRED["mapping"].values())
    for mordred_name, ours in MORDRED["mapping"].items():
        assert v[ours] == pytest.approx(ref["mordred"][mordred_name], rel=1e-10, abs=1e-12), (
            drug,
            ours,
        )


def test_mordred_fixture_is_reproducible(drug_graphs):
    """Recompute the fixture live, so a stale fixture cannot hide a bug."""
    mordred = pytest.importorskip("mordred")
    from mordred import ABCIndex, BalabanJ, Calculator, Chi, WienerIndex, ZagrebIndex

    del mordred
    calc = Calculator(
        [
            ABCIndex.ABCIndex(),
            BalabanJ.BalabanJ(),
            WienerIndex.WienerIndex(False),
            ZagrebIndex.ZagrebIndex(1, 1),
            ZagrebIndex.ZagrebIndex(2, 1),
            Chi.Chi("path", 1, "d", False),
        ]
    )
    names = [str(d) for d in calc.descriptors]
    for drug in ("aspirin", "imatinib", "vancomycin", "morphine", "memantine"):
        live = dict(zip(names, (float(x) for x in calc(drug_graphs[drug].mol))))
        for k, x in MORDRED["values"][drug]["mordred"].items():
            assert live[k] == pytest.approx(x, rel=1e-12), (drug, k)


# ---------------------------------------------------------------------------
# values printed in the source papers
# ---------------------------------------------------------------------------
HCTZ = "C1NC2=CC(=C(C=C2S(=O)(=O)N1)S(=O)(=O)N)Cl"
GEMCITABINE = "C1=CN(C(=O)N=C1N)[C@H]2C([C@@H]([C@H](O2)CO)O)(F)F"

# Godlin & Radha, Sci. Rep. 16 (2026) 8889: worked example for
# hydrochlorothiazide (G17), values rounded to 2 decimals; k = 0 is the plain
# degree version. ReZ1 = (a+b)/(ab), ReZ2 = ab/(a+b) (= ISI), SSD = SDD.
_PAPER_PHI = {
    "ABC": PHI["ABC"][0],
    "R": PHI["R"][0],
    "ReZ1": lambda a, b: (a + b) / (a * b),
    "ReZ2": PHI["ISI"][0],
    "ReZ3": PHI["ReZ3"][0],
    "SSD": PHI["SDD"][0],
    "AG": PHI["AG"][0],
    "HG": PHI["HG"][0],
    "BMG": PHI["BMG"][0],
    "GBM": PHI["GBM"][0],
}
GODLIN_HCTZ = {
    0: dict(ABC=13.43, R=7.72, ReZ1=17, ReZ2=20.51, ReZ3=618, SSD=50.08, AG=19.59, HG=3.11, BMG=82.91, GBM=3.94),
    1: dict(ABC=13.55, R=8.55, ReZ1=18.83, ReZ2=18.42, ReZ3=462, SSD=50.92, AG=19.69, HG=3.89, BMG=79.02, GBM=4.14),
    2: dict(ABC=12.26, R=7.89, ReZ1=16.5, ReZ2=23.39, ReZ3=934, SSD=41.08, AG=18.61, HG=3.85, BMG=85.07, GBM=3.95),
    3: dict(ABC=11.91, R=8.55, ReZ1=18.67, ReZ2=20.51, ReZ3=738, SSD=50.08, AG=19.59, HG=4.51, BMG=82.91, GBM=4.04),
}


def test_hydrochlorothiazide_edge_partition_matches_godlin_2026():
    g = standardize(HCTZ).graph
    part = {}
    for u, v in g.edges:
        key = tuple(sorted((g.degrees[u], g.degrees[v])))
        part[key] = part.get(key, 0) + 1
    assert part == {(1, 3): 1, (1, 4): 5, (2, 2): 2, (2, 3): 5, (2, 4): 1, (3, 3): 2, (3, 4): 2}


@pytest.mark.parametrize("k", sorted(GODLIN_HCTZ))
def test_hydrochlorothiazide_matches_godlin_2026(k):
    g = standardize(HCTZ).graph
    w = g.degrees if k == 0 else [modified_reverse_degree(d, k, g.max_degree) for d in g.degrees]
    for name, expected in GODLIN_HCTZ[k].items():
        phi = _PAPER_PHI[name]
        got = sum(phi(w[u], w[v]) for u, v in g.edges)
        assert round(got, 2) == pytest.approx(expected), (k, name, got)


def test_registry_reproduces_godlin_2026():
    v, _ = compute_indices(standardize(HCTZ).graph)
    assert round(v["ABC"], 2) == 13.43 and v["ReZ3"] == 618 and round(v["ISI"], 2) == 20.51
    assert round(v["AG"], 2) == 19.59 and round(v["HG"], 2) == 3.11
    assert round(v["BMG"], 2) == 82.91 and round(v["GBM"], 2) == 3.94 and round(v["SDD"], 2) == 50.08
    for k, abc, r, rez3 in [(1, 13.55, 8.55, 462), (2, 12.26, 7.89, 934), (3, 11.91, 8.55, 738)]:
        assert round(v[f"R{k}_ABC"], 2) == abc
        assert round(v[f"R{k}_R"], 2) == r
        assert v[f"R{k}_ReZ3"] == rez3


def test_gemcitabine_matches_arockiaraj_2025():
    """Arockiaraj et al., Sci. Rep. 15 (2025) 3639: ABC^d, ABC^nd, ABC^{m_k d}."""
    v, _ = compute_indices(standardize(GEMCITABINE).graph)
    assert round(v["ABC"], 2) == 13.84
    assert round(v["S_ABC"], 2) == 10.36
    assert round(v["R1_ABC"], 2) == 13.65
    assert round(v["R2_ABC"], 2) == 13.27
    assert round(v["R3_ABC"], 2) == 13.43


# ---------------------------------------------------------------------------
# modified reverse degree
# ---------------------------------------------------------------------------
def test_modified_reverse_degree_tables_from_the_sources():
    # Delta = 4 tables printed by Kalaam et al. 2024 (Sec. 3) and Godlin & Radha 2026
    table = {1: (4, 3, 2, 1), 2: (1, 4, 3, 2), 3: (2, 1, 4, 3)}
    for k, row in table.items():
        assert tuple(modified_reverse_degree(d, k, 4) for d in (1, 2, 3, 4)) == row


# Edge functions written out again, independently of indices.PHI, in the
# notation of Kalaam et al. 2024 Eqs. (1)-(9) and Godlin & Radha 2026.
_RK_FORMULAS = {
    "M1": lambda a, b: a + b,
    "M2": lambda a, b: a * b,
    "F": lambda a, b: a**2 + b**2,
    "HM": lambda a, b: (a + b) ** 2,  # "HZ"
    "ReZ3": lambda a, b: (a + b) * a * b,
    "BM": lambda a, b: a + b + a * b,
    "TM": lambda a, b: a**2 + b**2 + a * b,
    "BMH": lambda a, b: (a + b + a * b) * (a + b) / 2,
    "TMH": lambda a, b: (a**2 + b**2 + a * b) * (a + b) / 2,
    "R": lambda a, b: (a * b) ** -0.5,
    "ABC": lambda a, b: ((a + b - 2) / (a * b)) ** 0.5,
    "ISI": lambda a, b: a * b / (a + b),  # "ReZ_2"
    "RR": lambda a, b: (a * b) ** 0.5,  # "G"
}


@pytest.mark.parametrize("smiles", [HCTZ, GEMCITABINE, "CC(C)(C)c1ccc(O)cc1"])
def test_rk_indices_use_the_published_weight_tables(smiles):
    """All R_k indices on Delta = 4 molecules, with weights looked up in the
    printed tables (not computed by modified_reverse_degree)."""
    table = {1: (4, 3, 2, 1), 2: (1, 4, 3, 2), 3: (2, 1, 4, 3)}
    g = standardize(smiles).graph
    assert g.max_degree == 4
    v, why = compute_indices(g)
    assert not any(n.startswith("R") and n[1].isdigit() for n in why)
    rk_names = {s.name for s in REGISTRY if s.family == "reverse_degree" and s.name != "HCM2"}
    seen = set()
    for k, row in table.items():
        w = [row[d - 1] for d in g.degrees]
        for key, f in _RK_FORMULAS.items():
            name = f"R{k}_{key}"
            seen.add(name)
            assert v[name] == pytest.approx(sum(f(w[a], w[b]) for a, b in g.edges), rel=1e-12), name
        if k == 1:
            assert v["HCM2"] == sum((w[a] * w[b]) ** 2 for a, b in g.edges)
    assert seen == rk_names


@pytest.mark.parametrize("delta", range(1, 9))
def test_modified_reverse_degree_is_a_permutation_for_k_up_to_delta(delta):
    for k in range(1, delta + 1):
        img = sorted(modified_reverse_degree(d, k, delta) for d in range(1, delta + 1))
        assert img == list(range(1, delta + 1))
    # k = 1 is Kulli's reverse degree c = Delta - d + 1
    assert [modified_reverse_degree(d, 1, delta) for d in range(1, delta + 1)] == [
        delta - d + 1 for d in range(1, delta + 1)
    ]


@pytest.mark.parametrize("k", [0, 4, 5])
def test_modified_reverse_degree_outside_range_is_undefined(k):
    with pytest.raises(UndefinedIndex):
        modified_reverse_degree(1, k, 3)


# ---------------------------------------------------------------------------
# NaN, never zero
# ---------------------------------------------------------------------------
def _nan_keys(v):
    return {k for k, x in v.items() if isinstance(x, float) and math.isnan(x)}


def test_reasons_are_given_exactly_for_nan_values(drug_values):
    for name, (v, why) in drug_values.items():
        assert _nan_keys(v) == set(why), name


def test_drugs_have_no_undefined_indices_except_rk_beyond_delta(drug_graphs, drug_values):
    for name, (v, why) in drug_values.items():
        delta = drug_graphs[name].graph.max_degree
        expected = {n for n in INDEX_NAMES for k in RK_VALUES if k > delta and n.startswith(f"R{k}_")}
        assert set(why) == expected, name


def test_rk_beyond_delta_is_nan_with_reason():
    g = path(5)  # Delta = 2
    v, why = compute_indices(g)
    for k in RK_VALUES:
        for spec in REGISTRY:
            if spec.name.startswith(f"R{k}_"):
                if k > 2:
                    assert math.isnan(v[spec.name]) and "outside 1..Delta" in why[spec.name]
                else:
                    assert spec.name not in why
    assert v["R3_M1"] != 0


def test_division_by_zero_is_nan_with_reason():
    g = path(2)  # single edge (1,1): AZI has d_u + d_v - 2 = 0
    v, why = compute_indices(g)
    assert math.isnan(v["AZI"]) and "division by zero" in why["AZI"]
    assert v["ABC"] == 0.0 and "ABC" not in why  # sqrt(0/1) = 0 is a true zero
    assert v["M1"] == 2


def test_disconnected_graph_distance_indices_are_nan():
    g = MolGraph.from_edges(5, [(0, 1), (1, 2), (3, 4)])
    v, why = compute_indices(g)
    for name in ("W", "J", "Sz", "Sz_e", "PI", "PI_v", "Mo", "Mo_e", "W_e", "DD", "Gut"):
        assert math.isnan(v[name]) and "components" in why[name], name
    assert v["M1"] == 8 and v["Z"] == 3 * 2  # Z(P3) * Z(P2)
    assert "M1" not in why


def test_edgeless_graph():
    g = MolGraph.from_edges(3, [])
    v, why = compute_indices(g)
    assert math.isnan(v["I_deg"]) and "I_deg" in why
    assert v["Z"] == 1 and v["M1"] == 0  # true zero: no edges


def test_neighbourhood_sum_zero_is_nan():
    # isolated vertex has S = 0 but no incident edge, so no undefined term;
    # an edge between two leaves of K2 has S = 1 each.
    g = MolGraph.from_edges(3, [(0, 1)])
    v, why = compute_indices(g)
    assert v["S_R"] == 1.0 and "S_R" not in why


# ---------------------------------------------------------------------------
# identities and invariance
# ---------------------------------------------------------------------------
def test_algebraic_identities(drug_values, drug_graphs):
    for name, (v, _) in drug_values.items():
        g = drug_graphs[name].graph
        assert v["sigma"] == v["F"] - 2 * v["M2"]
        assert v["HM"] == v["F"] + 2 * v["M2"]
        # sum over vertices of d^2 equals the edge form of M1
        assert v["M1"] == sum(d * d for d in g.degrees)
        # R_1 uses Kulli's c = Delta - d + 1: CM1 = sum_v d_v c_v
        delta = g.max_degree
        assert v["R1_M1"] == sum(d * (delta - d + 1) for d in g.degrees)
        for k in RK_VALUES:
            if k > delta:
                continue
            assert v[f"R{k}_HM"] == v[f"R{k}_F"] + 2 * v[f"R{k}_M2"]
            assert v[f"R{k}_BM"] == v[f"R{k}_M1"] + v[f"R{k}_M2"]
            assert v[f"R{k}_TM"] == v[f"R{k}_F"] + v[f"R{k}_M2"]
        # redefined first Zagreb sum (1/d_u + 1/d_v) = n, which is why it is not registered
        assert sum(1 / g.degrees[a] + 1 / g.degrees[b] for a, b in g.edges) == pytest.approx(g.n)


@pytest.mark.parametrize("drug", ["aspirin", "morphine", "imatinib", "sitagliptin", "paclitaxel"])
def test_invariance_under_atom_renumbering(drug, drug_graphs):
    smi = next(d["smiles"] for d in DRUGS if d["name"] == drug)
    base, _ = compute_indices(drug_graphs[drug].graph)
    mol = Chem.MolFromSmiles(smi)
    rng = random.Random(drug)
    for _ in range(3):
        order = list(range(mol.GetNumAtoms()))
        rng.shuffle(order)
        shuffled = Chem.MolToSmiles(Chem.RenumberAtoms(mol, order), canonical=False)
        r = standardize(shuffled)
        assert r.inchikey == drug_graphs[drug].inchikey
        v, _ = compute_indices(r.graph)
        for k in INDEX_NAMES:
            if isinstance(base[k], int):
                assert v[k] == base[k], k
            else:
                assert v[k] == pytest.approx(base[k], rel=1e-12, nan_ok=True), k


def test_invariance_under_graph_relabelling(drug_graphs):
    g = drug_graphs["donepezil"].graph
    base, _ = compute_indices(g)
    perm = list(range(g.n))
    random.Random(7).shuffle(perm)
    v, _ = compute_indices(g.relabel(perm))
    for k in INDEX_NAMES:
        assert v[k] == pytest.approx(base[k], rel=1e-12, nan_ok=True), k


def test_no_two_indices_are_identical_on_the_drug_panel(drug_values):
    """Guards against re-introducing exact duplicates such as v0.x DS1 == DS2."""
    names = [n for n in INDEX_NAMES]
    cols = {n: [drug_values[d][0][n] for d in sorted(drug_values)] for n in names}
    for a, b in itertools.combinations(names, 2):
        xa, xb = cols[a], cols[b]
        same = all(
            (math.isnan(p) and math.isnan(q)) if (isinstance(p, float) and math.isnan(p)) else p == pytest.approx(q, rel=1e-12)
            for p, q in zip(xa, xb)
        )
        assert not same, f"{a} and {b} are identical on every drug"


def test_legacy_defects_are_fixed(drug_graphs):
    """Concrete v0.x bugs confirmed by running the legacy ml.py."""
    v, _ = compute_indices(drug_graphs["aspirin"].graph)
    # v0.x SS_ABC used (S_u + S_v + S_u S_v - 2) and gave 15.02 for aspirin
    assert round(v["S_ABC"], 2) == 7.51
    # v0.x Hosoya returned 1 for 19..25 edges and 1+m+m(m-1)/4 above 25
    assert hosoya_z(path(21)) == fib(22) == 17711
    v, _ = compute_indices(drug_graphs["imatinib"].graph)
    assert v["Z"] == 79_015_504


# ---------------------------------------------------------------------------
# regressions from the v1 verification round
# ---------------------------------------------------------------------------
def test_godlin_hctz_table_is_complete():
    # Godlin & Radha 2026, Table 5: 10 indices for degree and l = 1, 2, 3
    assert sum(len(row) for row in GODLIN_HCTZ.values()) == 40
    assert all(len(row) == 10 for row in GODLIN_HCTZ.values())


@pytest.mark.parametrize("n", [0, 1])
def test_rk_is_nan_on_graphs_without_edges_including_empty(n):
    v, why = compute_indices(MolGraph.from_edges(n, []))
    for name in ("R1_M1", "R3_ABC", "HCM2"):
        assert math.isnan(v[name]), (n, name, v[name])
        assert "Delta(G) >= 1" in why[name]


def _mg(G) -> MolGraph:
    G = nx.convert_node_labels_to_integers(G)
    return MolGraph.from_edges(G.number_of_nodes(), G.edges())


def test_hosoya_frontier_cap_gives_nan_with_reason_quickly():
    g = _mg(nx.random_regular_graph(3, 100, seed=0))  # width 18, ~17 s uncapped
    t0 = time.perf_counter()
    with pytest.raises(UndefinedIndex, match=r"frontier width 18 exceeds the cap 12"):
        hosoya_z(g)
    v, why = compute_indices(g, ["Z"])
    assert time.perf_counter() - t0 < 1.0
    assert math.isnan(v["Z"]) and "frontier width" in why["Z"]


def test_hosoya_cap_is_configurable_and_exact_when_lifted():
    g = _mg(nx.grid_2d_graph(3, 4))  # 12 vertices, 17 edges
    with pytest.raises(UndefinedIndex):
        hosoya_z(g, max_frontier=1)
    assert hosoya_z(g, max_frontier=None) == hosoya_z(g) == brute_force_matchings(g)


def test_hosoya_default_cap_admits_c60():
    c60 = standardize(
        "c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4"
        "c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41"
    ).graph
    assert HOSOYA_MAX_FRONTIER >= 10
    assert hosoya_z(c60) == 1417036634543488
