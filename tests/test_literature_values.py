"""Independent literature check of the reverse-degree and neighbourhood-degree-sum indices.

Every number asserted here is printed in a published paper that is *not* one of the
sources the indices were implemented from (kulli2018, arockiaraj2023ijqc, kalaam2024,
godlin2026, arockiaraj2025lung, mondal2021, arockiaraj2023sar), and most come from
groups other than Arockiaraj's. The provenance, the definitions as printed and the
classification of every comparison (agrees / different convention / error in the paper)
are in ``docs/LITERATURE_VALUE_CHECK.md``. Transcribed tables are in
``tests/fixtures/literature/*.json``.

Three levels of evidence are used:

* full pipeline: the molecule's SMILES goes through ``standardize`` and
  ``compute_indices`` (or ``PHI`` on ``GraphContext.nbr_sum`` for indices that are
  not registered), and the printed value is reproduced;
* explicit graph: a graph family the paper defines exactly (cacti, Cayley trees) is
  built and the paper's closed forms are reproduced for several parameter values;
* edge partition: where only the paper's edge partition is available, ``PHI`` and
  ``modified_reverse_degree`` reproduce the printed numbers from it.

Values the paper gets wrong are asserted as *discrepancies* (our value, the paper's
value and the reason), so the numbers quoted in the document are regression-tested too.
All tests are offline.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

import pytest

from topoqspr.graph import MolGraph, standardize
from topoqspr.indices import PHI, GraphContext, compute_indices, modified_reverse_degree

FIX = Path(__file__).resolve().parent / "fixtures" / "literature"


def _load(name: str) -> dict:
    return json.loads((FIX / name).read_text())


YOUSSEF = _load("youssef2024_fete2.json")
ALBADRANI = _load("albadrani2022_graphene.json")
RAO = _load("rao2024_reverse_zagreb.json")
ABUBAKAR = _load("abubakar2024_antitb.json")
ALMULLA = _load("almulla2026_antifungal.json")
RAVI2023 = _load("ravi2023_reduced_reverse.json")
RAVI2025 = _load("ravi2025_parkinson.json")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _decimals(x) -> int:
    if isinstance(x, int):
        return -1  # exact
    s = repr(float(x))
    return len(s.split(".")[1]) if "." in s else 0


def _matches_printed(ours: float, printed) -> bool:
    """True if ``ours`` rounds to ``printed`` at the number of decimals printed.

    Integers are compared exactly. A half-unit tolerance (plus 1e-9) is used for
    decimals, so both rounding and truncation by the authors are accepted only when
    they coincide with the printed digits.
    """
    d = _decimals(printed)
    if d < 0:
        return abs(ours - printed) < 1e-9
    return abs(ours - printed) <= 0.5 * 10 ** (-d) + 1e-9


def _within_last_digit(ours: float, printed) -> bool:
    """True if ``ours`` is within one unit of the last printed digit: the authors truncated instead
    of rounding, or rounded intermediate values; listed row by row below, never applied silently."""
    d = _decimals(printed)
    return d >= 0 and abs(ours - printed) <= 10 ** (-d) + 1e-9


# Printed values that differ from ours by at most one unit of the last printed digit (truncation
# or rounding of intermediates by the authors); every other printed value matches at its printed
# precision. Each entry is (table id, column).
_LAST_DIGIT = {("C13", "SC"), ("C24", "SC"), ("C19", "AG"), ("C28", "SC"),      # Albadrani 2022 Table 3
               ("Ethambutol", "NGA"), ("Linezolid", "NGA"), ("Linezolid", "NRR"), ("Pyrazinamide", "NGA")}  # Abubakar 2024 Table 3


def _check_printed(row_id: str, col: str, ours: float, printed, paper_errors: dict) -> None:
    if (row_id, col) in paper_errors:
        assert not _matches_printed(ours, printed)
        assert ours == pytest.approx(paper_errors[(row_id, col)], abs=5e-5)
    elif (row_id, col) in _LAST_DIGIT:
        assert not _matches_printed(ours, printed) and _within_last_digit(ours, printed), (row_id, col, ours, printed)
    else:
        assert _matches_printed(ours, printed), (row_id, col, ours, printed)


def _graph(smiles: str) -> MolGraph:
    res = standardize(smiles)
    assert res.graph is not None, (smiles, res)
    return res.graph


def _parse_pair(key: str) -> tuple[int, int]:
    a, b = key.strip("()").split(",")
    return tuple(sorted((int(a), int(b))))  # type: ignore[return-value]


def _pairs(g: MolGraph, w) -> Counter:
    return Counter(tuple(sorted((w[u], w[v]))) for u, v in g.edges)


def _nbr(g: MolGraph, key: str) -> float:
    """``sum_{uv} phi(S_u, S_v)`` with TopoQSPR's edge function ``PHI[key]`` and ``S`` from ``GraphContext``."""
    f, S = PHI[key][0], GraphContext(g).nbr_sum
    return sum(f(S[u], S[v]) for u, v in g.edges)


def _deg(g: MolGraph, key: str) -> float:
    f, d = PHI[key][0], g.degrees
    return sum(f(d[u], d[v]) for u, v in g.edges)


def _kulli(g: MolGraph, key: str) -> float:
    """``sum_{uv} phi(c_u, c_v)`` with ``c = M_1R(d) = Delta - d + 1`` from ``modified_reverse_degree``."""
    f, D = PHI[key][0], g.max_degree
    c = [modified_reverse_degree(x, 1, D) for x in g.degrees]
    return sum(f(c[u], c[v]) for u, v in g.edges)


def _from_partition(part: dict, key_or_func) -> float:
    f = PHI[key_or_func][0] if isinstance(key_or_func, str) else key_or_func
    return sum(n * f(a, b) for (a, b), n in part.items())


# ---------------------------------------------------------------------------
# Youssef et al., Sci. Rep. 14 (2024) 20533 -- FeTe2 network, Kulli reverse degree
# ---------------------------------------------------------------------------
def _fete2_degree_partition(m: int, n: int) -> dict:
    # Table 1 (degree pairs); quoted in the fixture
    return {(1, 2): 2, (1, 3): n + 1, (2, 2): 2 * (m - 1), (2, 3): 2 * (-2 + 2 * m + n),
            (3, 3): 3 * (2 * m * n - n + 1) - 5 * m}


@pytest.mark.parametrize("m", range(1, 9))
@pytest.mark.parametrize("n", range(1, 9))
def test_youssef_2024_partition_is_consistent_with_its_vertex_counts(m, n):
    part = _fete2_degree_partition(m, n)
    ends = Counter()
    for (a, b), c in part.items():
        ends[a] += c
        ends[b] += c
    # "the number of vertices of degree 1 are n + 3, degree 2 are 4m + n - 3, and degree 3 are 4mn - 2m - n + 1"
    assert ends[1] == n + 3
    assert ends[2] == 2 * (4 * m + n - 3)
    assert ends[3] == 3 * (4 * m * n - 2 * m - n + 1)
    assert sum(part.values()) == 6 * m * n + m  # "size ... 6mn + m"


_YOUSSEF_KEYS = {"RR1": "M2", "RRhalf": "RR", "RRmhalf": "R", "ABC": "ABC", "GA": "GA",
                 "M1": "M1", "M2": "M2", "HM": "HM", "F": "F"}


@pytest.mark.parametrize("row", YOUSSEF["rows_tables_3_4_5"], ids=lambda r: f"m=n={r['m']}")
def test_youssef_2024_fete2_tables_3_to_5(row):
    """Tables 3-5: reverse-degree indices of FeTe2 from the Table 1 partition and M_1R = Kulli c."""
    m, n = row["m"], row["n"]
    rev = {}
    for (a, b), c in _fete2_degree_partition(m, n).items():
        key = tuple(sorted((modified_reverse_degree(a, 1, 3), modified_reverse_degree(b, 1, 3))))
        rev[key] = rev.get(key, 0) + c
    for col, key in _YOUSSEF_KEYS.items():
        assert _from_partition(rev, key) == pytest.approx(row["values"][col], rel=2e-6, abs=2e-6), col
    # rho = -1 (not registered) as an arithmetic check of the same partition
    assert _from_partition(rev, lambda a, b: 1 / (a * b)) == pytest.approx(row["values"]["RRm1"], rel=2e-6)


def test_youssef_2024_reverse_partition_is_table_2():
    rev = Counter()
    for (a, b), c in _fete2_degree_partition(3, 2).items():
        rev[tuple(sorted((modified_reverse_degree(a, 1, 3), modified_reverse_degree(b, 1, 3))))] += c
    m, n = 3, 2
    assert dict(rev) == {(2, 3): 2, (1, 3): n + 1, (2, 2): 2 * (m - 1), (1, 2): 2 * (-2 + 2 * m + n),
                         (1, 1): 3 * (2 * m * n - n + 1) - 5 * m}


# ---------------------------------------------------------------------------
# Albadrani et al., Materials 15 (2022) 2889 -- graphene sheets, Kulli reverse degree
# ---------------------------------------------------------------------------
_ALB_REGISTRY = {"M1": "R1_M1", "M2": "R1_M2", "R": "R1_R", "HM": "R1_HM", "GA": "R1_GA", "AG": "R1_AG",
                 "ISI": "R1_ISI", "ReZ1": "R1_ReZ1", "ReZ3": "R1_ReZ3"}
_ALB_PHI_ONLY = {"SC": "SC", "H": "H"}  # not registered for R_k; checked with PHI on the same reverse degrees
# printed values that disagree with the paper's own partition (see the document)
_ALB_PAPER_ERRORS = {("C32", "GA"): 40.0849, ("C30", "SC"): 23.6037, ("C30", "H"): 29.8333}


@pytest.mark.parametrize("row", [r for r in ALBADRANI["rows_table_3"] if "smiles" in r], ids=lambda r: r["sheet"])
def test_albadrani_2022_table_3_full_pipeline(row):
    g = _graph(row["smiles"])
    assert g.n == int(row["sheet"][1:]) and len(g.edges) == g.n + row["rings"] - 1
    vals, reasons = compute_indices(g, list(_ALB_REGISTRY.values()))
    assert not reasons
    for col, printed in row["values"].items():
        ours = vals[_ALB_REGISTRY[col]] if col in _ALB_REGISTRY else _kulli(g, _ALB_PHI_ONLY[col])
        _check_printed(row["sheet"], col, ours, printed, _ALB_PAPER_ERRORS)


def _solve_reverse_partition(row) -> list[tuple[int, int, int]]:
    """Kulli reverse degrees for Delta = 3: deg 3 -> 1, deg 2 -> 2. Solve RM1, RM2 and m for (x11, x12, x22)."""
    m = int(row["sheet"][1:]) + row["rings"] - 1
    v = row["values"]
    return [(x11, x12, m - x11 - x12) for x11 in range(m + 1) for x12 in range(m + 1 - x11)
            if 2 * x11 + 3 * x12 + 4 * (m - x11 - x12) == v["M1"] and x11 + 2 * x12 + 4 * (m - x11 - x12) == v["M2"]]


@pytest.mark.parametrize("sheet", ["C19", "C22", "C28", "C30"])
def test_albadrani_2022_table_3_partition_level(sheet):
    """Sheets without an unambiguous name: the partition implied by RM1, RM2, m reproduces the other columns."""
    row = next(r for r in ALBADRANI["rows_table_3"] if r["sheet"] == sheet)
    sols = _solve_reverse_partition(row)
    assert len(sols) == 1
    x11, x12, x22 = sols[0]
    part = {(1, 1): x11, (1, 2): x12, (2, 2): x22}
    for col, printed in row["values"].items():
        ours = _from_partition(part, col if col in PHI else {"ReZ1": "ReZ1"}[col])
        _check_printed(sheet, col, ours, printed, _ALB_PAPER_ERRORS)
    # handshake: deg-2 vertices = (2 x22 + x12)/2 (reverse 2 = degree 2) must be an integer
    n2_twice = 2 * x22 + x12
    if sheet == "C22":
        # (x11, x12, x22) = (10, 9, 8) is the only partition with the printed RM1 = 79, RM2 = 60 and m = 27,
        # but it needs 12.5 vertices of degree 2: no graph has it (paper error in the C22 row).
        assert sols[0] == (10, 9, 8) and n2_twice % 2 == 1
    else:
        assert n2_twice % 2 == 0


def _graphene_case1(l: int, k: int) -> dict:
    # Theorem 1, case l > 1, k > 1 (also used by the paper for k = 1): reverse-degree partition
    return {(1, 1): 3 * l * k - 2 * k - l - 1, (1, 2): 4 * k + 2 * l - 4, (2, 2): l + 4}


@pytest.mark.parametrize("l", [2, 3, 5])
@pytest.mark.parametrize("k", [2, 3, 4])
def test_albadrani_2022_remark_1_closed_forms(l, k):
    p = _graphene_case1(l, k)
    assert _from_partition(p, "M1") == 6 * l * k + 8 * k + 8 * l + 2
    assert _from_partition(p, "M2") == 3 * l * k + 6 * k + 7 * l + 7
    assert _from_partition(p, "HM") == 12 * l * k + 28 * k + 30 * l + 24
    assert _from_partition(p, "ReZ3") == 6 * l * k + 20 * k + 26 * l + 38


@pytest.mark.parametrize("k", [2, 3, 7])
def test_albadrani_2022_remark_2_acenes(k):
    """l = 1, k > 1 (linear acenes); partition (k-1, 4k-4, 6). Two printed constants are wrong."""
    p = {(1, 1): k - 1, (1, 2): 4 * k - 4, (2, 2): 6}
    s2 = math.sqrt(2)
    # also the full pipeline on the acene itself: naphthalene (k=2), anthracene (k=3), ...
    smiles = "c1ccc2cc3ccccc3cc2c1" if k == 3 else ("c1ccc2ccccc2c1" if k == 2 else None)
    if smiles:
        g = _graph(smiles)
        assert _pairs(g, [modified_reverse_degree(d, 1, 3) for d in g.degrees]) == Counter(
            {key: v for key, v in p.items() if v})
    assert _from_partition(p, "M1") == 14 * k + 10
    assert _from_partition(p, "M2") == 9 * k + 15
    assert _from_partition(p, "HM") == 40 * k + 56
    assert _from_partition(p, "ISI") == pytest.approx(19 * k / 6 + 17 / 6)
    assert _from_partition(p, "AG") == pytest.approx((1 + 3 * s2) * k + 5 - 3 * s2)
    assert _from_partition(p, "GA") == pytest.approx((3 + 8 * s2) / 3 * k + 5 - 8 * s2 / 3)
    assert _from_partition(p, "ReZ3") == 26 * k + 70
    # Remark 2 (iii) prints RR = (1 + 2 sqrt2) k + (3 - 4/sqrt2); the partition gives constant 2 - 2 sqrt2
    rr = _from_partition(p, "R")
    assert rr == pytest.approx((1 + 2 * s2) * k + 2 - 2 * s2)
    assert rr - ((1 + 2 * s2) * k + 3 - 4 / s2) == pytest.approx(-1)
    # Remark 2 (x) prints RReZG1 = 8k - 1/2; the partition gives 8k - 2
    assert _from_partition(p, "ReZ1") == pytest.approx(8 * k - 2)


@pytest.mark.parametrize("l", [2, 3, 6])
def test_albadrani_2022_remark_3_first_reverse_zagreb_is_misprinted(l):
    p = _graphene_case1(l, 1)
    assert p == {(1, 1): 2 * l - 3, (1, 2): 2 * l, (2, 2): l + 4}  # the paper's Case 3 partition
    assert _from_partition(p, "M1") == 14 * l + 10  # printed: 14 l - 2
    assert _from_partition(p, "M2") == 10 * l + 13  # printed, agrees
    assert _from_partition(p, "HM") == 42 * l + 52  # printed, agrees
    assert _from_partition(p, "ReZ3") == 32 * l + 58  # printed, agrees


# ---------------------------------------------------------------------------
# Rao, Chen, Ahmad, Ahmad, ACS Omega 9 (2024) 31056 -- Kulli reverse Zagreb indices
# ---------------------------------------------------------------------------
_RAO_REG = {"M1": "M1", "M2": "M2", "HM1": "HM", "CM1": "R1_M1", "CM2": "R1_M2", "HCM1": "R1_HM", "HCM2": "HCM2"}


@pytest.mark.parametrize("comp", [c for c in RAO["compounds_table_3"] if "smiles" in c], ids=lambda c: c["name"])
def test_rao_2024_table_3_full_pipeline(comp):
    g = _graph(comp["smiles"])
    vals, reasons = compute_indices(g, list(_RAO_REG.values()))
    assert not reasons
    for col, printed in comp["values"].items():
        ours = vals[_RAO_REG[col]] if col in _RAO_REG else _deg(g, "HCM2")  # HM2 = sum (d_u d_v)^2
        assert ours == printed, (comp["name"], col, ours, printed)


def test_rao_2024_rows_with_odd_first_zagreb_are_impossible():
    # M1 = sum_v d_v^2 has the parity of sum_v d_v = 2m, so it is even for every graph
    odd = sorted(c["name"] for c in RAO["compounds_table_3"] if c["values"]["M1"] % 2)
    assert odd == ["Amathamide C", "Amathaspiramide D"]
    for c in RAO["compounds_table_3"]:
        if "smiles" in c:
            assert compute_indices(_graph(c["smiles"]), ["M1"])[0]["M1"] % 2 == 0


def _cycle_edges(vs):
    return [(vs[i], vs[(i + 1) % len(vs)]) for i in range(len(vs))]


def _pendants(edges, n0, at, k):
    for _ in range(k):
        edges.append((at, n0))
        n0 += 1
    return n0


def _rao_graph(theorem, n, k=0, m=0, sizes=()):
    if theorem in ("2.1", "2.2", "2.3", "2.4"):
        edges, nxt = _cycle_edges(list(range(n))), n
        for v in range(n):
            if theorem == "2.1":
                nxt = _pendants(edges, nxt, v, 1)
            elif theorem == "2.2" and v == 0:
                nxt = _pendants(edges, nxt, v, k)
            elif theorem == "2.3":
                nxt = _pendants(edges, nxt, v, k)
            elif theorem == "2.4":
                nxt = _pendants(edges, nxt, v, k if v == 0 else 1)
        return MolGraph.from_edges(nxt, edges)
    if theorem in ("2.5", "2.6", "2.7"):
        sizes = {"2.5": (n, n), "2.6": (n, m), "2.7": sizes}[theorem]
        edges, nxt = [], 1
        for s in sizes:  # every cycle passes through vertex 0
            vs = [0] + list(range(nxt, nxt + s - 1))
            nxt += s - 1
            edges += _cycle_edges(vs)
        return MolGraph.from_edges(nxt, edges)
    raise ValueError(theorem)


def _rao_closed_forms(theorem, n, k=0, m=0, sizes=()):
    a = sum(sizes)
    K = len(sizes)
    return {
        "2.1": dict(CM1=6 * n, CM2=4 * n, HCM1=20 * n, HCM2=10 * n, M1=10 * n, M2=12 * n),
        "2.2": dict(CM1=2 * k * n + 2 * n + k * k + k, CM2=n * k * k + 2 * k * n + n - k * k,
                    HCM1=4 * n * k * k + 8 * k * n + 4 * n + k ** 3 + k,
                    HCM2=(n - 2) * (k + 1) ** 4 + k * (k + 2) ** 2 + 2 * (k + 1) ** 2,
                    M1=4 * n + k * k + 5 * k, M2=4 * n + 2 * k * k + 6 * k),
        "2.3": dict(CM1=n * (k + 1) * (k + 2), CM2=n * (k + 1) ** 2, HCM1=n * (k + 4) * (k + 1) ** 2,
                    HCM2=n * (k + 1) * (k * k + 3 * k + 1), M1=n * k * k + 5 * k * n + 4 * n,
                    M2=2 * n * k * k + 6 * k * n + 4 * n),
        "2.4": dict(CM1=2 * k * n + 2 * n + 3 * k * k - k, CM2=n * k * k + 2 * k * n + k ** 3 - 2 * k * k + 2 * k,
                    HCM1=4 * n * k * k + 8 * k * n + 4 * n + 5 * k ** 3 - 4 * k * k + 5 * k - 2,
                    HCM2=n * k ** 4 + 4 * n * k ** 3 + 4 * n * k * k + k ** 5 - 3 * k ** 4 - 3 * k ** 3 + 2 * k * k + 4 * k,
                    M1=4 * n + k * k + 11 * k - 6, M2=3 * n + k * k + 17 * k - 9),
        "2.5": dict(CM1=12 * n - 8, CM2=18 * n - 24, HCM1=72 * n - 80, HCM2=162 * n - 288, M1=8 * n + 8, M2=8 * n + 16),
        "2.6": dict(CM1=6 * m + 6 * n - 8, CM2=9 * m + 9 * n - 24, HCM1=36 * m + 36 * n - 80,
                    HCM2=81 * m + 81 * n - 288, M1=4 * m + 4 * n + 8, M2=4 * m + 4 * n + 16),
        "2.7": dict(CM1=(4 * K - 2) * (a - 2 * K) + 4 * K * K, CM2=(a - 2 * K) * (2 * K - 1) ** 2 + K * (4 * K - 2),
                    HCM1=(a - 2 * K) * (4 * K - 2) ** 2 + 8 * K ** 3,
                    HCM2=(a - 2 * K) * (2 * K - 1) ** 4 + 2 * K * (2 * K - 1) ** 2,
                    M1=4 * K * K - 4 * K + 4 * a, M2=8 * K * K - 8 * K + 4 * a),
    }[theorem]


_RAO_CASES = (
    [("2.1", dict(n=n)) for n in (3, 4, 6, 9)]
    + [("2.2", dict(n=n, k=k)) for n in (3, 5, 8) for k in (1, 2, 4)]
    + [("2.3", dict(n=n, k=k)) for n in (3, 5, 8) for k in (1, 2, 4)]
    + [("2.5", dict(n=n)) for n in (3, 5, 8)]
    + [("2.6", dict(n=n, m=m)) for n, m in ((3, 4), (5, 8), (6, 3))]
    + [("2.7", dict(n=0, sizes=s)) for s in ((3, 4, 5), (4, 4, 4, 4), (3, 6), (5, 5, 6, 7, 3))]
)
_RAO_24_CASES = [dict(n=n, k=k) for n in (3, 5, 8) for k in (2, 3, 4)]


def _rao_ours(g) -> dict:
    vals, reasons = compute_indices(g, ["R1_M1", "R1_M2", "R1_HM", "HCM2", "M1", "M2"])
    assert not reasons
    return dict(CM1=vals["R1_M1"], CM2=vals["R1_M2"], HCM1=vals["R1_HM"], HCM2=vals["HCM2"], M1=vals["M1"], M2=vals["M2"])


@pytest.mark.parametrize("theorem,kw", _RAO_CASES, ids=lambda x: str(x))
def test_rao_2024_cactus_theorems_on_explicit_graphs(theorem, kw):
    """Theorems 2.1, 2.3, 2.5, 2.6, 2.7 reproduce exactly; Theorem 2.2 reproduces except for its M2
    (paper error: it prints 4n + 2k^2 + 6k; the graph it defines, a cycle with k pendant edges at one
    vertex, has M2 = 4n + k^2 + 6k, as its own M1 = 4n + k^2 + 5k already implies a hub of degree k + 2)."""
    ours = _rao_ours(_rao_graph(theorem, **kw))
    expected = _rao_closed_forms(theorem, **kw)
    if theorem == "2.2":
        n, k = kw["n"], kw["k"]
        assert ours["M2"] == 4 * n + k * k + 6 * k != expected["M2"]
        expected["M2"] = ours["M2"]
    assert ours == expected


@pytest.mark.parametrize("kw", _RAO_24_CASES, ids=lambda x: str(x))
def test_rao_2024_theorem_2_4_fits_no_single_graph(kw):
    """Theorem 2.4 ("k pendant edges at one vertex of C_n and n - 1 pendant edges at the other
    vertices") is inconsistent with its own formulas: on the graph it defines none of the six indices
    matches (except by coincidence at n = 3, k = 3); its M1 and CM1 are those of a graph with one
    pendant edge at only k - 1 other vertices, but its M2, CM2, HCM1 and HCM2 are not those of that
    graph either. Not a TopoQSPR discrepancy: the index values on explicit graphs are unambiguous."""
    n, k = kw["n"], kw["k"]
    expected = _rao_closed_forms("2.4", n=n, k=k)
    defined = _rao_ours(_rao_graph("2.4", n=n, k=k))
    if (n, k) == (3, 3):
        assert defined == expected
    else:
        assert all(defined[key] != expected[key] for key in defined)
    if k - 1 <= n - 1:
        edges, nxt = _cycle_edges(list(range(n))), n
        nxt = _pendants(edges, nxt, 0, k)
        for v in range(1, k):
            nxt = _pendants(edges, nxt, v, 1)
        alt = _rao_ours(MolGraph.from_edges(nxt, edges))
        assert alt["M1"] == expected["M1"] and alt["CM1"] == expected["CM1"]
        if (n, k) != (3, 3):
            assert alt != expected  # the products (M2, CM2) or the hyper indices (HCM1, HCM2) differ


# ---------------------------------------------------------------------------
# Abubakar et al., Heliyon 10 (2024) e28260 -- antituberculosis drugs, neighbourhood degree sum
# ---------------------------------------------------------------------------
_ABU_REG = {"NGA": "S_GA", "NABC": "S_ABC", "NR": "S_R", "NRR": "S_RR", "NHM": "S_H"}
_ABU_PHI = {"NM1": "M1", "NM2": "M2", "NHP": "HM"}
# drugs whose printed NM1 equals 2 M2 of the real molecular graph (so the paper used the same graph)
_ABU_SAME_GRAPH = ["Ethambutol", "Ethionamide", "Linezolid", "Pyrazinamide", "Terizidone"]
_ABU_PAPER_ERRORS = {("Ethambutol", "NABC"): 8.0289, ("Linezolid", "NHP"): 3336}  # printed 9.91 and 3226
_ABU_BY_NAME = {d["name"]: d for d in ABUBAKAR["drugs_table_3"]}


def _abubakar_ours(g):
    vals, reasons = compute_indices(g, list(_ABU_REG.values()))
    assert not reasons
    out = {col: vals[name] for col, name in _ABU_REG.items()}
    out.update({col: _nbr(g, key) for col, key in _ABU_PHI.items()})
    return out


def test_abubakar_2024_terizidone_partition_is_table_2():
    g = _graph(_ABU_BY_NAME["Terizidone"]["smiles"])
    expected = Counter({_parse_pair(k): v for k, v in ABUBAKAR["terizidone_partition_table_2"].items()})
    assert _pairs(g, GraphContext(g).nbr_sum) == expected


@pytest.mark.parametrize("name", _ABU_SAME_GRAPH)
def test_abubakar_2024_table_3_full_pipeline(name):
    d = _ABU_BY_NAME[name]
    ours = _abubakar_ours(_graph(d["smiles"]))
    for col, printed in d["values"].items():
        _check_printed(name, col, ours[col], printed, _ABU_PAPER_ERRORS)


def test_abubakar_2024_first_neighbourhood_zagreb_identifies_the_graph():
    """NM1 = sum_uv (S_u + S_v) = sum_v d_v S_v = 2 M2 for every graph, so NM1 is even and fixes M2."""
    odd, other_graph = [], []
    for d in ABUBAKAR["drugs_table_3"]:
        g = _graph(d["smiles"])
        nm1 = _nbr(g, "M1")
        assert nm1 == 2 * compute_indices(g, ["M2"])[0]["M2"]
        printed = d["values"]["NM1"]
        if printed % 2:
            odd.append(d["name"])
        elif printed != nm1:
            other_graph.append((d["name"], printed, nm1))
        else:
            assert d["name"] in _ABU_SAME_GRAPH
    assert odd == ["Amikacin", "Bedaquiline", "Clofazimine", "Delamanid", "Isoniazid", "Levofloxacin",
                   "4-Aminosalicylic acid"]
    assert other_graph == [("Imipenem", 286, 246), ("Moxifloxacin", 424, 414), ("Rifampin", 742, 760)]


# ---------------------------------------------------------------------------
# Almulla et al., Sci. Rep. 16 (2026) 25695 -- antifungal drugs, neighbourhood degree sum
# ---------------------------------------------------------------------------
_ALM_INCONSISTENT = {"Marbofloxacin", "Terbinafine", "Miconazole", "Voriconazole"}
_ALM_OK = [d for d in ALMULLA["drugs"] if d["name"] not in _ALM_INCONSISTENT]


def _alm_partition(d) -> Counter:
    return Counter({_parse_pair(k): v for k, v in d["partition_table_2"].items()})


@pytest.mark.parametrize("drug", _ALM_OK, ids=lambda d: d["name"])
def test_almulla_2026_tables_2_and_4_full_pipeline(drug):
    g = _graph(drug["smiles"])
    assert _pairs(g, GraphContext(g).nbr_sum) == _alm_partition(drug)
    v = drug["values_table_4"]
    assert _nbr(g, "M1") == v["NM1"] and _nbr(g, "M2") == v["NM2"] and _nbr(g, "F") == v["NF"]
    s_abc = compute_indices(g, ["S_ABC"])[0]["S_ABC"]
    assert _matches_printed(s_abc, v["NABC"]), (drug["name"], s_abc, v["NABC"])


@pytest.mark.parametrize("drug", [d for d in ALMULLA["drugs"] if d["name"] in _ALM_INCONSISTENT],
                         ids=lambda d: d["name"])
def test_almulla_2026_rows_whose_partition_is_not_the_molecule(drug):
    """The printed values follow from the printed partition, but the partition is not that of the molecule."""
    g = _graph(drug["smiles"])
    part = _alm_partition(drug)
    v = drug["values_table_4"]
    assert _from_partition(part, "M1") == v["NM1"] and _from_partition(part, "M2") == v["NM2"]
    assert _from_partition(part, "ABC") == pytest.approx(v["NABC"], abs=5e-5)
    assert _pairs(g, GraphContext(g).nbr_sum) != part
    n_edges = {"Marbofloxacin": (28, 29), "Terbinafine": (21, 23), "Miconazole": (27, 27), "Voriconazole": (28, 27)}
    assert (sum(part.values()), len(g.edges)) == n_edges[drug["name"]]
    if drug["name"] in ("Terbinafine", "Miconazole"):
        assert v["NM1"] % 2 == 1  # impossible for any graph (NM1 = 2 M2)


# ---------------------------------------------------------------------------
# Mondal, De, Pal, Eur. Phys. J. Plus 136 (2021) 303 -- Cayley tree dendrimers
# ---------------------------------------------------------------------------
def _cayley(s: int, t: int) -> MolGraph:
    """T_C(s, t): a centre with s children; every later vertex gets s - 1 children; t levels."""
    edges, nxt, frontier = [], 1, [0]
    for level in range(1, t + 1):
        new = []
        for v in frontier:
            for _ in range(s if level == 1 else s - 1):
                edges.append((v, nxt))
                new.append(nxt)
                nxt += 1
        frontier = new
    return MolGraph.from_edges(nxt, edges)


@pytest.mark.parametrize("s", [3, 4, 5])
@pytest.mark.parametrize("t", [3, 4, 5])
def test_mondal_2021_cayley_tree_closed_forms(s, t):
    g = _cayley(s, t)
    S = GraphContext(g).nbr_sum
    sig = sum((s - 1) ** (i - 1) for i in range(1, t + 1))
    sig2 = sum((s - 1) ** (i - 1) for i in range(1, t - 1))
    a = (s - 1) ** (t - 2)
    assert g.n == s * sig + 1 and len(g.edges) == s * sig
    # Corollaries 9, 10 (vertex sums), Theorems 6, 7, Corollaries 14, 15, 16
    assert sum(x * x for x in S) == s * a * (s ** 3 + 3 * s * s - 4 * s + 1) + s ** 5 * sig2 + s ** 4
    assert sum(x ** 3 for x in S) == s * a * (s ** 4 + 7 * s ** 3 - 12 * s * s + 6 * s - 1) + s ** 7 * sig2 + s ** 6
    assert _nbr(g, "F") == s * s * a * (-2 * s ** 4 + s ** 3 + 5 * s * s - 5 * s + 1) + 2 * s ** 5 * sig
    assert _nbr(g, "M2") == s * s * a * (-(s ** 4) + 4 * s * s - 4 * s + 1) + s ** 5 * sig
    assert _nbr(g, "HM") == s * s * a * (-4 * s ** 4 + s ** 3 + 13 * s * s - 13 * s + 3) + 4 * s ** 5 * sig
    assert _nbr(g, "M1") == 2 * s ** 3 * sig - 2 * s * s * (s - 1) ** t
    nd2 = s * (s - 1) ** (t - 1) / math.sqrt(3 * s - 1) + s * a / math.sqrt(s * s + 2 * s - 1) + (sig - s * a) / math.sqrt(2)
    assert compute_indices(g, ["S_SC"])[0]["S_SC"] == pytest.approx(nd2, rel=1e-12)


# ---------------------------------------------------------------------------
# Ravi & Desikan, Sci. Rep. 13 (2023) 3239 -- reduced reverse degree (different convention)
# ---------------------------------------------------------------------------
def test_ravi_2023_reduced_reverse_degree_is_kulli_plus_one():
    # RR(v) = Delta - d + 2 = c_v + 1: a different vertex weight from every TopoQSPR reverse-degree index
    for delta in range(2, 7):
        for d in range(1, delta + 1):
            assert delta - d + 2 == modified_reverse_degree(d, 1, delta) + 1


@pytest.mark.parametrize("row", RAVI2023["benzenoids_table_2_RRR_minus1"], ids=lambda r: r["name"])
def test_ravi_2023_printed_benzenoid_values_are_plain_degree_values(row):
    """Table 2 prints RRR_{-1} = sum 1/(RR_u RR_v), but the numbers are sum 1/(d_u d_v) (paper error)."""
    g = _graph(row["smiles"])
    plain = _deg(g, "M2") and sum(1 / (g.degrees[u] * g.degrees[v]) for u, v in g.edges)
    reduced = sum(1 / ((g.max_degree - g.degrees[u] + 2) * (g.max_degree - g.degrees[v] + 2)) for u, v in g.edges)
    assert _matches_printed(plain, row["value"])
    if row["name"] == "benzene":
        assert reduced == plain  # Delta = 2: RR = d
    else:
        assert not _matches_printed(reduced, row["value"])


def test_ravi_2023_hap_closed_forms_from_table_7():
    """Reduced-reverse edge partition of HA-paclitaxel (Table 7) with TopoQSPR's edge functions."""
    cf = RAVI2023["hap_closed_forms"]
    funcs = {k: PHI[k][0] for k in ("M1", "M2", "HM", "HCM2", "F", "ABC", "GA", "RR", "R")}
    funcs["Rm1"] = lambda a, b: 1 / (a * b)
    part_a = {_parse_pair(k): a for k, (a, b) in RAVI2023["hap_partition_table_7"].items()}
    part_b = {_parse_pair(k): b for k, (a, b) in RAVI2023["hap_partition_table_7"].items()}
    assert sum(part_a.values()) == 96 and sum(part_b.values()) == 0  # "|E((HAP)_n)| = 96n"
    for key, (slope, icpt) in cf.items():
        ours_slope, ours_icpt = _from_partition(part_a, funcs[key]), _from_partition(part_b, funcs[key])
        if key in ("GA", "R"):
            # printed 93.9919881117 n + 1.9836 and 26.1823 n - 0.0075; the paper's own partition gives
            # 94.4586 n + 0.0041 and 28.8590 n - 0.1484 (slopes differ by 0.47 and 2.68 per unit)
            assert abs(ours_slope - slope) > 0.4, key
        else:
            assert ours_slope == pytest.approx(slope, abs=2e-3) and ours_icpt == pytest.approx(icpt, abs=2e-3), key
    one = {_parse_pair(k): a + b for k, (a, b) in RAVI2023["hap_partition_table_7"].items()}
    assert _from_partition(one, funcs["GA"]) == pytest.approx(94.4626, abs=5e-5)
    assert _from_partition(one, funcs["R"]) == pytest.approx(28.7106, abs=5e-5)


# ---------------------------------------------------------------------------
# Ravi & Chidambaram, Sci. Rep. 15 (2025) 16941 -- Parkinson drugs, open neighbourhood sum
# ---------------------------------------------------------------------------
def test_ravi_2025_procyclidine_open_neighbourhood_partition_is_table_2():
    g = _graph(RAVI2025["procyclidine"]["smiles"])
    expected = Counter({_parse_pair(k): v for k, v in RAVI2025["procyclidine_open_partition_table_2"].items()})
    assert _pairs(g, GraphContext(g).nbr_sum) == expected


def test_fixtures_carry_doi_and_source():
    for f in sorted(FIX.glob("*.json")):
        d = json.loads(f.read_text())
        assert d["doi"].startswith("10.") and d["source"] and d["convention"], f.name
