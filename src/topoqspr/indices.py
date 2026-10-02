r"""Verified graph-theoretic topological indices.

Every index is an entry of :data:`REGISTRY` with

* ``name``    - the column name used throughout TopoQSPR,
* ``family``  - one of :data:`FAMILIES`,
* ``latex``   - the defining formula,
* ``ref``     - the BibTeX key (``docs/index_refs.bib``) of the source that
  defines it,
* ``func``    - a function of a :class:`GraphContext`.

:func:`compute_indices` evaluates the registry on a
:class:`~topoqspr.graph.MolGraph` and returns ``(values, reasons)``. An index
that is mathematically undefined for the graph (division by zero, a
disconnected graph for a distance index, ``k > \Delta`` for a modified
reverse degree) is ``NaN`` in ``values`` and its reason is in ``reasons``.
Nothing is ever replaced by 0.

Notation (G simple, hydrogen-suppressed, ``n`` vertices, ``m`` edges)
-------------------------------------------------------------------
``d_u``   degree of vertex u; ``\Delta`` the maximum degree.
``S_u``   neighbourhood degree sum, ``S_u = \sum_{w \sim u} d_w``.
``c_u``   Kulli reverse degree, ``c_u = \Delta - d_u + 1`` (Kulli 2018).
``M_kR(d_u)`` modified (generalised) reverse degree (Kalaam et al. 2024)::

        M_kR(d) = \Delta + k - d              if k <= d
                = (\Delta + k - d) mod \Delta   if k >  d

``D(u, v)`` shortest-path distance; ``D_u = \sum_v D(u, v)``.
For an edge ``e = uv``: ``n_u(e)`` = #vertices closer to u than to v,
``m_u(e)`` = #edges f with ``D(f, u) < D(f, v)`` where
``D(f, x) = min(D(a, x), D(b, x))`` for ``f = ab``.

Edge-additive indices are written as ``\sum_{uv \in E} \varphi(x_u, x_v)``
where the vertex weight ``x`` is ``d``, ``S`` or ``M_kR(d)`` depending on
the family; the same edge function ``\varphi`` is reused across families, as
in Godlin & Radha (2026) and Arockiaraj et al. (2025), Eqs. (2)-(4).

Integer-valued indices are returned as exact Python ``int`` (Hosoya Z in
particular can exceed 2**53); the others as ``float``.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass
from functools import cached_property
from typing import Callable, Iterable, Optional, Sequence, Union

from .graph import MolGraph

__all__ = [
    "FAMILIES",
    "IndexSpec",
    "REGISTRY",
    "INDEX_NAMES",
    "UndefinedIndex",
    "GraphContext",
    "compute_indices",
    "hosoya_z",
    "distance_matrix",
    "modified_reverse_degree",
    "RK_VALUES",
    "PHI",
]

Number = Union[int, float]

FAMILIES: tuple[str, ...] = (
    "degree",
    "reverse_degree",
    "neighbourhood_degree_sum",
    "distance",
    "matching",
    "entropy",
)

#: Values of k for which the modified reverse degree indices are registered
#: (the range used by Kalaam et al. 2024 and Godlin & Radha 2026).
RK_VALUES: tuple[int, ...] = (1, 2, 3)


class UndefinedIndex(Exception):
    """Raised by an index function when the index is undefined for a graph."""


# ---------------------------------------------------------------------------
# Graph context (lazily cached derived quantities)
# ---------------------------------------------------------------------------
def distance_matrix(graph: MolGraph) -> list[list[Optional[int]]]:
    """All-pairs shortest-path distances by BFS (``None`` = unreachable)."""
    n, adj = graph.n, graph.adjacency
    out: list[list[Optional[int]]] = []
    for s in range(n):
        dist: list[Optional[int]] = [None] * n
        dist[s] = 0
        q = deque([s])
        while q:
            v = q.popleft()
            dv = dist[v] + 1  # type: ignore[operator]
            for u in adj[v]:
                if dist[u] is None:
                    dist[u] = dv
                    q.append(u)
        out.append(dist)
    return out


def modified_reverse_degree(d: int, k: int, delta: int) -> int:
    r"""Modified reverse degree ``M_kR(d)`` (Kalaam, Greeni, Arockiaraj 2024).

    ``M_kR(d) = \Delta + k - d`` if ``k <= d``, else ``(\Delta + k - d) mod \Delta``.
    The ``mod`` applies to the whole sum, as written with explicit parentheses
    by Godlin & Radha (2026), Eq. (1), and as their worked tables confirm
    (Delta = 4: k = 2 maps d = 1, 2, 3, 4 to 1, 4, 3, 2).

    Raises :class:`UndefinedIndex` unless ``1 <= k <= Delta``. Inside that
    range the map is a permutation of ``{1, ..., Delta}`` (never 0). The
    sources only use ``k <= Delta`` ("the maximum degrees ... are 3 or 4;
    therefore ... l is limited from 1 to 3", Godlin & Radha 2026); for
    ``k > Delta`` the rule sends every vertex with ``d = k - Delta`` to 0,
    which is not a reverse degree, so we refuse rather than guess.
    """
    if delta < 1:
        raise UndefinedIndex("modified reverse degree needs Delta(G) >= 1")
    if not 1 <= k <= delta:
        raise UndefinedIndex(
            f"modified reverse degree with k={k} is outside 1..Delta(G)={delta}"
        )
    if k <= d:
        return delta + k - d
    return (delta + k - d) % delta


class GraphContext:
    """A graph plus lazily computed invariants shared by many indices."""

    def __init__(self, graph: MolGraph):
        self.g = graph

    @cached_property
    def deg(self) -> tuple[int, ...]:
        return self.g.degrees

    @cached_property
    def nbr_sum(self) -> tuple[int, ...]:
        """S_u = sum of the degrees of the neighbours of u."""
        d, adj = self.g.degrees, self.g.adjacency
        return tuple(sum(d[w] for w in adj[u]) for u in range(self.g.n))

    def rk(self, k: int) -> tuple[int, ...]:
        delta = self.g.max_degree
        return tuple(modified_reverse_degree(d, k, delta) for d in self.g.degrees)

    @cached_property
    def dist(self) -> list[list[int]]:
        """Distance matrix; raises UndefinedIndex if G is disconnected/empty."""
        if self.g.n == 0:
            raise UndefinedIndex("empty graph")
        if not self.g.is_connected:
            raise UndefinedIndex(
                f"distance-based index undefined: graph has {self.g.n_components} components"
            )
        return distance_matrix(self.g)  # type: ignore[return-value]

    @cached_property
    def dist_sums(self) -> list[int]:
        return [sum(row) for row in self.dist]

    @cached_property
    def edge_split(self) -> list[tuple[int, int, int, int]]:
        """For every edge e = uv (in graph order): (n_u, n_v, m_u, m_v)."""
        D = self.dist
        n, edges = self.g.n, self.g.edges
        out = []
        for u, v in edges:
            du, dv = D[u], D[v]
            n_u = sum(1 for x in range(n) if du[x] < dv[x])
            n_v = sum(1 for x in range(n) if dv[x] < du[x])
            m_u = m_v = 0
            for a, b in edges:
                fu = du[a] if du[a] < du[b] else du[b]
                fv = dv[a] if dv[a] < dv[b] else dv[b]
                if fu < fv:
                    m_u += 1
                elif fv < fu:
                    m_v += 1
            out.append((n_u, n_v, m_u, m_v))
        return out


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class IndexSpec:
    name: str
    family: str
    latex: str
    ref: str
    func: Callable[[GraphContext], Number]
    description: str = ""


def _edge_sum(
    phi: Callable[[int, int], Number], weights: Sequence[int], edges: Iterable[tuple[int, int]]
) -> Number:
    total: Number = 0
    for u, v in edges:
        a, b = weights[u], weights[v]
        try:
            total += phi(a, b)
        except ZeroDivisionError:
            raise UndefinedIndex(
                f"division by zero on edge ({u},{v}) with end weights ({a},{b})"
            ) from None
        except ValueError as exc:  # sqrt of a negative number
            raise UndefinedIndex(f"{exc} on edge ({u},{v}) with end weights ({a},{b})") from None
    return total


# Edge functions phi(a, b). Exact integer arithmetic where the value is integral.
_sqrt = math.sqrt
PHI: dict[str, tuple[Callable[[int, int], Number], str]] = {
    "M1": (lambda a, b: a + b, r"x_u + x_v"),
    "M2": (lambda a, b: a * b, r"x_u x_v"),
    "F": (lambda a, b: a * a + b * b, r"x_u^2 + x_v^2"),
    "HM": (lambda a, b: (a + b) ** 2, r"(x_u + x_v)^2"),
    "R": (lambda a, b: 1.0 / _sqrt(a * b), r"\frac{1}{\sqrt{x_u x_v}}"),
    "RR": (lambda a, b: _sqrt(a * b), r"\sqrt{x_u x_v}"),
    "SC": (lambda a, b: 1.0 / _sqrt(a + b), r"\frac{1}{\sqrt{x_u + x_v}}"),
    "H": (lambda a, b: 2.0 / (a + b), r"\frac{2}{x_u + x_v}"),
    "GA": (lambda a, b: 2.0 * _sqrt(a * b) / (a + b), r"\frac{2\sqrt{x_u x_v}}{x_u + x_v}"),
    "ABC": (lambda a, b: _sqrt((a + b - 2) / (a * b)), r"\sqrt{\frac{x_u + x_v - 2}{x_u x_v}}"),
    "AZI": (lambda a, b: (a * b / (a + b - 2)) ** 3, r"\left(\frac{x_u x_v}{x_u + x_v - 2}\right)^3"),
    "SO": (lambda a, b: _sqrt(a * a + b * b), r"\sqrt{x_u^2 + x_v^2}"),
    "ISI": (lambda a, b: a * b / (a + b), r"\frac{x_u x_v}{x_u + x_v}"),
    "SDD": (lambda a, b: (a * a + b * b) / (a * b), r"\frac{x_u^2 + x_v^2}{x_u x_v}"),
    "sigma": (lambda a, b: (a - b) ** 2, r"(x_u - x_v)^2"),
    "ABS": (lambda a, b: _sqrt((a + b - 2) / (a + b)), r"\sqrt{\frac{x_u + x_v - 2}{x_u + x_v}}"),
    "ReZ3": (lambda a, b: a * b * (a + b), r"x_u x_v (x_u + x_v)"),
    "AG": (lambda a, b: (a + b) / (2.0 * _sqrt(a * b)), r"\frac{x_u + x_v}{2\sqrt{x_u x_v}}"),
    "HG": (lambda a, b: 2.0 / (_sqrt(a * b) * (a + b)), r"\frac{2}{\sqrt{x_u x_v}\,(x_u + x_v)}"),
    "GBM": (lambda a, b: _sqrt(a * b) / (a + b + a * b), r"\frac{\sqrt{x_u x_v}}{x_u + x_v + x_u x_v}"),
    "BMG": (lambda a, b: (a + b + a * b) / _sqrt(a * b), r"\frac{x_u + x_v + x_u x_v}{\sqrt{x_u x_v}}"),
    "BM": (lambda a, b: a + b + a * b, r"x_u + x_v + x_u x_v"),
    "TM": (lambda a, b: a * a + b * b + a * b, r"x_u^2 + x_v^2 + x_u x_v"),
    # BMH/TMH are integers or half-integers; returned as float (exact below 2**52).
    "BMH": (lambda a, b: (a + b + a * b) * (a + b) / 2, r"\frac{(x_u + x_v + x_u x_v)(x_u + x_v)}{2}"),
    "TMH": (lambda a, b: (a * a + b * b + a * b) * (a + b) / 2, r"\frac{(x_u^2 + x_v^2 + x_u x_v)(x_u + x_v)}{2}"),
    "HCM2": (lambda a, b: (a * b) ** 2, r"(x_u x_v)^2"),
}

_INT_PHI = {"M1", "M2", "F", "HM", "sigma", "ReZ3", "BM", "TM", "HCM2"}


_BRACKET = {"M1", "F", "BM", "TM"}  # phi with a top-level sum: bracket it under the sigma


def _phi_latex(key: str, weight: str) -> str:
    body = PHI[key][1].replace("x_u", weight + "_u").replace("x_v", weight + "_v")
    if key in _BRACKET:
        body = r"\left[" + body + r"\right]"
    return r"\sum_{uv \in E(G)} " + body


# --- degree family -----------------------------------------------------------
_DEGREE: list[tuple[str, str, str]] = [
    # (name, phi key, ref)
    ("M1", "M1", "gutman1972"),
    ("M2", "M2", "gutman1975"),
    ("F", "F", "furtula2015"),
    ("HM", "HM", "shirdel2013"),
    ("R", "R", "randic1975"),
    ("RR", "RR", "gutman2014rr"),
    ("SC", "SC", "zhou2009"),
    ("H", "H", "fajtlowicz1987"),
    ("GA", "GA", "vukicevic2009"),
    ("ABC", "ABC", "estrada1998"),
    ("AZI", "AZI", "furtula2010"),
    ("SO", "SO", "gutman2021"),
    ("ISI", "ISI", "vukicevic2010"),
    ("SDD", "SDD", "vukicevic2010"),
    ("sigma", "sigma", "gutman2018sigma"),
    ("ABS", "ABS", "ali2022"),
    ("ReZ3", "ReZ3", "ranjini2013"),
    ("AG", "AG", "shegehall2015"),
    ("HG", "HG", "arockiaraj2023sar"),
    ("GBM", "GBM", "arockiaraj2023sar"),
    ("BMG", "BMG", "arockiaraj2023sar"),
    ("BMH", "BMH", "arockiaraj2023sar"),
    ("TMH", "TMH", "arockiaraj2023sar"),
]

# --- neighbourhood-degree-sum family -----------------------------------------
# (v1 name, phi key, ref, literature name)
_NBR: list[tuple[str, str, str, str]] = [
    ("S_ABC", "ABC", "ghorbani2010", "ABC_4"),
    ("S_GA", "GA", "graovac2011", "GA_5"),
    ("S_RR", "RR", "mondal2021", "ND_1"),
    ("S_SC", "SC", "mondal2021", "ND_2"),
    ("S_ReZ3", "ReZ3", "mondal2021", "ND_3"),
    ("S_R", "R", "mondal2021", "ND_4"),
    ("S_SDD", "SDD", "mondal2021", "ND_5"),
    ("S_H", "H", "arockiaraj2025lung", "H^{nd}"),
    ("S_GBM", "GBM", "arockiaraj2025lung", "GBM^{nd}"),
    ("S_TMH", "TMH", "arockiaraj2025lung", "TMH^{nd}"),
]

# --- modified reverse degree family (k = 1..3) --------------------------------
# Kalaam et al. 2024 Eqs. (1)-(9); R, ABC, ISI from Godlin & Radha 2026; RR from Arockiaraj et al. 2025.
_RK: list[tuple[str, str]] = [
    ("M1", "kalaam2024"),
    ("M2", "kalaam2024"),
    ("F", "kalaam2024"),
    ("HM", "kalaam2024"),  # "HZ" (hyper-Zagreb) in the source
    ("ReZ3", "kalaam2024"),
    ("BM", "kalaam2024"),
    ("TM", "kalaam2024"),
    ("BMH", "kalaam2024"),
    ("TMH", "kalaam2024"),
    ("R", "godlin2026"),
    ("ABC", "godlin2026"),
    ("ISI", "godlin2026"),  # "ReZ_2" (redefined second Zagreb) in the source
    ("RR", "arockiaraj2025lung"),  # "G" (geometric) in the source
]


def _as_int_if(key: str, value: Number) -> Number:
    return int(value) if key in _INT_PHI else float(value)


def _deg_func(key: str) -> Callable[[GraphContext], Number]:
    phi = PHI[key][0]

    def f(ctx: GraphContext) -> Number:
        return _as_int_if(key, _edge_sum(phi, ctx.deg, ctx.g.edges))

    return f


def _nbr_func(key: str) -> Callable[[GraphContext], Number]:
    phi = PHI[key][0]

    def f(ctx: GraphContext) -> Number:
        return _as_int_if(key, _edge_sum(phi, ctx.nbr_sum, ctx.g.edges))

    return f


def _rk_func(key: str, k: int) -> Callable[[GraphContext], Number]:
    phi = PHI[key][0]

    def f(ctx: GraphContext) -> Number:
        return _as_int_if(key, _edge_sum(phi, ctx.rk(k), ctx.g.edges))

    return f


# --- distance family ---------------------------------------------------------
def _wiener(ctx: GraphContext) -> int:
    """W = sum over unordered vertex pairs of D(u, v)."""
    return sum(ctx.dist_sums) // 2


def _balaban_j(ctx: GraphContext) -> float:
    """J = m/(mu+1) * sum_{uv in E} (D_u D_v)^(-1/2), mu = m - n + 1."""
    g = ctx.g
    if g.m == 0:
        raise UndefinedIndex("Balaban J undefined for a graph without edges")
    mu = g.m - g.n + 1
    s = ctx.dist_sums
    return g.m / (mu + 1) * sum(1.0 / math.sqrt(s[u] * s[v]) for u, v in g.edges)


def _szeged(ctx: GraphContext) -> int:
    return sum(nu * nv for nu, nv, _, _ in ctx.edge_split)


def _edge_szeged(ctx: GraphContext) -> int:
    return sum(mu * mv for _, _, mu, mv in ctx.edge_split)


def _pi_edge(ctx: GraphContext) -> int:
    return sum(mu + mv for _, _, mu, mv in ctx.edge_split)


def _pi_vertex(ctx: GraphContext) -> int:
    return sum(nu + nv for nu, nv, _, _ in ctx.edge_split)


def _mostar(ctx: GraphContext) -> int:
    return sum(abs(nu - nv) for nu, nv, _, _ in ctx.edge_split)


def _edge_mostar(ctx: GraphContext) -> int:
    return sum(abs(mu - mv) for _, _, mu, mv in ctx.edge_split)


def _edge_wiener(ctx: GraphContext) -> int:
    """W_e = W(L(G)); for e = ab != f = xy, D_L(e, f) = min D(endpoints) + 1."""
    D = ctx.dist
    edges = ctx.g.edges
    total = 0
    for i, (a, b) in enumerate(edges):
        Da, Db = D[a], D[b]
        for x, y in edges[i + 1:]:
            total += min(Da[x], Da[y], Db[x], Db[y]) + 1
    return total


def _degree_distance(ctx: GraphContext) -> int:
    """DD = sum_{u<v} (d_u + d_v) D(u, v) = sum_u d_u D_u."""
    return sum(d * s for d, s in zip(ctx.deg, ctx.dist_sums))


def _gutman(ctx: GraphContext) -> int:
    """Gut = sum_{u<v} d_u d_v D(u, v)."""
    D, d, n = ctx.dist, ctx.deg, ctx.g.n
    return sum(d[u] * d[v] * D[u][v] for u in range(n) for v in range(u + 1, n))


# --- Hosoya Z ----------------------------------------------------------------
def _elimination_order(graph: MolGraph) -> list[int]:
    """Vertex order with a small 'frontier' for the matching-count DP.

    Greedy: start each component at a pseudo-peripheral vertex (end of a
    double BFS sweep from its smallest vertex) and repeatedly take the
    unprocessed vertex, adjacent to the processed set, that leaves the fewest
    processed vertices with unprocessed neighbours (ties: more processed
    neighbours, then smaller index). Deterministic.
    """
    n, adj = graph.n, graph.adjacency

    def farthest(src: int, allowed: set[int]) -> int:
        dist = {src: 0}
        q = deque([src])
        last = src
        while q:
            v = q.popleft()
            last = v
            for u in adj[v]:
                if u in allowed and u not in dist:
                    dist[u] = dist[v] + 1
                    q.append(u)
        # among vertices at maximal distance take the smallest index
        far = max(dist.values())
        return min(v for v, d in dist.items() if d == far) if dist else last

    done = [False] * n
    remaining = [len(adj[v]) for v in range(n)]  # unprocessed neighbours
    order: list[int] = []
    frontier: set[int] = set()
    for start in range(n):
        if done[start]:
            continue
        # collect the component of `start`
        comp = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for u in adj[v]:
                if u not in comp:
                    comp.add(u)
                    q.append(u)
        root = farthest(farthest(start, comp), comp)
        candidates = {root}
        while candidates:
            best, best_key = -1, None
            for v in candidates:
                n_proc = len(adj[v]) - remaining[v]
                closes = sum(1 for u in adj[v] if done[u] and remaining[u] == 1)
                new_size = len(frontier) - closes + (1 if remaining[v] - 0 > 0 else 0)
                key = (new_size, -n_proc, v)
                if best_key is None or key < best_key:
                    best, best_key = v, key
            v = best
            candidates.discard(v)
            done[v] = True
            order.append(v)
            for u in adj[v]:
                remaining[u] -= 1
                if done[u]:
                    if remaining[u] == 0:
                        frontier.discard(u)
                else:
                    candidates.add(u)
            if remaining[v] > 0:
                frontier.add(v)
    return order


def hosoya_z(graph: MolGraph) -> int:
    """Hosoya index Z(G): the number of matchings of G, the empty one included.

    Z(G) = sum_k p(G, k), p(G, k) = number of k-element sets of pairwise
    disjoint edges, p(G, 0) = 1 (Hosoya 1971). Exact Python int.

    Algorithm: dynamic programming over a vertex order (a path-decomposition
    style "frontier" DP). After processing a prefix of the order, the state is
    the set of processed vertices that still have unprocessed neighbours and
    are not yet matched; each state carries the number of partial matchings
    leading to it. Adding vertex v either leaves v unmatched or matches it to
    one free processed neighbour; vertices whose neighbours are all processed
    leave the state. The cost is O(n * 2^w * Delta) for frontier width w,
    which is small for molecular graphs (w <= ~8 for drugs with <= 150 heavy
    atoms), so this is exact and fast. Verified against brute-force
    enumeration and the closed forms Z(P_n) = F_{n+1}, Z(C_n) = L_n,
    Z(K_{1,n-1}) = n.
    """
    n, adj = graph.n, graph.adjacency
    order = _elimination_order(graph)
    done = [False] * n
    remaining = [len(adj[v]) for v in range(n)]
    states: dict[int, int] = {0: 1}  # bitmask of free frontier vertices -> count
    for v in order:
        bit_v = 1 << v
        proc_nbrs = [u for u in adj[v] if done[u]]
        new: dict[int, int] = {}
        for free, cnt in states.items():
            key = free | bit_v  # v left unmatched (for now)
            new[key] = new.get(key, 0) + cnt
            for u in proc_nbrs:
                bu = 1 << u
                if free & bu:
                    key = free & ~bu  # match v with u
                    new[key] = new.get(key, 0) + cnt
        done[v] = True
        # vertices that leave the frontier: v (if it has no unprocessed
        # neighbour) and processed neighbours whose last neighbour was v
        drop = 0
        for u in proc_nbrs:
            remaining[u] -= 1
            if remaining[u] == 0:
                drop |= 1 << u
        remaining_v = sum(1 for u in adj[v] if not done[u])
        remaining[v] = remaining_v
        if remaining_v == 0:
            drop |= bit_v
        if drop:
            keep = ~drop
            merged: dict[int, int] = {}
            for free, cnt in new.items():
                key = free & keep
                merged[key] = merged.get(key, 0) + cnt
            new = merged
        states = new
    return sum(states.values())


# --- entropy -----------------------------------------------------------------
def _degree_entropy(ctx: GraphContext) -> float:
    """I_d(G) = -sum_v (d_v / 2m) log2(d_v / 2m) (Cao, Dehmer, Shi 2014, t = 1)."""
    two_m = 2 * ctx.g.m
    if two_m == 0:
        raise UndefinedIndex("degree entropy undefined for a graph without edges")
    h = 0.0
    for d in ctx.deg:
        if d:
            p = d / two_m
            h -= p * math.log2(p)
    return h


def _build_registry() -> tuple[IndexSpec, ...]:
    specs: list[IndexSpec] = []
    for name, key, ref in _DEGREE:
        specs.append(IndexSpec(name, "degree", _phi_latex(key, "d"), ref, _deg_func(key)))
    for k in RK_VALUES:
        for key, ref in _RK:
            specs.append(
                IndexSpec(
                    f"R{k}_{key}",
                    "reverse_degree",
                    _phi_latex(key, "r").replace("r_u", rf"M_{{{k}}}R(d_u)").replace(
                        "r_v", rf"M_{{{k}}}R(d_v)"
                    ),
                    ref,
                    _rk_func(key, k),
                    description=f"phi = {key} evaluated on modified reverse degrees, k = {k}",
                )
            )
    specs.append(
        IndexSpec(
            "HCM2",
            "reverse_degree",
            r"\sum_{uv \in E(G)} (c_u c_v)^2,\ c_u = \Delta - d_u + 1",
            "kulli2018",
            _rk_func("HCM2", 1),
            description="second reverse hyper-Zagreb index (c_u = M_1R(d_u))",
        )
    )
    for name, key, ref, lit in _NBR:
        specs.append(
            IndexSpec(
                name,
                "neighbourhood_degree_sum",
                _phi_latex(key, "S"),
                ref,
                _nbr_func(key),
                description=f"literature name {lit}",
            )
        )
    dist = [
        ("W", r"\sum_{\{u,v\} \subseteq V} D(u,v)", "wiener1947", _wiener),
        (
            "J",
            r"\frac{m}{\mu + 1} \sum_{uv \in E} (D_u D_v)^{-1/2},\ \mu = m - n + 1",
            "balaban1982",
            _balaban_j,
        ),
        ("Sz", r"\sum_{uv \in E} n_u(e) n_v(e)", "gutman1994sz", _szeged),
        ("Sz_e", r"\sum_{uv \in E} m_u(e) m_v(e)", "gutman2008sze", _edge_szeged),
        ("PI", r"\sum_{uv \in E} [m_u(e) + m_v(e)]", "khadikar2000", _pi_edge),
        ("PI_v", r"\sum_{uv \in E} [n_u(e) + n_v(e)]", "khalifeh2008pi", _pi_vertex),
        ("Mo", r"\sum_{uv \in E} |n_u(e) - n_v(e)|", "doslic2018", _mostar),
        ("Mo_e", r"\sum_{uv \in E} |m_u(e) - m_v(e)|", "arockiaraj2019", _edge_mostar),
        ("W_e", r"W(L(G)) = \sum_{\{e,f\} \subseteq E} D_{L(G)}(e,f)", "iranmanesh2009", _edge_wiener),
        ("DD", r"\sum_{\{u,v\} \subseteq V} (d_u + d_v) D(u,v)", "dobrynin1994", _degree_distance),
        ("Gut", r"\sum_{\{u,v\} \subseteq V} d_u d_v D(u,v)", "gutman1994schultz", _gutman),
    ]
    for name, latex, ref, fn in dist:
        specs.append(IndexSpec(name, "distance", latex, ref, fn))
    specs.append(
        IndexSpec(
            "Z",
            "matching",
            r"Z(G) = \sum_{k \ge 0} p(G, k)",
            "hosoya1971",
            lambda ctx: hosoya_z(ctx.g),
            description="number of k-matchings summed over k, empty matching included",
        )
    )
    specs.append(
        IndexSpec(
            "I_deg",
            "entropy",
            r"-\sum_{v \in V} \frac{d_v}{2m} \log_2 \frac{d_v}{2m}",
            "cao2014",
            _degree_entropy,
        )
    )
    names = [s.name for s in specs]
    assert len(names) == len(set(names)), "duplicate index names"
    return tuple(specs)


REGISTRY: tuple[IndexSpec, ...] = _build_registry()
INDEX_NAMES: tuple[str, ...] = tuple(s.name for s in REGISTRY)
_BY_NAME = {s.name: s for s in REGISTRY}


def compute_indices(
    graph: MolGraph, names: Optional[Iterable[str]] = None
) -> tuple[dict[str, Number], dict[str, str]]:
    """Evaluate registry indices on ``graph``.

    Returns ``(values, reasons)``. ``values`` has one entry per requested
    index, in registry order; an undefined index is ``float('nan')`` and
    ``reasons[name]`` says why. Unexpected exceptions are not swallowed.
    """
    if names is None:
        specs = REGISTRY
    else:
        wanted = set(names)
        unknown = sorted(wanted - _BY_NAME.keys())
        if unknown:
            raise KeyError(f"unknown index names: {unknown}")
        specs = tuple(s for s in REGISTRY if s.name in wanted)
    ctx = GraphContext(graph)
    values: dict[str, Number] = {}
    reasons: dict[str, str] = {}
    for spec in specs:
        try:
            values[spec.name] = spec.func(ctx)
        except UndefinedIndex as exc:
            values[spec.name] = float("nan")
            reasons[spec.name] = str(exc)
    return values, reasons
