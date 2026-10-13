"""Independent reference implementation of the 104 TopoQSPR topological indices.

Written from the mathematical definitions (docs/INDEX_VERIFICATION.md, the
LaTeX strings of the registry, and the primary definitions summarised there),
NOT from src/topoqspr/indices.py.  No code of the package is imported here.

Number systems
--------------
* Rational-valued quantities are computed exactly with fractions.Fraction.
* Irrational quantities (square roots, logarithms) are computed with
  decimal.Decimal at 60 significant digits, which is used as the "exact"
  reference (the package works in binary64, i.e. ~16 digits).

An undefined value is represented by the Undefined class (with a reason).
"""
from __future__ import annotations

from collections import deque
from decimal import Decimal, getcontext
from fractions import Fraction
from functools import lru_cache
from itertools import combinations

import networkx as nx

getcontext().prec = 60
LN2 = Decimal(2).ln()


class Undefined:
    def __init__(self, reason: str):
        self.reason = reason

    def __repr__(self):
        return f"Undefined({self.reason!r})"


def D(x) -> Decimal:
    """Exact rational -> Decimal at context precision."""
    if isinstance(x, Decimal):
        return x
    x = Fraction(x)
    return Decimal(x.numerator) / Decimal(x.denominator)


def SQ(x):
    """Square root. Raises TypeError for Fractions so that a phi using it is
    recognised as irrational and re-evaluated in Decimal."""
    if isinstance(x, Decimal):
        if x < 0:
            raise ArithmeticError("sqrt of negative")
        return x.sqrt()
    raise TypeError("irrational")


# ---------------------------------------------------------------------------
# Edge functions phi(a, b), typed from the formulas (a = x_u, b = x_v)
# ---------------------------------------------------------------------------
PHI = {
    "M1": lambda a, b: a + b,
    "M2": lambda a, b: a * b,
    "F": lambda a, b: a * a + b * b,
    "HM": lambda a, b: (a + b) ** 2,
    "R": lambda a, b: 1 / SQ(a * b),
    "RR": lambda a, b: SQ(a * b),
    "SC": lambda a, b: 1 / SQ(a + b),
    "H": lambda a, b: 2 / (a + b),
    "GA": lambda a, b: 2 * SQ(a * b) / (a + b),
    "ABC": lambda a, b: SQ((a + b - 2) / (a * b)),
    "AZI": lambda a, b: (a * b / (a + b - 2)) ** 3,
    "SO": lambda a, b: SQ(a * a + b * b),
    "ISI": lambda a, b: a * b / (a + b),
    "ReZ1": lambda a, b: (a + b) / (a * b),
    "mM2": lambda a, b: 1 / (a * b),
    "SDD": lambda a, b: (a * a + b * b) / (a * b),
    "sigma": lambda a, b: (a - b) ** 2,
    "ABS": lambda a, b: SQ((a + b - 2) / (a + b)),
    "ReZ3": lambda a, b: a * b * (a + b),
    "AG": lambda a, b: (a + b) / (2 * SQ(a * b)),
    "HG": lambda a, b: 2 / (SQ(a * b) * (a + b)),
    "GBM": lambda a, b: SQ(a * b) / (a + b + a * b),
    "BMG": lambda a, b: (a + b + a * b) / SQ(a * b),
    "BMH": lambda a, b: (a + b + a * b) * (a + b) / 2,
    "TMH": lambda a, b: (a * a + b * b + a * b) * (a + b) / 2,
    "BM": lambda a, b: a + b + a * b,
    "TM": lambda a, b: a * a + b * b + a * b,
    "GH": lambda a, b: SQ(a * b) * (a + b) / 2,
    "GTM": lambda a, b: SQ(a * b) / (a * a + b * b + a * b),
    "HBM": lambda a, b: 2 / ((a + b + a * b) * (a + b)),
    "HTM": lambda a, b: 2 / ((a * a + b * b + a * b) * (a + b)),
    "TMG": lambda a, b: (a * a + b * b + a * b) / SQ(a * b),
    "HCM2": lambda a, b: (a * b) ** 2,
}

DEGREE_KEYS = ["M1", "M2", "F", "HM", "R", "RR", "SC", "H", "GA", "ABC", "AZI", "SO", "ISI",
               "ReZ1", "mM2", "SDD", "sigma", "ABS", "ReZ3", "AG", "HG", "GBM", "BMG", "BMH",
               "TMH", "BM", "TM", "GH", "GTM", "HBM", "HTM", "TMG"]
RK_KEYS = ["M1", "M2", "F", "HM", "ReZ3", "BM", "TM", "BMH", "TMH", "R", "ABC", "ISI", "RR",
           "GA", "AG", "ReZ1"]
S_KEYS = ["ABC", "GA", "RR", "SC", "ReZ3", "R", "SDD", "H", "GBM", "TMH"]


def edge_sum(key, weights, edges):
    """sum over edges of PHI[key](w_u, w_v); exact Fraction if rational, else Decimal."""
    phi = PHI[key]
    try:
        tot = Fraction(0)
        for u, v in edges:
            tot += phi(Fraction(weights[u]), Fraction(weights[v]))
        return tot
    except TypeError:
        pass
    except ZeroDivisionError:
        return Undefined(f"division by zero in {key}")
    try:
        terms = [phi(D(weights[u]), D(weights[v])) for u, v in edges]
    except ArithmeticError as exc:
        return Undefined(f"{key}: {exc!r}")
    tot = Decimal(0)
    for t in terms:
        tot += t
    return tot


def reverse_degree(d: int, k: int, delta: int):
    """M_kR(d) = Delta - d + k if k <= d, else (Delta - d + k) mod Delta (Kalaam 2024;
    mod applied to the whole sum, Godlin & Radha 2026 Eq. 1). Undefined for k > Delta."""
    if delta == 0 or k > delta:
        return None
    return delta - d + k if k <= d else (delta - d + k) % delta


# ---------------------------------------------------------------------------
# Distance machinery
# ---------------------------------------------------------------------------
def bfs(G, s):
    dist = {s: 0}
    q = deque([s])
    while q:
        x = q.popleft()
        for y in G[x]:
            if y not in dist:
                dist[y] = dist[x] + 1
                q.append(y)
    return dist


def all_dist(G):
    return {v: bfs(G, v) for v in G}


def wiener_of(G):
    dm = all_dist(G)
    nodes = list(G)
    return sum(dm[u][v] for u, v in combinations(nodes, 2))


# ---------------------------------------------------------------------------
# Hosoya Z: two independent counting methods
# ---------------------------------------------------------------------------
def hosoya_bruteforce(G) -> int:
    """Enumerate matchings edge by edge (include / exclude); small graphs only."""
    edges = list(G.edges())

    def rec(i, used):
        if i == len(edges):
            return 1
        u, v = edges[i]
        c = rec(i + 1, used)
        if u not in used and v not in used:
            c += rec(i + 1, used | {u, v})
        return c

    return rec(0, frozenset())


def hosoya_recursive(G) -> int:
    """Z(G) = Z(G - v) + sum_{w ~ v} Z(G - v - w), on induced subgraphs, memoised on the
    vertex set, multiplicative over connected components, pivot = min-degree vertex."""
    adj = {v: frozenset(G[v]) for v in G}

    @lru_cache(maxsize=None)
    def z(vs: frozenset) -> int:
        if len(vs) <= 1:
            return 1
        # connected components
        comps, seen = [], set()
        for s in vs:
            if s in seen:
                continue
            comp, stack = {s}, [s]
            seen.add(s)
            while stack:
                x = stack.pop()
                for y in adj[x]:
                    if y in vs and y not in seen:
                        seen.add(y)
                        comp.add(y)
                        stack.append(y)
            comps.append(frozenset(comp))
        if len(comps) > 1:
            out = 1
            for c in comps:
                out *= z(c)
            return out
        v = min(vs, key=lambda x: (len(adj[x] & vs), x))
        rest = vs - {v}
        total = z(rest)
        for w in adj[v] & vs:
            total += z(rest - {w})
        return total

    return z(frozenset(G.nodes()))


# ---------------------------------------------------------------------------
# All indices
# ---------------------------------------------------------------------------
def reference_indices(G: nx.Graph, names, brute_z: bool = False) -> dict:
    """Reference values for networkx graph G (vertices arbitrary hashables)."""
    G = nx.convert_node_labels_to_integers(G)
    n, m = G.number_of_nodes(), G.number_of_edges()
    edges = list(G.edges())
    deg = dict(G.degree())
    delta = max(deg.values(), default=0)
    S = {v: sum(deg[w] for w in G[v]) for v in G}
    out = {}

    for key in DEGREE_KEYS:
        out[key] = edge_sum(key, deg, edges)
    for k in (1, 2, 3):
        mk = {v: reverse_degree(deg[v], k, delta) for v in G}
        for key in RK_KEYS:
            name = f"R{k}_{key}"
            if k > delta or delta == 0:
                out[name] = Undefined(f"k={k} > Delta={delta}")
            else:
                out[name] = edge_sum(key, mk, edges)
    c = {v: delta - deg[v] + 1 for v in G}
    out["HCM2"] = edge_sum("HCM2", c, edges)
    for key in S_KEYS:
        out[f"S_{key}"] = edge_sum(key, S, edges)

    connected = n > 0 and nx.is_connected(G)
    dist_names = ["W", "J", "Sz", "Sz_e", "PI", "PI_v", "Mo", "Mo_e", "W_e", "DD", "Gut"]
    if not connected:
        for nm in dist_names:
            out[nm] = Undefined("disconnected or empty")
    else:
        dm = all_dist(G)
        pairs = list(combinations(range(n), 2))
        out["W"] = Fraction(sum(dm[u][v] for u, v in pairs))
        out["DD"] = Fraction(sum((deg[u] + deg[v]) * dm[u][v] for u, v in pairs))
        out["Gut"] = Fraction(sum(deg[u] * deg[v] * dm[u][v] for u, v in pairs))
        Dsum = {u: sum(dm[u].values()) for u in G}
        mu = m - n + 1
        if m == 0:
            out["J"] = Undefined("no edges")
        else:
            s = Decimal(0)
            for u, v in edges:
                s += 1 / D(Dsum[u] * Dsum[v]).sqrt()
            out["J"] = D(Fraction(m, mu + 1)) * s
        sz = sze = pi = piv = mo = moe = 0
        for u, v in edges:
            nu = sum(1 for w in G if dm[w][u] < dm[w][v])
            nv = sum(1 for w in G if dm[w][v] < dm[w][u])
            mu_e = mv_e = 0
            for x, y in edges:
                du = min(dm[u][x], dm[u][y])
                dv = min(dm[v][x], dm[v][y])
                if du < dv:
                    mu_e += 1
                elif dv < du:
                    mv_e += 1
            sz += nu * nv
            piv += nu + nv
            mo += abs(nu - nv)
            sze += mu_e * mv_e
            pi += mu_e + mv_e
            moe += abs(mu_e - mv_e)
        out.update(Sz=Fraction(sz), Sz_e=Fraction(sze), PI=Fraction(pi), PI_v=Fraction(piv),
                   Mo=Fraction(mo), Mo_e=Fraction(moe))
        L = nx.line_graph(G)
        out["W_e"] = Fraction(wiener_of(L)) if L.number_of_nodes() > 0 else Fraction(0)

    out["Z"] = Fraction(hosoya_bruteforce(G) if brute_z else hosoya_recursive(G))
    if m == 0:
        out["I_deg"] = Undefined("no edges")
    else:
        tot = Decimal(0)
        for v in G:
            if deg[v] == 0:
                continue  # 0 log 0 := 0 (limit)
            p = D(Fraction(deg[v], 2 * m))
            tot -= p * p.ln() / LN2
        out["I_deg"] = tot
    missing = set(names) - set(out)
    extra = set(out) - set(names)
    if missing or extra:
        raise RuntimeError(f"name mismatch: missing {missing} extra {extra}")
    return {nm: out[nm] for nm in names}
