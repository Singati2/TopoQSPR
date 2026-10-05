"""Bond-wise (edge) entropies of the edge-additive topological indices.

For an index T(G) = sum_{uv in E} w(uv) with non-negative edge contributions w, the entropy
of the edge-weight distribution is

    ENT_T(G) = - sum_{uv in E, w(uv) > 0} (w(uv) / T) ln(w(uv) / T)

(Shannon entropy with natural logarithm; edges with w = 0 contribute 0). This is the
definition implemented by ChemGraphX (Jacob et al., J. Cheminform. 18:39, 2026;
``server/basic_values.py: entropy_sum``) and used for degree-based indices in the graph-entropy
literature (Dehmer & Mowshowitz 2011). ENT_T is undefined (NaN with a reason) when T = 0 or
when an edge contribution is negative or not finite.

These entropies are an **opt-in** extra for the web interface and the command line. They are
deliberately *not* part of :data:`topoqspr.indices.REGISTRY`, so the frozen study
(``plan-freeze-v1``), whose ``TI`` feature set is the full registry, is unchanged.
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

from .graph import MolGraph
from .indices import (_DEGREE, _NBR, _RK, PHI, RK_VALUES, GraphContext, UndefinedIndex)

__all__ = ["ENTROPY_NAMES", "edge_entropy", "index_entropies"]


def edge_entropy(contributions: Sequence[float]) -> float:
    """-sum p ln p over the positive contributions, p = w / sum(w)."""
    ws = [float(w) for w in contributions]
    if any(not math.isfinite(w) or w < 0 for w in ws):
        raise UndefinedIndex("an edge contribution is negative or not finite")
    total = math.fsum(ws)
    if total == 0:
        raise UndefinedIndex("index value is 0, so the edge distribution is undefined")
    return -math.fsum((w / total) * math.log(w / total) for w in ws if w > 0)


def _specs() -> list[tuple[str, str, Callable[[GraphContext], Sequence[int]]]]:
    out = []
    for name, key, _ref in _DEGREE:
        out.append((f"ENT_{name}", key, lambda c: c.deg))
    for k in RK_VALUES:
        for key, _ref in _RK:
            out.append((f"ENT_R{k}_{key}", key, (lambda kk: (lambda c: c.rk(kk)))(k)))
    out.append(("ENT_HCM2", "HCM2", lambda c: c.rk(1)))  # Kulli c_v = Delta - d_v + 1 = M_1R
    for name, key, _ref, _alias in _NBR:
        out.append((f"ENT_{name}", key, lambda c: c.nbr_sum))
    return out


_SPECS = _specs()
ENTROPY_NAMES: tuple[str, ...] = tuple(s[0] for s in _SPECS)


def index_entropies(graph: MolGraph) -> tuple[dict[str, float], dict[str, str]]:
    """Edge entropies for every edge-additive index. Returns (values, reasons for NaN)."""
    ctx = GraphContext(graph)
    values: dict[str, float] = {}
    reasons: dict[str, str] = {}
    for name, key, weights_of in _SPECS:
        phi = PHI[key][0]
        try:
            w = weights_of(ctx)
            contrib = []
            for u, v in graph.edges:
                try:
                    contrib.append(phi(w[u], w[v]))
                except ZeroDivisionError:
                    raise UndefinedIndex(f"division by zero on edge ({u},{v})") from None
                except ValueError as exc:
                    raise UndefinedIndex(f"{exc} on edge ({u},{v})") from None
            values[name] = edge_entropy(contrib)
        except UndefinedIndex as exc:
            values[name] = float("nan")
            reasons[name] = str(exc)
    return values, reasons
