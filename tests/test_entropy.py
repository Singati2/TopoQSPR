"""Edge entropies (opt-in, not in the frozen registry), checked against ChemGraphX's code."""

import json
import math
from pathlib import Path

import pytest

from topoqspr.entropy import ENTROPY_NAMES, edge_entropy, index_entropies
from topoqspr.graph import MolGraph
from topoqspr.indices import INDEX_NAMES

FIX = json.loads((Path(__file__).parent / "fixtures" / "entropy" / "chemgraphx_entropies.json").read_text())


def test_entropies_are_not_in_the_frozen_registry():
    assert not set(ENTROPY_NAMES) & set(INDEX_NAMES)
    assert len(ENTROPY_NAMES) == 89


def test_regular_graph_entropy_is_ln_m():
    g = MolGraph.from_edges(6, [(i, (i + 1) % 6) for i in range(6)])
    v, why = index_entropies(g)
    for k, x in v.items():
        if math.isnan(x):
            assert k in why  # e.g. sigma = 0 on a regular graph
        else:
            assert x == pytest.approx(math.log(6), rel=1e-14)


def test_undefined_cases_are_nan_with_reason():
    with pytest.raises(Exception):
        edge_entropy([0.0, 0.0])
    with pytest.raises(Exception):
        edge_entropy([1.0, -1.0])
    assert edge_entropy([1.0, 0.0, 1.0]) == pytest.approx(math.log(2))


@pytest.mark.parametrize("entry", FIX["graphs"], ids=lambda e: e["name"])
def test_matches_chemgraphx(entry):
    g = MolGraph.from_edges(entry["n"], [tuple(e) for e in entry["edges"]])
    ours, _ = index_entropies(g)
    for name, ref in entry["chemgraphx"].items():
        assert ours[name] == pytest.approx(ref, rel=1e-10), name
