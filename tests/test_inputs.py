"""Structure/graph inputs: SMILES, MOL/SDF, PDB, adjacency matrix and list."""

import pytest
from rdkit import Chem

from topoqspr.graph import standardize
from topoqspr.indices import compute_indices
from topoqspr.inputs import FORMATS, parse_input

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


def _m1(rec):
    return compute_indices(rec.graph)[0]["M1"]


def test_every_molecular_route_gives_the_smiles_graph():
    ref = compute_indices(standardize(ASPIRIN).graph)[0]
    mol = Chem.MolFromSmiles(ASPIRIN)
    for fmt, text in [("smiles", ASPIRIN), ("mol", Chem.MolToMolBlock(mol)),
                      ("sdf", Chem.MolToMolBlock(mol) + "$$$$\n"), ("pdb", Chem.MolToPDBBlock(mol))]:
        (rec,) = parse_input(text, fmt)
        vals = compute_indices(rec.graph)[0]
        assert vals["M1"] == ref["M1"] and vals["W"] == ref["W"], fmt


def test_adjacency_inputs_closed_forms():
    (c4,) = parse_input("A = [[0,1,0,1],[1,0,1,0],[0,1,0,1],[1,0,1,0]]", "adjacency_matrix")
    assert c4.kind == "graph" and _m1(c4) == 16
    (c6,) = parse_input('{"1":[2,6],"2":[1,3],"3":[2,4],"4":[3,5],"5":[4,6],"6":[5,1]}', "adjacency_list")
    assert _m1(c6) == 24 and compute_indices(c6.graph)[0]["W"] == 27  # W(C6) = 27
    (p3,) = parse_input("0 1 0\n1 0 1\n0 1 0", "adjacency_matrix")
    assert compute_indices(p3.graph)[0]["W"] == 4  # (n^3 - n)/6 = 4


@pytest.mark.parametrize("text,fmt,reason", [
    ("0 1\n1 1", "adjacency_matrix", "self-loops"),
    ("0 1 0\n0 0 1\n0 1 0", "adjacency_matrix", "not symmetric"),
    ("0 1\n1 0 1", "adjacency_matrix", "not square"),
    ("0 2\n2 0", "adjacency_matrix", "only 0 and 1"),
    ("1: 2\n3: 4", "adjacency_list", "disconnected"),
    ("C1CC", "smiles", "unparsable"),
])
def test_rejections_carry_reasons(text, fmt, reason):
    (rec,) = parse_input(text, fmt)
    assert rec.graph is None and reason in rec.reason


def test_sdf_multirecord_and_names():
    blk = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))
    named = "ethanol\n" + blk.split("\n", 1)[1]
    recs = parse_input(named + "$$$$\n" + blk + "$$$$\n", "sdf")
    assert [r.label for r in recs] == ["ethanol", "input_2"] and all(r.graph for r in recs)


def test_unknown_format():
    with pytest.raises(ValueError):
        parse_input("x", "xyz")
    assert "smiles" in FORMATS


def test_adjacency_list_accepts_edge_lines_and_names_the_formats():
    """One edge per line (the most common hand-typed form) gives the same graph as the other two."""
    (a,) = parse_input("1 2\n2 3\n3 4\n4 5\n5 6\n6 1", "adjacency_list")
    (b,) = parse_input("1: 2 6\n2: 3\n3: 4\n4: 5\n5: 6", "adjacency_list")
    (c,) = parse_input('{"1":[2,6],"2":[1,3],"3":[2,4],"4":[3,5],"5":[4,6],"6":[5,1]}', "adjacency_list")
    assert a.graph is not None and b.graph is not None
    assert sorted(a.graph.edges) == sorted(b.graph.edges) == sorted(c.graph.edges) and a.graph.n == 6
    (bad,) = parse_input("1 2 3\n", "adjacency_list")
    assert bad.graph is None and "accepted:" in bad.reason and "neither" in bad.reason
