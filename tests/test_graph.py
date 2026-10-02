"""Tests for topoqspr.graph: standardisation, rejections and the H-suppressed graph."""

from __future__ import annotations

import pytest
from rdkit import Chem

from topoqspr.graph import (
    ORGANIC_ELEMENTS,
    REJECTION_CODES,
    MolGraph,
    Rejection,
    graph_from_mol,
    standardize,
)

ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"
ASPIRIN_KEY = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


# ---------------------------------------------------------------------------
# Accepted structures
# ---------------------------------------------------------------------------
def test_aspirin_graph_and_provenance():
    r = standardize(ASPIRIN)
    assert r.ok and r.rejection is None
    assert r.canonical_smiles == "CC(=O)Oc1ccccc1C(=O)O"
    assert r.inchikey == ASPIRIN_KEY and len(r.inchikey) == 27
    g = r.graph
    assert (g.n, g.m) == (13, 13)
    assert sum(g.degrees) == 2 * g.m
    assert sorted(g.degrees) == sorted(a.GetDegree() for a in Chem.MolFromSmiles(ASPIRIN).GetAtoms())
    assert list(g.edges) == sorted(g.edges) and all(i < j for i, j in g.edges)
    assert g.is_connected and g.max_degree == 3


@pytest.mark.parametrize(
    "salt, parent_key",
    [
        ("CC(=O)OC1=CC=CC=C1C(=O)[O-].[Na+]", ASPIRIN_KEY),  # sodium salt
        ("CC(=O)OC1=CC=CC=C1C(=O)O.O.O", ASPIRIN_KEY),  # dihydrate
        ("CN(C)C(=N)N=C(N)N.Cl", "XZWYZXLIPXDOLR-UHFFFAOYSA-N"),  # metformin HCl
    ],
)
def test_salts_and_solvents_are_stripped_to_the_parent(salt, parent_key):
    r = standardize(salt)
    assert r.ok, r.rejection
    assert r.inchikey == parent_key
    assert r.graph.n_components == 1


def test_salt_and_free_base_give_identical_graphs():
    a = standardize("CN(C)C(=N)N=C(N)N").graph
    b = standardize("CN(C)C(=N)N=C(N)N.Cl").graph
    assert a == b


def test_hydrogens_are_never_vertices():
    plain = standardize("CC(=O)O").graph
    explicit = standardize("[H]C([H])([H])C(=O)O[H]").graph
    deuterated = standardize("[2H]C([2H])([2H])C(=O)O").graph
    assert plain == explicit == deuterated
    assert plain.n == 4


def test_selenium_boron_silicon_are_organic_by_default():
    for smi in ("CC[Se]CC", "OB(O)c1ccccc1", "C[Si](C)(C)C"):
        assert standardize(smi).ok, smi


# ---------------------------------------------------------------------------
# Rejections
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "smiles, code",
    [
        ("", "unparsable"),
        ("not a smiles", "unparsable"),
        ("C1CC", "unparsable"),  # unclosed ring
        ("c1cccc1", "unparsable"),  # cannot be kekulised
        ("[Na+].[Cl-]", "disallowed_element"),  # inorganic salt
        ("[Li+].[Li+].[O-]C([O-])=O", "disallowed_element"),  # lithium carbonate
        ("C[Pt](Cl)(N)N", "disallowed_element"),  # metal complex
        ("OS(=O)(=O)O", "no_carbon"),  # inorganic acid
        ("[H][H]", "no_carbon"),
        ("c1ccccc1.c1ccncc1", "mixture"),  # two organic components
        ("CC(=O)[O-].[NH4+]", "mixture"),  # ammonia is not in the salt list
        ("CO", "too_few_heavy_atoms"),
        ("C", "too_few_heavy_atoms"),
    ],
)
def test_rejection_reasons(smiles, code):
    r = standardize(smiles)
    assert not r.ok
    assert r.graph is None
    assert r.rejection.code == code
    assert r.rejection.detail


def test_rejection_still_records_parent_provenance():
    r = standardize("c1ccccc1.c1ccncc1")
    assert r.rejection.code == "mixture"
    assert r.canonical_smiles == "c1ccccc1.c1ccncc1"
    assert r.inchikey is not None


def test_allowed_elements_and_min_heavy_atoms_are_configurable():
    assert standardize("CC[Se]CC", allowed_elements=ORGANIC_ELEMENTS - {"Se"}).rejection.code == (
        "disallowed_element"
    )
    assert standardize("C[Pt](Cl)(N)N", allowed_elements=ORGANIC_ELEMENTS | {"Pt"}).ok
    assert standardize("CO", min_heavy_atoms=2).ok
    assert standardize("CCC", min_heavy_atoms=4).rejection.code == "too_few_heavy_atoms"


def test_every_rejection_code_is_known():
    with pytest.raises(ValueError):
        Rejection("made_up", "x")
    for code in REJECTION_CODES:
        assert Rejection(code, "x").code == code


def test_disconnected_heavy_atom_graph_is_detected():
    # A molecule whose heavy atoms are linked only through hydrogen; built by
    # hand because no sane SMILES produces it.
    rw = Chem.RWMol()
    for sym in ("C", "C", "H", "C"):
        rw.AddAtom(Chem.Atom(sym))
    rw.AddBond(0, 1, Chem.BondType.SINGLE)
    rw.AddBond(1, 2, Chem.BondType.SINGLE)
    rw.AddBond(2, 3, Chem.BondType.SINGLE)
    g = graph_from_mol(rw.GetMol())
    assert g.n == 3 and g.edges == ((0, 1),)
    assert not g.is_connected and g.n_components == 2


# ---------------------------------------------------------------------------
# MolGraph
# ---------------------------------------------------------------------------
def test_molgraph_normalises_edges():
    g1 = MolGraph.from_edges(4, [(3, 2), (0, 1), (1, 2)])
    g2 = MolGraph.from_edges(4, [(1, 0), (2, 1), (2, 3)])
    assert g1 == g2
    assert g1.edges == ((0, 1), (1, 2), (2, 3))
    assert g1.degrees == (1, 2, 2, 1)
    assert g1.adjacency == ((1,), (0, 2), (1, 3), (2,))


@pytest.mark.parametrize(
    "n, edges",
    [(3, [(0, 0)]), (3, [(0, 1), (1, 0)]), (3, [(0, 3)]), (2, [(-1, 1)])],
)
def test_molgraph_rejects_invalid_edges(n, edges):
    with pytest.raises(ValueError):
        MolGraph.from_edges(n, edges)


def test_molgraph_rejects_inconsistent_degrees():
    with pytest.raises(ValueError):
        MolGraph(n=2, edges=((0, 1),), degrees=(1, 2))


def test_relabel_is_an_isomorphism():
    g = standardize(ASPIRIN).graph
    perm = list(reversed(range(g.n)))
    h = g.relabel(perm)
    assert sorted(h.degrees) == sorted(g.degrees)
    assert h.relabel(perm) == g
    with pytest.raises(ValueError):
        g.relabel([0] * g.n)


def test_renumbered_molecule_gives_isomorphic_graph():
    mol = Chem.MolFromSmiles(ASPIRIN)
    order = [5, 3, 0, 12, 1, 7, 2, 4, 6, 8, 9, 11, 10]
    g1 = standardize(ASPIRIN).graph
    g2 = standardize(Chem.MolToSmiles(Chem.RenumberAtoms(mol, order), canonical=False)).graph
    assert sorted(g1.degrees) == sorted(g2.degrees)
    assert g1.m == g2.m


def test_standardize_is_deterministic():
    a, b = standardize(ASPIRIN), standardize(ASPIRIN)
    assert a.graph == b.graph and a.canonical_smiles == b.canonical_smiles and a.inchikey == b.inchikey
