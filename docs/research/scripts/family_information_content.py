"""TARGET-BLIND information content of each descriptor family on the Bradley isomer set.

For the parents in formula groups with >= 5 constitutions (same inclusion rule as
AMD-2026-10-13-E3), count, for every feature family, how many within-formula constitution
pairs it DISTINGUISHES (feature vectors differ beyond a relative tolerance of 1e-9, which
absorbs summation-order noise; F9). Composition distinguishes none by construction.
Theoretical bound for any unlabeled-graph invariant: pairs whose unlabeled heavy-atom graphs
are isomorphic can never be distinguished.

Families: the registry families used in Experiment B's exploratory run (degree 32,
reverse-degree 49, neighbourhood-degree-sum 10, distance 11, matching 1 [Hosoya Z],
entropy 1), all 104 TI, the 9 degree-pair edge counts + max degree (the information that the
81 degree / reverse-degree / HCM2 indices are linear in), RDKit 2D, ECFP counts, size.

    python -I docs/research/scripts/family_information_content.py docs/research/out
"""
import json
import sys
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from topoqspr import expdata
from topoqspr import features as ft
from topoqspr.experiments import index_table, prepare_structures

RDLogger.DisableLog("rdApp.*")
REPO = Path(__file__).resolve().parents[3]


def main(out: Path) -> None:
    agg = expdata.get_target("mp_bradley", offline=True)
    keys = set(agg.inchikey)
    del agg  # target-blind from here on
    recs = expdata.load_records("mp_bradley", offline=True)[["inchikey", "canonical_smiles", "status"]]
    recs = recs[recs.inchikey.isin(keys) & (recs.status == "ok")].drop_duplicates("inchikey")
    st = prepare_structures(recs.inchikey.tolist(), recs.canonical_smiles.tolist())
    mols = dict(zip(st.ids, st.mols))
    meta = pd.DataFrame({"formula": {i: CalcMolFormula(m) for i, m in mols.items()}})
    meta["constitution"] = meta.index.str[:14]
    meta = meta.reset_index().drop_duplicates("constitution").set_index("index")
    n_cons = meta.groupby("formula").size()
    meta = meta[meta.formula.map(n_cons) >= 5]
    ids = list(meta.index)
    kmols = [mols[i] for i in ids]
    idx = pd.Index(ids)
    st_k = prepare_structures(ids, [Chem.MolToSmiles(m) for m in kmols])
    ti, _ = index_table(st_k)
    ti = ti.loc[ids]
    fams = json.loads((REPO / "results/expB/exploratory_families/study.json").read_text())["feature_sets"]
    sets = {f: ti[[c.split(".", 1)[1] for c in cols]] for f, cols in fams.items() if f.startswith("TI_")}
    sets["TI_all_104"] = ti

    def etc(m):
        d = [a.GetDegree() for a in m.GetAtoms()]
        c = {}
        for b in m.GetBonds():
            k = tuple(sorted((d[b.GetBeginAtomIdx()], d[b.GetEndAtomIdx()])))
            c[f"e{k[0]}{k[1]}"] = c.get(f"e{k[0]}{k[1]}", 0) + 1
        c["maxdeg"] = max(d) if d else 0
        return c
    sets["edge_types_9_plus_maxdeg"] = pd.DataFrame([etc(m) for m in kmols], index=idx).fillna(0)
    sets["size"] = ft.size_features(kmols, idx)
    sets["composition"] = ft.composition_features(kmols, idx)
    sets["rdkit2d"] = ft.rdkit2d_features(kmols, None, idx)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    sets["ECFP_counts"] = pd.DataFrame(np.vstack([gen.GetCountFingerprintAsNumPy(m) for m in kmols]),
                                       index=idx).astype(float)

    pairs = [(a, b) for _, g in meta.groupby("formula") for a, b in combinations(g.index, 2)]

    def graph(m):
        g = nx.Graph(); g.add_nodes_from(range(m.GetNumAtoms()))
        g.add_edges_from((x.GetBeginAtomIdx(), x.GetEndAtomIdx()) for x in m.GetBonds())
        return g
    G = {i: graph(m) for i, m in zip(ids, kmols)}
    iso = sum(nx.is_isomorphic(G[a], G[b]) for a, b in pairs)
    rows = [dict(family="unlabeled-graph isomorphism (bound for any unlabeled invariant)",
                 n_features=None, pairs=len(pairs), distinguished=len(pairs) - iso)]
    for name, X in sets.items():
        A = X.to_numpy(float)
        pos = {i: k for k, i in enumerate(ids)}
        dist = 0
        for a, b in pairs:
            u, v = A[pos[a]], A[pos[b]]
            both_nan = np.isnan(u) & np.isnan(v)
            diff = np.abs(u - v) > 1e-9 * np.maximum(1.0, np.maximum(np.abs(u), np.abs(v)))
            diff = diff & ~both_nan | (np.isnan(u) ^ np.isnan(v))
            dist += bool(diff.any())
        rows.append(dict(family=name, n_features=X.shape[1], pairs=len(pairs), distinguished=dist))
    df = pd.DataFrame(rows)
    df["share_distinguished"] = df.distinguished / df.pairs
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "family_information_content.csv", index=False)
    print(f"{len(ids)} constitutions, {meta.formula.nunique()} formula groups, {len(pairs)} within-formula pairs")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
