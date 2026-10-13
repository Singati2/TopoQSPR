"""Gate C: TARGET-BLIND census of molecular-formula isomer groups in full measured datasets.

Uses the package's own record loading, standardisation and aggregation
(expdata.get_target, offline, default rules incl. dropping replicate-disagreeing parents),
then immediately discards the `value` and `sd` columns: only structure identifiers and the
fact that a parent has a usable measurement are kept. No target value is printed, summarised
or written.

Per endpoint: standardised parents; formula groups with >= 2/3/5/10 distinct constitutions
(InChIKey first block); molecules in those groups; stereo-only duplicates (same first block,
different full key); unlabeled heavy-atom graph classes inside formula groups (isomorphic
=> every unlabeled invariant collides); single-source share.

    python -I docs/research/scripts/isomer_census.py docs/research/out
"""
import json
import sys
from pathlib import Path

import networkx as nx
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from topoqspr import expdata

RDLogger.DisableLog("rdApp.*")
ENDPOINTS = ["mp_bradley", "logp_opera", "logs_aqsoldb"]


def graph(smiles: str) -> nx.Graph:
    m = Chem.MolFromSmiles(smiles)
    g = nx.Graph()
    g.add_nodes_from(range(m.GetNumAtoms()))
    g.add_edges_from((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds())
    return g


def census(name: str) -> tuple[dict, pd.DataFrame]:
    agg = expdata.get_target(name, offline=True)
    keys = set(agg["inchikey"])
    n_measured = len(keys)
    recs = expdata.load_records(name, offline=True)[["inchikey", "canonical_smiles", "status"]]
    agg = agg.drop(columns=[c for c in ("value", "sd") if c in agg.columns])  # target-blind
    del agg
    recs = recs[recs["inchikey"].isin(keys)]
    smi_col = next(c for c in ("canonical_smiles", "smiles_std", "std_smiles", "smiles") if c in recs.columns)
    rep = recs.drop_duplicates("inchikey").set_index("inchikey")[smi_col]
    df = pd.DataFrame({"smiles": rep})
    df["formula"] = [CalcMolFormula(Chem.MolFromSmiles(s)) for s in df.smiles]
    df["constitution"] = df.index.str[:14]
    stereo_dup = int(df.groupby("constitution").size().gt(1).sum())
    cons = df.reset_index().drop_duplicates("constitution")
    sizes = cons.groupby("formula").size()
    out = dict(endpoint=name, n_measured_parents=n_measured, n_constitutions=len(cons),
               constitutions_with_stereo_variants=stereo_dup,
               formula_groups_ge2=int((sizes >= 2).sum()), ge3=int((sizes >= 3).sum()),
               ge5=int((sizes >= 5).sum()), ge10=int((sizes >= 10).sum()),
               molecules_in_ge5=int(sizes[sizes >= 5].sum()),
               molecules_in_ge2=int(sizes[sizes >= 2].sum()),
               largest_group=int(sizes.max()), largest_formula=str(sizes.idxmax()))
    # unlabeled-graph collision classes inside formula groups >= 2
    n_pairs = n_iso_pairs = 0
    classes_rows = []
    for f, grp in cons[cons.formula.isin(sizes[sizes >= 2].index)].groupby("formula"):
        reps, cls = [], []
        for s in grp.smiles:
            g = graph(s)
            for k, r in enumerate(reps):
                if nx.is_isomorphic(g, r):
                    cls.append(k); break
            else:
                reps.append(g); cls.append(len(reps) - 1)
        counts = pd.Series(cls).value_counts()
        n = len(cls)
        n_pairs += n * (n - 1) // 2
        n_iso_pairs += int(sum(c * (c - 1) // 2 for c in counts))
        classes_rows.append(dict(formula=f, n_constitutions=n, n_unlabeled_graphs=len(reps),
                                 largest_collision_class=int(counts.max())))
    out["within_formula_pairs"] = n_pairs
    out["within_formula_pairs_unlabeled_isomorphic"] = n_iso_pairs
    cls = pd.DataFrame(classes_rows)
    out["formula_groups_ge5_with_ge2_distinct_unlabeled_graphs"] = int(
        ((cls.n_constitutions >= 5) & (cls.n_unlabeled_graphs >= 2)).sum())
    return out, cls.assign(endpoint=name)


def main(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    rows, classes = [], []
    for e in ENDPOINTS:
        o, c = census(e)
        print(json.dumps(o), flush=True)
        rows.append(o); classes.append(c)
    pd.DataFrame(rows).to_csv(outdir / "isomer_census.csv", index=False)
    pd.concat(classes).to_csv(outdir / "isomer_census_formula_groups.csv", index=False)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
