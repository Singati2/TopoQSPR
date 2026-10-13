"""Value-by-value comparison of TopoQSPR's 104 indices with an independent implementation.

Usage (from anywhere):
    ~/.venvs/topoqspr/bin/python -I docs/research/index_verification/verify_indices.py \
        ~/topoqspr_runs/main_confirmatory/molecules.csv docs/research/index_verification/out

The reference implementation is refimpl.py (same directory), written from the
definitions only. The package is used through its public API:
  * topoqspr.graph.MolGraph.from_edges + topoqspr.indices.compute_indices (abstract graphs)
  * topoqspr.experiments.prepare_structures + index_table (SMILES route)
Only REGISTRY names are read from topoqspr.indices (no computation code is reused).
"""
from __future__ import annotations

import csv
import json
import math
import os
import sys
import time
from collections import defaultdict
from decimal import Decimal
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import networkx as nx  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402

import refimpl as ref  # noqa: E402
from refimpl import D, Undefined  # noqa: E402
from topoqspr.experiments import index_table, prepare_structures, size_table  # noqa: E402
from topoqspr.graph import MolGraph  # noqa: E402
from topoqspr.indices import REGISTRY, compute_indices  # noqa: E402

RDLogger.DisableLog("rdApp.*")
NAMES = [s.name for s in REGISTRY]
REL_TOL = 1e-12
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


# ---------------------------------------------------------------------------
# Graph helpers
# ---------------------------------------------------------------------------
def nx_to_molgraph(G):
    G = nx.convert_node_labels_to_integers(G)
    return MolGraph.from_edges(G.number_of_nodes(), list(G.edges()))


def molgraph_to_nx(g):
    G = nx.Graph()
    G.add_nodes_from(range(g.n))
    G.add_edges_from(g.edges)
    return G


def nx_to_smiles(G):
    """Realise a graph (max degree 4) as a saturated hydrocarbon."""
    G = nx.convert_node_labels_to_integers(G)
    rw = Chem.RWMol()
    for _ in G.nodes():
        rw.AddAtom(Chem.Atom(6))
    for u, v in G.edges():
        rw.AddBond(int(u), int(v), Chem.BondType.SINGLE)
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(mol)


def smiles_to_nx(smi):
    """My own hydrogen-suppressed graph of a SMILES (heavy atoms, heavy-heavy bonds)."""
    mol = Chem.MolFromSmiles(smi)
    G = nx.Graph()
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    G.add_nodes_from(heavy)
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in G and j in G:
            G.add_edge(i, j)
    return G


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------
def is_nan(x):
    return isinstance(x, float) and math.isnan(x)


def compare(refv, pkgv, pkg_reason):
    """Return (status, absdiff, reldiff)."""
    if isinstance(refv, Undefined):
        if is_nan(pkgv):
            return ("both_undefined" if pkg_reason else "both_undefined_NO_REASON"), 0.0, 0.0
        return "ref_undefined_pkg_value", math.inf, math.inf
    if is_nan(pkgv):
        return "pkg_nan_ref_value", math.inf, math.inf
    if isinstance(refv, Fraction):
        if isinstance(pkgv, int) and not isinstance(pkgv, bool):
            diff = abs(Fraction(pkgv) - refv)
        else:
            diff = abs(Fraction(float(pkgv)) - refv)
        a = float(diff)
        r = float(diff / abs(refv)) if refv != 0 else (0.0 if diff == 0 else math.inf)
    else:
        diff = abs(Decimal(float(pkgv)) - refv) if not isinstance(pkgv, int) else abs(D(pkgv) - refv)
        a = float(diff)
        r = float(diff / abs(refv)) if refv != 0 else (0.0 if diff == 0 else math.inf)
    ok = (r <= REL_TOL) or (a <= 1e-12 and refv == 0)
    return ("agree" if ok else "DISAGREE"), a, r


def fmt(v):
    if isinstance(v, Undefined):
        return f"undefined({v.reason})"
    if isinstance(v, Fraction):
        return str(v)
    if isinstance(v, Decimal):
        return f"{v:.25g}"
    return repr(v)


# ---------------------------------------------------------------------------
# Test sets
# ---------------------------------------------------------------------------
def atlas_graphs():
    out, counts = [], defaultdict(lambda: [0, 0])
    for i, G in enumerate(nx.graph_atlas_g()):
        n = G.number_of_nodes()
        if n < 2 or n > 7 or not nx.is_connected(G):
            continue
        counts[n][0] += 1
        if max(d for _, d in G.degree()) > 4:
            continue
        counts[n][1] += 1
        out.append((f"atlas{i}", G))
    return out, dict(counts)


def closed_form_graphs(nmax=30):
    gs = []
    for n in range(2, nmax + 1):
        gs.append((f"P{n}", "path", n, nx.path_graph(n)))
    for n in range(3, nmax + 1):
        gs.append((f"C{n}", "cycle", n, nx.cycle_graph(n)))
    for s in range(2, nmax + 1):
        gs.append((f"S{s}", "star", s, nx.star_graph(s)))
    return gs


def fib(k):
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def lucas(k):
    a, b = 2, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def log2D(x):
    return D(x).ln() / ref.LN2


def hand_closed_forms(kind, p):
    """Closed forms derived by hand (not via refimpl edge functions)."""
    F_ = Fraction
    c = {}
    if kind == "path" and p >= 3:
        n = p
        W = F_(n * (n * n - 1), 6)
        c.update(M1=F_(4 * n - 6), M2=F_(4 * n - 8), F=F_(8 * n - 14), sigma=F_(2), W=W, Sz=W,
                 PI_v=F_(n * (n - 1)), Mo=F_((n - 1) ** 2 // 2), PI=F_((n - 1) * (n - 2)),
                 Mo_e=F_((n - 1) ** 2 // 2),
                 W_e=F_((n - 1) * ((n - 1) ** 2 - 1), 6), DD=4 * W - n * (n - 1),
                 Gut=4 * W - (2 * n - 1) * (n - 1), Z=F_(fib(n + 1)),
                 R=D(2).sqrt() + D(F_(n - 3, 2)), H=F_(4, 3) + F_(n - 3, 2),
                 ReZ1=F_(n), HCM2=F_(2 * 2 ** 2 + (n - 3) * 1))
        m = n - 1
        c["I_deg"] = -(2 * D(F_(1, 2 * m)) * log2D(F_(1, 2 * m))
                       + (n - 2) * D(F_(2, 2 * m)) * log2D(F_(2, 2 * m)))
    if kind == "cycle":
        n = p
        W = F_(n ** 3, 8) if n % 2 == 0 else F_(n * (n * n - 1), 8)
        Dv = 2 * W / n
        c.update(M1=F_(4 * n), M2=F_(4 * n), F=F_(8 * n), sigma=F_(0), R=F_(n, 2), W=W,
                 Sz=F_(n ** 3, 4) if n % 2 == 0 else F_(n * (n - 1) ** 2, 4),
                 PI_v=F_(n * n) if n % 2 == 0 else F_(n * (n - 1)),
                 PI=F_(n * (n - 2)) if n % 2 == 0 else F_(n * (n - 1)),
                 Mo=F_(0), Mo_e=F_(0), Z=F_(lucas(n)), W_e=W, ReZ1=F_(n), HCM2=F_(n),
                 DD=4 * W, Gut=4 * W,
                 J=D(F_(n * n, 2)) / D(Dv), I_deg=log2D(n))
    if kind == "star":
        s = p
        c.update(M1=F_(s * (s + 1)), M2=F_(s * s), F=F_(s + s ** 3), sigma=F_(s * (s - 1) ** 2),
                 W=F_(s * s), Sz=F_(s * s), PI_v=F_((s + 1) * s), Mo=F_(s * (s - 1)),
                 PI=F_(s * (s - 1)), Mo_e=F_(s * (s - 1)), W_e=F_(s * (s - 1), 2), Z=F_(s + 1),
                 R=D(s).sqrt(), ReZ1=F_(s + 1),
                 J=D(s) * D(s).sqrt() / D(2 * s - 1).sqrt(),
                 I_deg=D(F_(1, 2)) + D(F_(1, 2)) * log2D(2 * s),
                 HCM2=F_(s * s ** 2))  # c(centre)=1, c(leaf)=s
        # DD = sum_pairs (du+dv) D: centre-leaf s*(s+1)*1, leaf-leaf C(s,2)*2*2
        c["DD"] = F_(s * (s + 1) + 4 * (s * (s - 1) // 2))
        c["Gut"] = F_(s * s + 2 * (s * (s - 1) // 2))
    return c


def undefined_probes():
    probes = []
    probes.append(("K2", MolGraph.from_edges(2, [(0, 1)])))
    probes.append(("P3", MolGraph.from_edges(3, [(0, 1), (1, 2)])))
    probes.append(("2K2_disconnected", MolGraph.from_edges(4, [(0, 1), (2, 3)])))
    probes.append(("P3+K1_isolated", MolGraph.from_edges(4, [(0, 1), (1, 2)])))
    probes.append(("K1", MolGraph.from_edges(1, [])))
    probes.append(("empty", MolGraph.from_edges(0, [])))
    probes.append(("3K1_edgeless", MolGraph.from_edges(3, [])))
    probes.append(("K5_delta4", MolGraph.from_edges(5, [(i, j) for i in range(5) for j in range(i + 1, 5)])))
    probes.append(("K6_delta5", MolGraph.from_edges(6, [(i, j) for i in range(6) for j in range(i + 1, 6)])))
    return probes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(mol_csv, outdir):
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    rows = []          # one row per (graph, route, index)
    refvals = {}       # graph id -> {name: ref}
    meta = {}          # graph id -> dict(n, m, degseq, delta, group)
    pkgvals = {}       # (graph id, route) -> {name: pkg}

    def record(gid, group, route, G, vals, reasons):
        if gid not in refvals:
            refvals[gid] = ref.reference_indices(G, NAMES, brute_z=(G.number_of_nodes() <= 10))
            degs = sorted((d for _, d in G.degree()), reverse=True)
            meta[gid] = dict(group=group, n=G.number_of_nodes(), m=G.number_of_edges(),
                             degseq=degs, delta=max(degs, default=0),
                             connected=G.number_of_nodes() > 0 and nx.is_connected(G))
        pkgvals[(gid, route)] = vals
        for nm in NAMES:
            st, a, r = compare(refvals[gid][nm], vals[nm], reasons.get(nm))
            rows.append(dict(index=nm, graph=gid, group=group, route=route, status=st, absdiff=a,
                             reldiff=r, ref=fmt(refvals[gid][nm]), pkg=repr(vals[nm]),
                             pkg_type=type(vals[nm]).__name__, reason=reasons.get(nm, "")))

    # ---- (a) exhaustive atlas -------------------------------------------------
    atlas, counts = atlas_graphs()
    log(f"atlas: connected graphs per n (all, max-degree<=4): {counts}; total used {len(atlas)}")
    smiles_ids, smiles_list, smiles_graph = [], [], {}
    sanit_fail = []
    for gid, G in atlas:
        g = nx_to_molgraph(G)
        vals, reasons = compute_indices(g)
        record(gid, "atlas", "direct", G, vals, reasons)
        smi = nx_to_smiles(G)
        if smi is None:
            sanit_fail.append(gid)
            continue
        smiles_ids.append(gid)
        smiles_list.append(smi)
        smiles_graph[gid] = G
    log(f"atlas: RDKit sanitisation failures: {len(sanit_fail)} {sanit_fail[:10]}")

    # ---- (b) closed-form families -------------------------------------------
    closed = closed_form_graphs(30)
    for gid, kind, p, G in closed:
        vals, reasons = compute_indices(nx_to_molgraph(G))
        record(gid, kind, "direct", G, vals, reasons)
        if kind != "star" or p <= 4:
            smi = nx_to_smiles(G)
            if smi is not None:
                smiles_ids.append(gid)
                smiles_list.append(smi)
                smiles_graph[gid] = G

    # ---- SMILES route for atlas + closed forms ------------------------------
    st = prepare_structures(smiles_ids, smiles_list, min_heavy_atoms=2)
    tab = st.table
    rej = tab[tab.status != "ok"]
    log(f"SMILES route: {len(smiles_ids)} submitted, {len(st.ids)} accepted, {len(rej)} rejected")
    for _, r in rej.iterrows():
        log(f"  rejected {r['id']}: {r['reason']} {r['detail']}")
    T, reasons_list = index_table(st)
    reasons_by_id = defaultdict(dict)
    for rr in reasons_list:
        reasons_by_id[rr["row"]][rr["feature"]] = rr["reason"]
    iso_fail = []
    sz = size_table(st)
    ring_fail = []
    for gid, g in zip(st.ids, st.graphs):
        G = smiles_graph[gid]
        if not nx.is_isomorphic(G, molgraph_to_nx(g)):
            iso_fail.append(gid)
        vals = {nm: T.loc[gid, nm] for nm in NAMES}
        vals = {k: (v.item() if hasattr(v, "item") else v) for k, v in vals.items()}
        group = "atlas" if gid.startswith("atlas") else meta[gid]["group"]
        record(gid, group, "smiles", G, vals, reasons_by_id.get(gid, {}))
        if sz.loc[gid, "n_rings"] != G.number_of_edges() - G.number_of_nodes() + 1:
            ring_fail.append(gid)
    log(f"SMILES route: graph not isomorphic to intended graph: {len(iso_fail)} {iso_fail[:10]}")
    log(f"SMILES route: n_rings != m-n+1: {len(ring_fail)} {ring_fail[:10]}")
    log(f"index_table dtypes: {sorted(set(str(T[c].dtype) for c in NAMES))}")

    # ---- (c) drugs -----------------------------------------------------------
    import pandas as pd
    mols = pd.read_csv(mol_csv)
    ok = mols[mols.status == "ok"].copy()
    must = ["CHEMBL506569", "CHEMBL1751", "CHEMBL1360"]
    ok = ok.sort_values(["n_heavy", "id"]).reset_index(drop=True)
    # full cohort through the package (used for the Z > 2^53 selection and panel-wide checks)
    full = prepare_structures(list(ok["id"]), list(ok["smiles"]))
    FT, _ = index_table(full)
    fsz = size_table(full)
    big_z = list(FT["Z"].sort_values(ascending=False).index[FT["Z"].sort_values(ascending=False) > 2 ** 53])
    log(f"cohort: {len(full.ids)} molecules; package Z > 2^53 for {len(big_z)}: {big_z}")
    picks = list(must) + [b for b in big_z if b not in must]
    step = max(1, len(ok) // 24)
    for i in range(0, len(ok), step):
        if len(picks) >= 36:
            break
        rid = ok.loc[i, "id"]
        if rid not in picks:
            picks.append(rid)
    sel = ok.set_index("id").loc[picks]
    dst = prepare_structures(list(sel.index), list(sel["smiles"]))
    log(f"drugs: {len(picks)} selected, {len(dst.ids)} accepted")
    DT, dreasons = index_table(dst)
    dreasons_by_id = defaultdict(dict)
    for rr in dreasons:
        dreasons_by_id[rr["row"]][rr["feature"]] = rr["reason"]
    dsz = size_table(dst)
    drug_info = []
    for gid, g, mol in zip(dst.ids, dst.graphs, dst.mols):
        G = smiles_to_nx(sel.loc[gid, "smiles"])
        iso = nx.is_isomorphic(G, molgraph_to_nx(g))
        sssr = mol.GetRingInfo().NumRings()
        cyc = G.number_of_edges() - G.number_of_nodes() + 1
        t1 = time.time()
        vals = {nm: (DT.loc[gid, nm].item() if hasattr(DT.loc[gid, nm], "item") else DT.loc[gid, nm])
                for nm in NAMES}
        # Hosoya Z directly from the package as a Python int (no DataFrame round trip)
        pv, pr = compute_indices(g, ["Z"])
        record(gid, "drug", "smiles", G, vals, dreasons_by_id.get(gid, {}))
        tz = time.time() - t1
        zref = refvals[gid]["Z"]
        drug_info.append(dict(id=gid, n=G.number_of_nodes(), m=G.number_of_edges(), iso=iso,
                              n_rings_pkg=int(dsz.loc[gid, "n_rings"]), cyclomatic=cyc, rdkit_sssr=sssr,
                              Z_ref=str(zref), Z_pkg_direct=str(pv["Z"]),
                              Z_pkg_table=repr(vals["Z"]), Z_type_table=type(vals["Z"]).__name__,
                              Z_exact_direct=(isinstance(pv["Z"], int) and Fraction(pv["Z"]) == zref),
                              Z_exact_table=(not is_nan(vals["Z"]) and Fraction(vals["Z"]) == zref),
                              delta=meta[gid]["delta"], ref_time_s=round(tz, 3)))
    log("drugs: id n m iso rings_pkg cyclomatic sssr Z_ref Z_pkg(direct) Z_table_type exact_direct exact_table")
    for d in drug_info:
        log("  ", d["id"], d["n"], d["m"], d["iso"], d["n_rings_pkg"], d["cyclomatic"], d["rdkit_sssr"],
            d["Z_ref"], d["Z_pkg_direct"], d["Z_type_table"], d["Z_exact_direct"], d["Z_exact_table"])

    # ---- undefined probes (direct MolGraph) ---------------------------------
    probe_rows = []
    for gid, g in undefined_probes():
        vals, reasons = compute_indices(g)
        G = molgraph_to_nx(g)
        record(gid, "probe", "direct", G, vals, reasons)
        for nm in NAMES:
            if is_nan(vals[nm]) or isinstance(refvals[gid][nm], Undefined):
                probe_rows.append(dict(graph=gid, index=nm, pkg=repr(vals[nm]), reason=reasons.get(nm, ""),
                                       ref=fmt(refvals[gid][nm])))

    # ---- closed forms: hand formula vs reference and vs package -------------
    cf_rows = []
    for gid, kind, p, G in closed:
        hand = hand_closed_forms(kind, p)
        for nm, hv in hand.items():
            rv = refvals[gid][nm]
            pv = pkgvals[(gid, "direct")][nm]
            s2, a2, r2 = compare(hv, pv, "")
            # reference vs hand at high precision
            hd = D(hv) if isinstance(hv, Fraction) else hv
            rd = D(rv) if isinstance(rv, Fraction) else rv
            hi = abs(hd - rd) <= Decimal("1e-40") * max(Decimal(1), abs(hd))
            cf_rows.append(dict(graph=gid, index=nm, hand=fmt(hv), ref_matches_hand_1e40=hi,
                                pkg_vs_hand=s2, pkg_rel=r2))
    log(f"closed forms: {len(cf_rows)} (graph, index) checks; ref==hand (1e-40): "
        f"{sum(r['ref_matches_hand_1e40'] for r in cf_rows)}; pkg agrees with hand: "
        f"{sum(r['pkg_vs_hand'] == 'agree' for r in cf_rows)}")
    for r in cf_rows:
        if not r["ref_matches_hand_1e40"] or r["pkg_vs_hand"] != "agree":
            log("   closed-form mismatch:", r)

    # ---- Hosoya: brute force vs recursive cross-check on atlas --------------
    zb = 0
    for gid, G in atlas:
        if ref.hosoya_bruteforce(G) != ref.hosoya_recursive(G):
            zb += 1
    log(f"Hosoya brute-force vs recursion disagreements on atlas: {zb} of {len(atlas)}")

    # ---- summaries -----------------------------------------------------------
    per = {}
    for nm in NAMES:
        rs = [r for r in rows if r["index"] == nm]
        compared = [r for r in rs if r["status"] in ("agree", "DISAGREE")]
        bad = [r for r in rs if r["status"] not in ("agree", "both_undefined")]
        mx_a = max((r["absdiff"] for r in compared), default=0.0)
        mx_r = max((r["reldiff"] for r in compared), default=0.0)
        worst = max(compared, key=lambda r: r["reldiff"], default=None)
        types = sorted({r["pkg_type"] for r in rs if r["route"] == "direct"})
        per[nm] = dict(index=nm, n_rows=len(rs),
                       n_graphs=len({r["graph"] for r in rs}),
                       n_value_comparisons=len(compared),
                       n_both_undefined=sum(r["status"] == "both_undefined" for r in rs),
                       n_problem=len(bad), max_abs=mx_a, max_rel=mx_r,
                       worst_graph=worst["graph"] if worst else "", direct_types="/".join(types),
                       problem_statuses="/".join(sorted({r["status"] for r in bad})),
                       first_problem=(f"{bad[0]['graph']}[{bad[0]['route']}]: ref={bad[0]['ref']} "
                                      f"pkg={bad[0]['pkg']}") if bad else "")
    with open(os.path.join(outdir, "summary_per_index.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(next(iter(per.values())).keys()))
        w.writeheader()
        w.writerows(per.values())
    with open(os.path.join(outdir, "problems.csv"), "w", newline="") as fh:
        bad = [r for r in rows if r["status"] not in ("agree", "both_undefined")]
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(bad)
    with open(os.path.join(outdir, "closed_forms.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cf_rows[0].keys()))
        w.writeheader()
        w.writerows(cf_rows)
    with open(os.path.join(outdir, "undefined_probes.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(probe_rows[0].keys()))
        w.writeheader()
        w.writerows(probe_rows)
    with open(os.path.join(outdir, "drugs.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(drug_info[0].keys()))
        w.writeheader()
        w.writerows(drug_info)

    n_ok = sum(1 for v in per.values() if v["n_problem"] == 0)
    log(f"rows: {len(rows)}; graphs: {len(refvals)}; indices with no problem: {n_ok}/{len(NAMES)}")
    for v in per.values():
        if v["n_problem"]:
            log("  PROBLEM", v["index"], v["problem_statuses"], v["n_problem"], v["first_problem"])
    log(f"max rel diff over all value comparisons: {max(v['max_rel'] for v in per.values()):.3e}")
    log(f"total value comparisons: {sum(v['n_value_comparisons'] for v in per.values())}; "
        f"both-undefined: {sum(v['n_both_undefined'] for v in per.values())}")

    # ---- dependencies ---------------------------------------------------------
    import dependencies
    dep = dependencies.run(refvals, meta, pkgvals, NAMES, outdir, log)

    # ---- cohort-wide checks on package values --------------------------------
    import pandas as pd
    deltas = pd.Series({gid: g.max_degree for gid, g in zip(full.ids, full.graphs)})
    meta_df = fsz.copy()
    meta_df["delta"] = deltas
    ring_bad = 0
    sssr_bad = 0
    for gid, g, mol in zip(full.ids, full.graphs, full.mols):
        cyc = g.m - g.n + 1
        if int(fsz.loc[gid, "n_rings"]) != cyc or not g.is_connected:
            ring_bad += 1
        if mol.GetRingInfo().NumRings() != cyc:
            sssr_bad += 1
    log(f"[panel] n_rings (size_table) != m - n + 1 on {ring_bad} of {len(full.ids)} cohort molecules; "
        f"RDKit SSSR ring count != m - n + 1 on {sssr_bad}")
    dependencies.full_panel(FT, meta_df, NAMES, dep, outdir, log)

    # ---- Hosoya Z beyond 2^53: exact int (direct) vs float64 (index_table) ---
    zrows = []
    for gid, g in zip(full.ids, full.graphs):
        if gid not in big_z:
            continue
        zi = compute_indices(g, ["Z"])[0]["Z"]
        zr = int(refvals[gid]["Z"]) if gid in refvals else ref.hosoya_recursive(molgraph_to_nx(g))
        zt = FT.loc[gid, "Z"]
        zrows.append(dict(id=gid, n=g.n, Z_ref=zr, Z_direct=zi, direct_exact=(zi == zr),
                          Z_table=repr(float(zt)), table_abs_err=int(Fraction(float(zt)) - zr)))
    for nn in (80, 100):
        G = nx.path_graph(nn)
        zr = fib(nn + 1)
        zi = compute_indices(nx_to_molgraph(G), ["Z"])[0]["Z"]
        st1 = prepare_structures([f"P{nn}"], ["C" * nn])
        zt = index_table(st1)[0].loc[f"P{nn}", "Z"]
        zrows.append(dict(id=f"P{nn} (alkane, Z=F_{nn+1})", n=nn, Z_ref=zr, Z_direct=zi, direct_exact=(zi == zr),
                          Z_table=repr(float(zt)), table_abs_err=int(Fraction(float(zt)) - zr)))
    log("[Z>2^53] id n Z_ref Z_direct exact table_float table_err")
    for r in zrows:
        log("   ", r["id"], r["n"], r["Z_ref"], r["Z_direct"], r["direct_exact"], r["Z_table"], r["table_abs_err"])
    with open(os.path.join(outdir, "hosoya_large.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(zrows[0].keys()))
        w.writeheader()
        w.writerows(zrows)

    log(f"elapsed {time.time() - t0:.1f} s")
    with open(os.path.join(outdir, "run_log.txt"), "w") as fh:
        fh.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main(os.path.expanduser(sys.argv[1]), sys.argv[2])
