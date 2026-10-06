"""Head-to-head of the TopoQSPR web API and the live ChemGraphX API on one graph.

Usage: python tools/compare_sites.py [--topoqspr http://127.0.0.1:8000] [--chemgraphx URL] [--smiles SMILES]

Builds the hydrogen-suppressed graph of the molecule with RDKit, sends it to both services in the
format each accepts (TopoQSPR: adjacency matrix text; ChemGraphX: its JSON adjacency list in a
``textareaData`` form field, as its client does), times three requests per endpoint, and compares
the indices both compute under a fixed name map.  Prints a Markdown table; nothing is typed by hand.
ChemGraphX's endpoint URL is the one hard-coded in its client (github.com/aditya-jindal/chemgraphx,
client/src/App.jsx).
"""

from __future__ import annotations

import argparse
import json
import math
import time
import urllib.request
import uuid

CHEMGRAPHX = "https://vpqebjeug2dtqoi3exx6nlfdta0jnxyl.lambda-url.us-east-1.on.aws"
# ChemGraphX display name -> TopoQSPR registry name (degree block, distance block)
NAME_MAP = {
    "First Zagreb": "M1", "Second Zagreb": "M2", "Forgotten": "F", "Hyper Zagreb": "HM",
    "Randic": "R", "Reciprocal Randic": "RR", "Sum Connectivity": "SC", "Harmonic": "H",
    "Geometric Arithmetic": "GA", "Atom Bond Connectivity": "ABC", "Augmented Zagreb": "AZI",
    "Randić": "R", "Somber": "SO", "Symmetric Division Degree": "SDD", "Sigma": "sigma",
    "Wiener Index": "W", "Szeged": "Sz", "Pi": "PI_v", "Edge Pi": "PI", "Mostar": "Mo", "Edge Mostar": "Mo_e",
    "Edge Szeged": "Sz_e", "Edge Wiener": "W_e",
}
# ChemGraphX's edge Wiener counts adjacent edge pairs as distance 0; TopoQSPR uses W(L(G)) (adjacent
# edges at distance 1), so ChemGraphX = TopoQSPR - m(m-1)/2 (documented in tools/crosscheck_chemgraphx.py)
CONVENTION = {"Edge Wiener": lambda ours, m: ours - m * (m - 1) / 2}


def _post_json(url: str, payload: dict, timeout: int = 120) -> tuple[dict, float]:
    body = json.dumps(payload).encode()
    t0 = time.perf_counter()
    with urllib.request.urlopen(urllib.request.Request(url, body, {"Content-Type": "application/json"}), timeout=timeout) as r:
        data = json.load(r)
    return data, time.perf_counter() - t0


def _post_form(url: str, field: str, text: str, timeout: int = 120) -> tuple[dict, float]:
    boundary = uuid.uuid4().hex
    body = (f"--{boundary}\r\nContent-Disposition: form-data; name=\"{field}\"\r\n\r\n{text}\r\n--{boundary}--\r\n").encode()
    req = urllib.request.Request(url, body, {"Content-Type": f"multipart/form-data; boundary={boundary}"})
    t0 = time.perf_counter()
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = json.load(r)
    return data, time.perf_counter() - t0


def graph_inputs(smiles: str) -> tuple[str, str, int, int]:
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    A = Chem.GetAdjacencyMatrix(mol)
    matrix = "\n".join(" ".join(str(int(x)) for x in row) for row in A) + "\n"
    adj = {str(i + 1): [j + 1 for j, x in enumerate(row) if x] for i, row in enumerate(A)}
    return matrix, json.dumps(adj), int(len(A)), int(A.sum() // 2)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--topoqspr", default="http://127.0.0.1:8000")
    ap.add_argument("--chemgraphx", default=CHEMGRAPHX)
    ap.add_argument("--smiles", default="CC(=O)Oc1ccccc1C(=O)O", help="aspirin by default")
    ap.add_argument("--repeats", type=int, default=3)
    a = ap.parse_args()
    matrix, adj_json, n, m = graph_inputs(a.smiles)
    print(f"Graph: {a.smiles}, {n} vertices, {m} edges\n")

    # TopoQSPR: one request returns every index and entropy
    times, ours = [], None
    for _ in range(a.repeats):
        d, dt = _post_json(f"{a.topoqspr}/api/indices", {"format": "adjacency_matrix", "text": matrix, "entropies": True, "view": True})
        times.append(dt)
        ours = dict(zip(d["columns"], d["rows"][0]))
    n_ours = sum(1 for k in ours if k not in ("id", "input_smiles", "status", "reason", "smiles", "inchikey", "kind", "n_vertices", "n_edges", "edge_partition"))
    print(f"TopoQSPR  /api/indices: {n_ours} values in one request; wall {min(times):.2f}-{max(times):.2f} s over {a.repeats} requests")

    # ChemGraphX: three endpoints
    cgx, ctimes = {}, {}
    for ep in ("graph_information", "degree_based", "distance_based"):
        ts = []
        for _ in range(a.repeats):
            d, dt = _post_form(f"{a.chemgraphx}/{ep}", "textareaData", adj_json)
            ts.append(dt)
        ctimes[ep] = (min(ts), max(ts))
        cgx[ep] = d["data"]
    deg = cgx["degree_based"]["four_columns"]
    dist = cgx["distance_based"]
    values = {}
    values.update(deg[" Degree Based Indices"])
    values.update(dist["Distance Indices"])
    n_cgx = sum(len(v) for v in deg.values()) + len(dist["Distance Indices"]) + len(dist["Distance Entropies"])
    for ep, (lo, hi) in ctimes.items():
        print(f"ChemGraphX /{ep}: wall {lo:.2f}-{hi:.2f} s")
    print(f"ChemGraphX total: {n_cgx} values over three requests\n")

    # agreement on the shared, name-mapped indices
    print("| Index (ChemGraphX name) | ChemGraphX | TopoQSPR | agree |\n|---|---|---|---|")
    n_ok = n_cmp = 0
    for cname, oname in NAME_MAP.items():
        if cname not in values or oname not in ours or ours[oname] is None:
            continue
        x, y = float(values[cname]), float(ours[oname])
        y_conv = CONVENTION[cname](y, m) if cname in CONVENTION else y
        ok = math.isclose(x, y_conv, rel_tol=1e-9, abs_tol=1e-9)
        n_cmp += 1
        n_ok += ok
        note = " (after the documented convention shift)" if cname in CONVENTION else ""
        print(f"| {cname} → `{oname}` | {x:.10g} | {y:.10g}{note} | {'yes' if ok else 'NO'} |")
    ent_cgx = deg[" Degree Based Entropies"]
    e_ok = e_cmp = 0
    for cname, oname in NAME_MAP.items():
        if cname in ent_cgx and f"ENT_{oname}" in ours and ours[f"ENT_{oname}"] is not None:
            e_cmp += 1
            e_ok += math.isclose(float(ent_cgx[cname]), float(ours[f"ENT_{oname}"]), rel_tol=1e-9, abs_tol=1e-9)
    print(f"\nShared indices compared: {n_cmp}, agreeing: {n_ok}; shared degree entropies compared: {e_cmp}, agreeing: {e_ok}")
    print(f"Graph information agrees: vertices {cgx['graph_information']['Vertices Count']} = {int(ours['n_vertices'])}, "
          f"edges {cgx['graph_information']['Edge Count']} = {int(ours['n_edges'])}")


if __name__ == "__main__":
    main()
