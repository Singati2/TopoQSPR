"""Algebraic-dependency analysis on the reference values (called by verify_indices.py)."""
from __future__ import annotations

import json
import os
from collections import defaultdict
from decimal import Decimal
from fractions import Fraction

import numpy as np

import refimpl as ref
from refimpl import D, Undefined

TOL_DEC = Decimal("1e-35")


def _is_def(v):
    return not isinstance(v, Undefined)


def _eq(a, b):
    if isinstance(a, Fraction) and isinstance(b, Fraction):
        return a == b
    a, b = D(a), D(b)
    return abs(a - b) <= Decimal("1e-40") * max(Decimal(1), abs(a), abs(b))


def rref(M, field):
    """Column-ordered RREF. field = 'Q' (Fractions, exact) or 'R' (Decimals, tolerance).
    Returns (pivot_cols, reduced matrix)."""
    A = [list(r) for r in M]
    nr, nc = len(A), len(A[0])
    zero = (lambda x: x == 0) if field == "Q" else (lambda x: abs(x) < TOL_DEC)
    piv, r = [], 0
    for j in range(nc):
        if r == nr:
            break
        best = max(range(r, nr), key=lambda i: abs(A[i][j]))
        if zero(A[best][j]):
            for i in range(r, nr):  # clean numerical dust
                if field == "R":
                    A[i][j] = A[i][j] * 0
            continue
        A[r], A[best] = A[best], A[r]
        p = A[r][j]
        A[r] = [x / p for x in A[r]]
        for i in range(nr):
            if i != r and not zero(A[i][j]) if field == "R" else (i != r and A[i][j] != 0):
                f = A[i][j]
                Ai, Ar = A[i], A[r]
                for k in range(j, nc):
                    Ai[k] -= f * Ar[k]
        piv.append(j)
        r += 1
    return piv, A


def _nice(c):
    if isinstance(c, Fraction):
        return str(c)
    f = Fraction(str(c)).limit_denominator(10000)
    if abs(D(f) - c) < Decimal("1e-30"):
        return str(f)
    return f"{c:.15g}"


def relations_from(piv, A, cols, scale=None):
    rels = {}
    pivset = set(piv)
    for j, name in enumerate(cols):
        if j in pivset:
            continue
        terms = []
        for k, pj in enumerate(piv):
            if pj > j:
                break
            a = A[k][j]
            if scale is not None:
                a = a * scale[j] / scale[pj]
            if (a != 0) if isinstance(a, Fraction) else (abs(a) > Decimal("1e-30")):
                terms.append((cols[pj], a))
        rels[name] = terms
    return rels


def check_relation(name, terms, vals_by_graph):
    """Return (n_checked, n_fail, worst_rel, failing witness)."""
    nchk = nfail = 0
    worst = Decimal(0)
    wit = None
    for gid, vals in vals_by_graph.items():
        if not _is_def(vals.get(name)) or any(not _is_def(vals.get(t)) for t, _ in terms):
            continue
        lhs = vals[name]
        exact = isinstance(lhs, Fraction) and all(isinstance(vals[t], Fraction) and isinstance(c, Fraction)
                                                  for t, c in terms)
        if exact:
            res = lhs - sum((c * vals[t] for t, c in terms), Fraction(0))
            bad = res != 0
            rel = Decimal(0) if not bad else Decimal(1)
        else:
            rhs = sum((D(c) * D(vals[t]) for t, c in terms), Decimal(0))
            scale = max(Decimal(1), abs(D(lhs)), sum((abs(D(c) * D(vals[t])) for t, c in terms), Decimal(0)))
            rel = abs(D(lhs) - rhs) / scale
            bad = rel > Decimal("1e-25")
        nchk += 1
        if bad:
            nfail += 1
            if wit is None:
                wit = gid
        worst = max(worst, rel)
    return nchk, nfail, float(worst), wit


def run(refvals, meta, pkgvals, NAMES, outdir, log):
    out = {}
    base = {}
    for gid, v in refvals.items():
        mt = meta[gid]
        vv = dict(v)
        vv["const"] = Fraction(1)
        vv["n"] = Fraction(mt["n"])
        vv["m"] = Fraction(mt["m"])
        base[gid] = vv
    connected = [g for g in refvals if meta[g]["connected"] and meta[g]["n"] >= 2]
    atlas = [g for g in connected if meta[g]["group"] == "atlas"]
    others = {g: base[g] for g in connected if meta[g]["group"] != "atlas"}

    # ---- 1. functions of (n), (n, m), (n, m, degree sequence) ----------------
    keyfns = {"n": lambda mt: (mt["n"],), "n,m": lambda mt: (mt["n"], mt["m"]),
              "n,m,degseq": lambda mt: (mt["n"], mt["m"], tuple(mt["degseq"]))}
    func_of = {}
    for label, kf in keyfns.items():
        groups = defaultdict(list)
        for g in connected:
            groups[kf(meta[g])].append(g)
        nontriv = sum(1 for v in groups.values() if len(v) > 1)
        hits = []
        for nm in NAMES:
            ok = True
            for gs in groups.values():
                vals = [base[g][nm] for g in gs if _is_def(base[g][nm])]
                if any(not _eq(vals[0], x) for x in vals[1:]):
                    ok = False
                    break
            if ok:
                hits.append(nm)
        func_of[label] = dict(nontrivial_groups=nontriv, indices=hits)
        log(f"[dep] indices that are exact functions of ({label}) on {len(connected)} connected graphs "
            f"({nontriv} groups with >1 graph): {hits}")
    out["function_of"] = func_of

    # ---- 2. exact rational RREF (atlas, rows with every index defined) --------
    rows = [g for g in atlas if all(_is_def(base[g][nm]) for nm in NAMES)]
    log(f"[dep] atlas rows with all 104 indices defined (Delta >= 3): {len(rows)} of {len(atlas)}")
    rational = [nm for nm in NAMES if all(isinstance(base[g][nm], Fraction) for g in connected
                                          if _is_def(base[g][nm]))]
    cols_q = ["const", "n", "m"] + rational
    Mq = [[base[g][c] for c in cols_q] for g in rows]
    piv, A = rref(Mq, "Q")
    rels_q = relations_from(piv, A, cols_q)
    log(f"[dep] rational-valued indices: {len(rational)}; exact rank of [1,n,m,rational] over atlas: "
        f"{len(piv)} of {len(cols_q)} columns")
    qres = {}
    regime = {g: v for g, v in others.items() if 3 <= meta[g]["delta"] <= 4}
    for nm, terms in rels_q.items():
        nchk, nfail, worst, wit = check_relation(nm, terms, regime)
        uchk, ufail, _, uwit = check_relation(nm, terms, others)
        s = " + ".join(f"({_nice(c)})*{t}" for t, c in terms) or "0"
        status = (f"holds on {nchk} non-atlas graphs with 3<=Delta<=4" if nfail == 0
                  else f"FAILS on {nfail}/{nchk} non-atlas graphs with 3<=Delta<=4 (e.g. {wit})")
        status += "; UNIVERSAL on all connected test graphs" if ufail == 0 else f"; not universal (e.g. {uwit})"
        qres[nm] = dict(relation=f"{nm} = {s}", checked_in_regime=nchk, failures_in_regime=nfail,
                        witness=wit, universal=ufail == 0, universal_witness=uwit)
        log(f"[dep:Q] {nm} = {s}   [{status}]")
    out["rational_rref"] = dict(columns=cols_q, rank=len(piv), relations=qres)

    # ---- 3. real RREF over all indices (60-digit Decimals) -------------------
    def real_rref(rowset, label):
        cols = ["const", "n", "m"] + NAMES
        M = [[D(base[g][c]) for c in cols] for g in rowset]
        scale = [max(abs(M[i][j]) for i in range(len(M))) or Decimal(1) for j in range(len(cols))]
        Mn = [[M[i][j] / scale[j] for j in range(len(cols))] for i in range(len(M))]
        piv, A = rref(Mn, "R")
        rels = relations_from(piv, A, cols, scale)
        res = {}
        nhold = 0
        for nm, terms in rels.items():
            other = {g: v for g, v in others.items()
                     if (3 <= meta[g]["delta"] <= 4 if label == "all"
                         else meta[g]["delta"] == int(label.split("=")[1]))}
            nchk, nfail, worst, wit = check_relation(nm, terms, other)
            uchk, ufail, _, uwit = check_relation(nm, terms, others)
            nhold += nfail == 0
            res[nm] = dict(depends_on=[t for t, _ in terms], coefficients=[str(c) for _, c in terms],
                           relation=f"{nm} = " + " + ".join(f"({_nice(c)})*{t}" for t, c in terms),
                           checked_outside_atlas=nchk, failures=nfail, worst_rel=worst, witness=wit,
                           universal=ufail == 0, universal_witness=uwit)
        log(f"[dep:R {label}] rows {len(rowset)}; numerical rank of [1,n,m,104 indices] = {len(piv)} "
            f"of {len(cols)}; dependent columns {len(rels)}, of which {nhold} hold on every non-atlas "
            f"test graph of that stratum ({'3<=Delta<=4' if label == 'all' else label})")
        for nm, r in res.items():
            flag = "holds" if r["failures"] == 0 else f"FAILS {r['failures']}/{r['checked_outside_atlas']} e.g. {r['witness']}"
            flag += "; universal" if r["universal"] else f"; not universal, e.g. {r['universal_witness']}"
            if len(r["depends_on"]) <= 6:
                log(f"   {r['relation']}  [{flag}, {r['checked_outside_atlas']} graphs]")
            else:
                log(f"   {nm} in span of {len(r['depends_on'])} earlier columns  [{flag}, "
                    f"{r['checked_outside_atlas']} graphs]")
        return dict(rank=len(piv), pivots=[cols[j] for j in piv], relations=res)

    out["real_rref_all"] = real_rref(rows, "all")
    d4 = [g for g in rows if meta[g]["delta"] == 4]
    out["real_rref_delta4"] = real_rref(d4, "delta=4")

    # ---- 4. degree-family coefficient rank over the 10 edge types -------------
    types = [(a, b) for a in range(1, 5) for b in range(a, 5)]
    fam = {}
    for key in ref.DEGREE_KEYS:
        fam[key] = {t: t for t in types}
    def coef_row(key, wmap):
        row = []
        for a, b in types:
            try:
                v = ref.PHI[key](Fraction(wmap[a]), Fraction(wmap[b]))
                row.append(D(v))
            except TypeError:
                row.append(ref.PHI[key](D(wmap[a]), D(wmap[b])))
            except ZeroDivisionError:
                row.append(None)
        return row
    ident = {d: d for d in range(1, 5)}
    deg_rows = {k: coef_row(k, ident) for k in ref.DEGREE_KEYS}
    rk_rows = {}
    for k in (1, 2, 3):
        wm = {d: ref.reverse_degree(d, k, 4) for d in range(1, 5)}
        for key in ref.RK_KEYS:
            rk_rows[f"R{k}_{key}"] = coef_row(key, wm)
    rk_rows["HCM2"] = coef_row("HCM2", {d: 4 - d + 1 for d in range(1, 5)})
    # AZI is undefined on type (1,1) (only K2); drop that type for the rank
    tt = [i for i, t in enumerate(types) if t != (1, 1)]
    def rank_of(rowdict):
        M = [[r[i] for i in tt] for r in rowdict.values()]
        # rank of the transpose (columns = indices) equals rank
        Mt = [list(x) for x in zip(*M)]
        p, _ = rref(Mt, "R")
        return len(p)
    r_deg = rank_of(deg_rows)
    r_all = rank_of({**deg_rows, **rk_rows})
    nm_rows = {"n": [D(Fraction(1, a) + Fraction(1, b)) for a, b in types],
               "m": [Decimal(1)] * len(types)}
    r_deg_nm = rank_of({**deg_rows, **nm_rows})
    log(f"[dep] Delta<=4 degree-based indices are sum_(a<=b) N_ab phi(a,b) over 9 edge types "
        f"(1,2)..(4,4) (excluding (1,1) = K2). Exact-coefficient rank: 32 degree indices -> {r_deg}; "
        f"with n and m -> {r_deg_nm}; 32 degree + 48 R_k + HCM2 (Delta = 4 weights) -> {r_all}")
    out["edge_type_rank"] = dict(types=[f"{a}{b}" for a, b in types if (a, b) != (1, 1)],
                                 degree=r_deg, degree_plus_n_m=r_deg_nm, degree_rk_hcm2_delta4=r_all)

    # package values on the drug panel: numerical rank (float64 SVD) of the degree family
    drug_ids = [g for g in connected if meta[g]["group"] == "drug" and meta[g]["delta"] == 4]
    fam_names = list(ref.DEGREE_KEYS) + [f"R{k}_{key}" for k in (1, 2, 3) for key in ref.RK_KEYS] + ["HCM2"]
    X = np.array([[float(pkgvals[(g, "smiles")][nm]) for nm in fam_names] for g in drug_ids])
    Xs = X / np.abs(X).max(axis=0)
    sv = np.linalg.svd(Xs, compute_uv=False)
    num_rank = int((sv > sv[0] * 1e-10).sum())
    log(f"[dep] package float values, {len(drug_ids)} Delta=4 drugs x {len(fam_names)} degree/R_k/HCM2 "
        f"indices: numerical rank {num_rank}; singular values (scaled) {np.array2string(sv / sv[0], precision=2)}")
    out["drug_panel_numerical_rank"] = dict(n_drugs=len(drug_ids), n_indices=len(fam_names), rank=num_rank,
                                            sv=[float(x) for x in sv / sv[0]])

    # ---- 5. pairwise exact equality / proportionality -------------------------
    pairs = []
    for i, a in enumerate(NAMES):
        for b in NAMES[i + 1:]:
            ratio = None
            ok = True
            nn = 0
            for g in connected:
                va, vb = base[g][a], base[g][b]
                if not (_is_def(va) and _is_def(vb)):
                    continue
                nn += 1
                if D(vb) == 0:
                    if D(va) != 0:
                        ok = False
                        break
                    continue
                r = D(va) / D(vb)
                if ratio is None:
                    ratio = r
                elif abs(r - ratio) > Decimal("1e-35") * max(Decimal(1), abs(ratio)):
                    ok = False
                    break
            if ok and nn > 0:
                pairs.append((a, b, str(ratio)))
    log(f"[dep] exactly proportional index pairs on all {len(connected)} connected test graphs: {pairs}")
    out["proportional_pairs"] = pairs

    with open(os.path.join(outdir, "dependencies.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return out


def full_panel(T, meta_df, NAMES, dep, outdir, log):
    """Checks on the package's float values for every accepted molecule of the cohort."""
    res = {}
    n = meta_df["n_heavy"].astype(int)
    m = meta_df["n_bonds"].astype(int)
    res["ReZ1_equals_n"] = int(((T["ReZ1"] - n).abs() <= 1e-9 * n).sum())
    res["n_molecules"] = len(T)
    log(f"[panel] {len(T)} cohort molecules: ReZ1 == n_heavy (|diff| <= 1e-9 n) for {res['ReZ1_equals_n']}")
    fam = list(ref.DEGREE_KEYS) + [f"R{k}_{key}" for k in (1, 2, 3) for key in ref.RK_KEYS] + ["HCM2"]
    for delta in (3, 4):
        sub = T[meta_df["delta"] == delta]
        cols = [c for c in fam if sub[c].notna().all()]
        X = sub[cols].to_numpy(float)
        X = X / np.abs(X).max(axis=0)
        sv = np.linalg.svd(X, compute_uv=False)
        rk = int((sv > sv[0] * 1e-10).sum())
        gap = float(sv[rk] / sv[0]) if rk < len(sv) else 0.0
        res[f"rank_delta{delta}"] = dict(n=len(sub), n_indices=len(cols), rank=rk,
                                         last_kept=float(sv[rk - 1] / sv[0]), first_dropped=gap)
        log(f"[panel] Delta={delta}: {len(sub)} molecules x {len(cols)} degree/R_k/HCM2 indices: numerical rank "
            f"{rk} (sv[{rk-1}]/sv[0] = {sv[rk-1]/sv[0]:.2e}, sv[{rk}]/sv[0] = {gap:.2e})")
    # relations found on the atlas, checked on the package's floats
    for label, key, cond in (("3<=Delta<=4", "real_rref_all", meta_df["delta"].between(3, 4)),
                             ("Delta=4", "real_rref_delta4", meta_df["delta"] == 4)):
        rels = dep[key]["relations"]
        sub = T[cond].copy()
        sub["const"] = 1.0
        sub["n"] = n[cond].astype(float)
        sub["m"] = m[cond].astype(float)
        nfail_rel = 0
        worst_all = 0.0
        for nm, r in rels.items():
            if r["failures"]:
                continue
            terms = list(zip(r["depends_on"], (float(Decimal(c)) for c in r["coefficients"])))
            ok = sub[[nm] + [t for t, _ in terms]].notna().all(axis=1)
            lhs = sub.loc[ok, nm].to_numpy(float)
            rhs = sum(c * sub.loc[ok, t].to_numpy(float) for t, c in terms)
            scale = np.maximum(1.0, sum(np.abs(c * sub.loc[ok, t].to_numpy(float)) for t, c in terms))
            rel = np.abs(lhs - rhs) / scale
            worst = float(rel.max()) if len(rel) else 0.0
            worst_all = max(worst_all, worst)
            if worst > 1e-9:
                nfail_rel += 1
                log(f"[panel]   relation for {nm} ({label}) violated on package values: worst rel {worst:.2e}")
        log(f"[panel] {label}: {sum(1 for r in rels.values() if not r['failures'])} atlas relations re-checked on "
            f"{int(cond.sum())} cohort molecules (package floats): {nfail_rel} violated; worst rel residual "
            f"{worst_all:.2e}")
        res[f"relations_{key}"] = dict(n_molecules=int(cond.sum()), violated=nfail_rel, worst_rel=worst_all)
    with open(os.path.join(outdir, "full_panel.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    return res
