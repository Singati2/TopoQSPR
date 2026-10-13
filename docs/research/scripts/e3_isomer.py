"""AMD-2026-10-13-E3 (exploratory): see docs/research/06_PROSPECTIVE_STUDY_AMENDMENT.md.

    python -I docs/research/scripts/e3_isomer.py --smoke /tmp/e3_smoke
    python -I docs/research/scripts/e3_isomer.py --seeds 0,1,2,3,4 ~/topoqspr_runs/research_2026-10-13_E3
    python -I docs/research/scripts/e3_isomer.py --negative-control ~/topoqspr_runs/research_2026-10-13_E3
"""
import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
import sklearn
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from topoqspr import evaluate as ev
from topoqspr import expdata
from topoqspr import features as ft
from topoqspr import stats
from topoqspr.experiments import index_table, prepare_structures

RDLogger.DisableLog("rdApp.*")
FAMILIES = ["ridge", "elasticnet", "pls", "knn", "svr", "extratrees"]
CELLS = [("mp_bradley", "E3a"), ("mp_bradley", "E3b"), ("logp_opera", "E3a"), ("logp_opera", "E3b")]
CONTRASTS = [("TIa_composition", "composition"), ("ECFP_TIa", "ECFP"),
             ("rdkit2d_asinh_TIa", "rdkit2d_asinh"), ("ECFP", "TIa_composition"),
             ("rdkit2d_asinh", "TIa_composition"),
             ("distance_composition", "composition"), ("edgetypes_composition", "composition")]
MARGIN = 0.05


def hsh(df: pd.DataFrame) -> str:
    return hashlib.sha256(np.ascontiguousarray(df.to_numpy(dtype=float)).tobytes()
                          + "|".join(map(str, df.columns)).encode()).hexdigest()


REPO = Path(__file__).resolve().parents[3]
DIST = [c.split(".", 1)[1] for c in json.loads(
    (REPO / "results/expB/exploratory_families/study.json").read_text())["feature_sets"]["TI_distance"]]


def parents(endpoint: str, drop_disagreeing: bool = True):
    """Standardised parents with a usable measurement: (ids, smiles, y, agg)."""
    agg = expdata.get_target(endpoint, offline=True, drop_disagreeing=drop_disagreeing)
    recs = expdata.load_records(endpoint, offline=True)[["record_id", "inchikey", "canonical_smiles", "status"]]
    recs = recs[recs.inchikey.isin(set(agg.inchikey)) & (recs.status == "ok")]
    return recs, agg


def physprop_free_bradley_keys() -> set:
    x = pd.read_excel(Path("~/.cache/topoqspr/expdata/bradley/BradleyMeltingPointDataset.xlsx").expanduser())
    phys = set(x.loc[x.source.astype(str).str.contains("PHYSPROP", case=False), "key"].astype(str))
    recs = expdata.load_records("mp_bradley", offline=True)[["record_id", "inchikey", "status"]]
    recs = recs[recs.status == "ok"]
    has = recs.assign(p=recs.record_id.astype(str).isin(phys)).groupby("inchikey").p.any()
    return set(has[~has].index)


def featurise(keys, smiles):
    st = prepare_structures(list(keys), list(smiles))
    ids = list(st.ok.index)
    mols = dict(zip(st.ids, st.mols))
    idx = pd.Index(ids, name="id")
    kmols = [mols[i] for i in ids]
    ti, _ = index_table(st)
    ti = ti.loc[ids]
    size = ft.size_features(kmols, idx)
    comp = ft.composition_features(kmols, idx)
    rd = ft.rdkit2d_features(kmols, None, idx)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    E = pd.DataFrame(np.vstack([gen.GetCountFingerprintAsNumPy(m).astype(float) for m in kmols]),
                     index=idx, columns=[f"ecfp.{i}" for i in range(2048)])

    def etc(m):
        d = [a.GetDegree() for a in m.GetAtoms()]
        c = {f"e{i}{j}": 0 for i in range(1, 5) for j in range(i, 5) if (i, j) != (1, 1)}
        for b in m.GetBonds():
            i, j = sorted((d[b.GetBeginAtomIdx()], d[b.GetEndAtomIdx()]))
            if f"e{i}{j}" in c:
                c[f"e{i}{j}"] += 1
        c["maxdeg"] = max(d)
        return c
    et = pd.DataFrame([etc(m) for m in kmols], index=idx).astype(float)
    tia = np.arcsinh(ti).add_prefix("ti.")
    rda = np.arcsinh(rd).add_prefix("rd.")
    sets = {"size": size, "composition": comp, "TIa": tia,
            "TIa_composition": pd.concat([tia, comp], axis=1),
            "distance_composition": pd.concat([np.arcsinh(ti[DIST]).add_prefix("dist."), comp], axis=1),
            "edgetypes_composition": pd.concat([et, comp], axis=1),
            "rdkit2d_asinh": rda, "rdkit2d_asinh_TIa": pd.concat([rda, tia], axis=1),
            "ECFP": E, "ECFP_TIa": pd.concat([E, tia], axis=1)}
    meta = pd.DataFrame({"formula": [CalcMolFormula(m) for m in kmols],
                         "constitution": [i[:14] for i in ids]}, index=idx)
    return sets, meta, list(ti.columns)


def include_ge5(meta: pd.DataFrame) -> pd.Index:
    n_cons = meta.drop_duplicates("constitution").groupby("formula").size()
    return meta.index[meta.formula.map(n_cons) >= 5]


def build(endpoint: str, drop_disagreeing: bool = True, restrict: set | None = None):
    """Inclusion uses structures only (plus the package's pre-specified replicate filter);
    y is attached last."""
    recs, agg = parents(endpoint, drop_disagreeing)
    if restrict is not None:
        recs = recs[recs.inchikey.isin(restrict)]
    rep = recs.drop_duplicates("inchikey")
    sets, meta, reg = featurise(rep.inchikey, rep.canonical_smiles)
    keep = include_ge5(meta)
    sets = {k: v.loc[keep] for k, v in sets.items()}
    meta = meta.loc[keep]
    info = dict(endpoint=endpoint, drop_disagreeing=drop_disagreeing, measured_parents=len(agg),
                candidates=len(rep), included=len(keep),
                formula_groups=int(meta.formula.nunique()), registry=reg,
                registry_sha256=hashlib.sha256("|".join(reg).encode()).hexdigest(),
                feature_sha256={k: hsh(v) for k, v in sets.items()})
    y = agg.set_index("inchikey").loc[keep, "value"].astype(float)  # attached last
    return sets, meta, y, info


def paired(res, a, b, alpha=0.05):
    fa = res[a].folds.set_index(["repeat", "fold"])
    fb = res[b].folds.set_index(["repeat", "fold"])
    assert fa[["n_train", "n_test"]].equals(fb[["n_train", "n_test"]]), (a, b)
    out = dict(set=a, reference=b)
    for k, d in (("d_r2", fa.r2 - fb.r2), ("d_mae", fb.mae - fa.mae)):
        t = stats.corrected_ttest(d, fa.n_train, fa.n_test, alpha=alpha)
        t90 = stats.corrected_ttest(d, fa.n_train, fa.n_test, alpha=0.10)
        out.update({f"{k}_mean": t["mean"], f"{k}_ci_low": t["ci_low"], f"{k}_ci_high": t["ci_high"],
                    f"{k}_p": t["p"], f"{k}_ci90_low": t90["ci_low"], f"{k}_ci90_high": t90["ci_high"],
                    f"{k}_n_folds_positive": int((d > 0).sum())})
    return out


def within_formula_r2(pred: pd.DataFrame, formula: pd.Series) -> float:
    vals = []
    for _, d in pred.groupby("repeat"):
        d = d.assign(f=formula.loc[d.row].values)
        d = d[d.groupby(["fold", "f"]).row.transform("size") >= 2]
        den = ((d.y_true - d.groupby(["fold", "f"]).y_true.transform("mean")) ** 2).sum()
        vals.append(1 - ((d.y_true - d.y_pred) ** 2).sum() / den if den > 0 else np.nan)
    return float(np.mean(vals))


def run_cell(sets, meta, y, estimand, seeds, n_jobs, set_names=None):
    groups = (meta.constitution if estimand == "E3a" else meta.formula).to_numpy(dtype=object)
    names = set_names or list(sets)
    res = {k: ev.nested_cv(sets[k], y, FAMILIES, groups=groups, seeds=seeds, n_jobs=n_jobs)
           for k in names}
    # integrity: identical folds, no group straddles folds
    ref = res[names[0]].predictions
    for k in names[1:]:
        p = res[k].predictions
        assert (p.sort_values(["repeat", "row"]).fold.values ==
                ref.sort_values(["repeat", "row"]).fold.values).all()
    g = pd.Series(groups, index=y.index)
    straddle = int(ref.assign(g=g.loc[ref.row].values).groupby(["repeat", "g"]).fold.nunique().gt(1).sum())
    assert straddle == 0, straddle
    return res


def external(out: Path, n_jobs: int) -> None:
    """06b R3: train on PHYSPROP-free Bradley isomer set; test on OPERA (PHYSPROP) parents."""
    from joblib import Parallel, delayed
    train_keys = physprop_free_bradley_keys()
    sets, meta, y, info = build("mp_bradley", restrict=train_keys)
    recs, agg_o = parents("mp_opera")
    rep = recs.drop_duplicates("inchikey")
    rep = rep[~rep.inchikey.str[:14].isin(set(meta.constitution))]
    tsets, tmeta, _ = featurise(rep.inchikey, rep.canonical_smiles)
    known = tmeta.formula.isin(set(meta.formula))
    tests = {"known_formula": tmeta.index[known], "unseen_formula": tmeta.index[~known]}
    yt = agg_o.set_index("inchikey")["value"].astype(float)
    groups = meta.constitution.to_numpy(dtype=object)

    def fit(k):
        sel = ev.select_and_refit(sets[k].to_numpy(float), y.to_numpy(), FAMILIES, inner=3,
                                  random_state=0, groups=groups)
        cols = list(sets[k].columns)
        X = tsets[k].reindex(columns=cols).to_numpy(float)
        return k, sel.family, json.dumps(sel.params, default=str), sel.estimator.predict(X)
    fits = Parallel(n_jobs=min(n_jobs, len(sets)))(delayed(fit)(k) for k in sets)
    pred = pd.DataFrame({k: p for k, _, _, p in fits}, index=tmeta.index)
    fam = {k: (f, prm) for k, f, prm, _ in fits}
    rng = np.random.default_rng(2026)
    rows, prow = [], []
    for name, ids in tests.items():
        yy, P, F = yt.loc[ids].to_numpy(), pred.loc[ids], tmeta.loc[ids, "formula"].to_numpy()
        uf = np.unique(F); members = [np.flatnonzero(F == f) for f in uf]
        boots = [np.concatenate([members[j] for j in rng.integers(0, uf.size, uf.size)]) for _ in range(2000)]
        def r2(a, b): return 1 - ((a - b) ** 2).sum() / ((a - a.mean()) ** 2).sum()
        for k in sets:
            rows.append(dict(test=name, set=k, n=len(ids), n_formulas=int(uf.size), family=fam[k][0],
                             r2=r2(yy, P[k].to_numpy()), mae=float(np.abs(yy - P[k]).mean())))
        for a, b in CONTRASTS:
            pa, pb = P[a].to_numpy(), P[b].to_numpy()
            d_r2 = r2(yy, pa) - r2(yy, pb); d_mae = np.abs(yy - pb).mean() - np.abs(yy - pa).mean()
            br = np.array([r2(yy[i], pa[i]) - r2(yy[i], pb[i]) for i in boots])
            bm = np.array([np.abs(yy[i] - pb[i]).mean() - np.abs(yy[i] - pa[i]).mean() for i in boots])
            prow.append(dict(test=name, set=a, reference=b, d_r2=d_r2,
                             d_r2_ci_low=np.percentile(br, 2.5), d_r2_ci_high=np.percentile(br, 97.5),
                             d_mae=d_mae, d_mae_ci_low=np.percentile(bm, 2.5), d_mae_ci_high=np.percentile(bm, 97.5),
                             boot_p_two_sided=float(min(1, 2 * min((br <= 0).mean(), (br >= 0).mean())))))
    d = out / "external"
    d.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(d / "summary.csv", index=False)
    pd.DataFrame(prow).to_csv(d / "paired.csv", index=False)
    pred.assign(y=yt.loc[pred.index], formula=tmeta.formula).to_csv(d / "predictions.csv.gz")
    (d / "info.json").write_text(json.dumps(dict(train={k: v for k, v in info.items() if k != "registry"},
                                                 train_physprop_free=True,
                                                 test_sizes={k: len(v) for k, v in tests.items()}), indent=1))
    print("EXTERNAL done", {k: len(v) for k, v in tests.items()}, flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--negative-control", action="store_true")
    ap.add_argument("--external", action="store_true")
    ap.add_argument("--sensitivity-noflt", action="store_true")
    ap.add_argument("--cells", default="mp_bradley:E3a,mp_bradley:E3b,logp_opera:E3a,logp_opera:E3b")
    ap.add_argument("--seeds", default="0,1,2,3,4")
    ap.add_argument("--n-jobs", type=int, default=14)
    a = ap.parse_args()
    out = Path(a.out).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    seeds = tuple(int(s) for s in a.seeds.split(","))
    t_start = time.time()
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                            cwd=Path(__file__).parent).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"],
                                capture_output=True, text=True, cwd=Path(__file__).parent).stdout.strip())
    info_all = dict(amendment="AMD-2026-10-13-E3", commit_at_start=commit, dirty_at_start=dirty,
                    python=platform.python_version(), rdkit=rdkit.__version__,
                    sklearn=sklearn.__version__, seeds=list(seeds), cells={})
    built = {}

    def get(endpoint):
        if endpoint not in built:
            built[endpoint] = build(endpoint)
            sets, meta, y, info = built[endpoint]
            meta.to_csv(out / f"inclusion_{endpoint}.csv")
            print(f"built {endpoint}: {json.dumps({k: v for k, v in info.items() if k != 'registry'})[:300]}", flush=True)
        return built[endpoint]

    if a.smoke:
        sets, meta, y, info = get("mp_bradley")
        y = pd.Series(np.random.default_rng(12345).permutation(y.to_numpy()), index=y.index)
        t0 = time.time()
        run_cell(sets, meta, y, "E3a", (0,), a.n_jobs)
        el = time.time() - t0
        n_total = sum(len(get(e)[2]) for e in ("mp_bradley", "logp_opera"))
        proj = el * 5 * 2 * n_total / len(y)  # 5 seeds, 2 estimands, both endpoints (size-linear guess)
        print(f"SMOKE done: one seed, MP E3a, {el:.0f}s; projected 5-seed total ~{proj / 3600:.1f} h "
              f"(linear in n; numbers discarded)", flush=True)
        (out / "smoke_timing.json").write_text(json.dumps(dict(seconds_one_seed_mp_e3a=el,
                                                               projected_hours_5_seeds=proj / 3600)))
        return

    if a.external:
        external(out, a.n_jobs)
        return

    if a.sensitivity_noflt:
        sets, meta, y, info = build("mp_bradley", drop_disagreeing=False)
        res = run_cell(sets, meta, y, "E3a", (0,), a.n_jobs, ["composition", "TIa_composition"])
        row = paired(res, "TIa_composition", "composition")
        pd.DataFrame([dict(included=len(y), **row)]).to_csv(out / "sensitivity_noflt_paired.csv", index=False)
        print("SENS", json.dumps({k: round(v, 4) if isinstance(v, float) else v for k, v in row.items()}))
        return

    if a.negative_control:
        sets, meta, y, info = get("mp_bradley")
        rng = np.random.default_rng(777)
        yp = y.copy()
        for f, idx in meta.groupby("formula").groups.items():
            yp.loc[idx] = rng.permutation(y.loc[idx].to_numpy())
        res = run_cell(sets, meta, yp, "E3a", (0,), a.n_jobs, ["composition", "TIa_composition"])
        row = paired(res, "TIa_composition", "composition")
        pd.DataFrame([row]).to_csv(out / "negative_control_paired.csv", index=False)
        print("NEGCTRL", json.dumps({k: round(v, 4) if isinstance(v, float) else v for k, v in row.items()}))
        return

    for endpoint, estimand in [tuple(c.split(":")) for c in a.cells.split(",")]:
        sets, meta, y, info = get(endpoint)
        info_all["cells"].setdefault(endpoint, {k: v for k, v in info.items() if k != "registry"})
        cdir = out / "cells" / f"{endpoint}_{estimand}"
        cdir.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        res = run_cell(sets, meta, y, estimand, seeds, a.n_jobs)
        rows = []
        for k, r in res.items():
            s = r.summary()
            rows.append(dict(endpoint=endpoint, estimand=estimand, set=k, n=len(y), p=r.meta["p"],
                             r2_fold_mean=s["r2_fold_mean"], r2_pooled_mean=s["r2_pooled_mean"],
                             mae_fold_mean=s["mae_fold_mean"], min_fold_r2=float(r.folds.r2.min()),
                             within_formula_r2=within_formula_r2(r.predictions, meta.formula),
                             n_failed_folds=s["n_failed_folds"],
                             family_counts=json.dumps(s["family_counts"]), reason=s["reason"]))
        pd.DataFrame(rows).to_csv(cdir / "summary.csv", index=False)
        pd.DataFrame([dict(endpoint=endpoint, estimand=estimand, **paired(res, x, ref))
                      for x, ref in CONTRASTS]).to_csv(cdir / "paired.csv", index=False)
        pd.concat([r.folds.assign(set=k) for k, r in res.items()]).to_csv(cdir / "folds.csv", index=False)
        pd.concat([r.predictions.assign(set=k) for k, r in res.items()]).to_csv(
            cdir / "predictions.csv.gz", index=False)
        print(f"{endpoint} {estimand} n={len(y)} {time.time() - t0:.0f}s", flush=True)
        info_all["wall_s"] = round(time.time() - t_start, 1)
        (out / "run_info.json").write_text(json.dumps(info_all, indent=1) + "\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
