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
             ("rdkit2d_asinh", "TIa_composition")]
MARGIN = 0.05


def hsh(df: pd.DataFrame) -> str:
    return hashlib.sha256(np.ascontiguousarray(df.to_numpy(dtype=float)).tobytes()
                          + "|".join(map(str, df.columns)).encode()).hexdigest()


def build(endpoint: str):
    """Target-blind inclusion, then features; y is attached last."""
    agg = expdata.get_target(endpoint, offline=True)
    recs = expdata.load_records(endpoint, offline=True)[["inchikey", "canonical_smiles", "status"]]
    recs = recs[recs.inchikey.isin(set(agg.inchikey)) & (recs.status == "ok")].drop_duplicates("inchikey")
    st = prepare_structures(recs.inchikey.tolist(), recs.canonical_smiles.tolist())
    ok = st.ok
    ids = list(ok.index)
    mols = {i: m for i, m in zip(st.ids, st.mols)}
    form = pd.Series({i: CalcMolFormula(mols[i]) for i in ids})
    cons = pd.Series({i: i[:14] for i in ids})
    n_cons = pd.DataFrame({"f": form, "c": cons}).drop_duplicates("c").groupby("f").size()
    keep = [i for i in ids if n_cons[form[i]] >= 5]
    info = dict(endpoint=endpoint, measured_parents=len(agg), standardised_ok=len(ids),
                rejected=int(len(recs) - len(ids)), included=len(keep),
                formula_groups=int((n_cons >= 5).sum()))
    st_k = prepare_structures(keep, [Chem.MolToSmiles(mols[i]) for i in keep])
    assert list(st_k.ok.index) == keep
    idx = pd.Index(keep, name="id")
    ti, _ = index_table(st_k)
    ti = ti.loc[keep]
    kmols = [mols[i] for i in keep]
    size = ft.size_features(kmols, idx)
    comp = ft.composition_features(kmols, idx)
    rd = ft.rdkit2d_features(kmols, None, idx)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    E = pd.DataFrame(np.vstack([gen.GetCountFingerprintAsNumPy(m).astype(float) for m in kmols]),
                     index=idx, columns=[f"ecfp.{i}" for i in range(2048)])
    tia = np.arcsinh(ti).add_prefix("ti.")
    rda = np.arcsinh(rd).add_prefix("rd.")
    sets = {"size": size, "composition": comp, "TIa": tia,
            "TIa_composition": pd.concat([tia, comp], axis=1),
            "rdkit2d_asinh": rda, "rdkit2d_asinh_TIa": pd.concat([rda, tia], axis=1),
            "ECFP": E, "ECFP_TIa": pd.concat([E, tia], axis=1)}
    meta = pd.DataFrame({"formula": form.loc[keep], "constitution": cons.loc[keep]}, index=idx)
    info["registry"] = list(ti.columns)
    info["registry_sha256"] = hashlib.sha256("|".join(ti.columns).encode()).hexdigest()
    info["feature_sha256"] = {k: hsh(v) for k, v in sets.items()}
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--negative-control", action="store_true")
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

    for endpoint, estimand in CELLS:
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
