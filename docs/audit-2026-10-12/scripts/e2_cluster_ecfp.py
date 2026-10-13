"""AMD-2026-10-12-E2 (exploratory): see docs/audit-2026-10-12/EXPERIMENT_AMENDMENT.md.

    # smoke test (permuted target, logd74_az, cluster split, seed 0; numbers discarded)
    python -I docs/audit-2026-10-12/scripts/e2_cluster_ecfp.py --smoke \
        ~/topoqspr_runs/main_confirmatory /tmp/e2_smoke
    # the run
    python -I docs/audit-2026-10-12/scripts/e2_cluster_ecfp.py \
        ~/topoqspr_runs/main_confirmatory ~/topoqspr_runs/audit_2026-10-12_E2
"""
import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
import sklearn
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

from topoqspr import evaluate as ev
from topoqspr import stats

FAMILIES = ["ridge", "elasticnet", "pls", "knn", "svr", "extratrees"]
TARGETS = ["mp_bradley", "logp_opera", "logd74_az", "logs_aqsoldb"]
SPLITS = ["cluster", "random", "scaffold"]
E1 = {"size": "size", "composition": "composition", "TIa_composition": "TI_asinh_composition",
      "TI_composition": "TI_composition", "rdkit2d_asinh": "rdkit2d_asinh"}
CONTRASTS = [("TIa_composition", "composition"), ("ECFP", "composition"),
             ("ECFP", "TIa_composition"), ("rdkit2d_asinh", "TIa_composition"),
             ("TIa_composition", "size")]
MARGIN = 0.05


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ecfp(smiles: pd.Series) -> tuple[pd.DataFrame, list]:
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    counts, bits = [], []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        counts.append(gen.GetCountFingerprintAsNumPy(m).astype(float))
        bits.append(gen.GetFingerprint(m))
    X = pd.DataFrame(np.vstack(counts), index=smiles.index,
                     columns=[f"ecfp.{i}" for i in range(2048)])
    return X, bits


def butina(bits: list, ids: pd.Index, threshold: float = 0.4) -> pd.Series:
    dists = []
    for i in range(1, len(bits)):
        sims = DataStructs.BulkTanimotoSimilarity(bits[i], bits[:i])
        dists.extend(1.0 - x for x in sims)
    clusters = Butina.ClusterData(dists, len(bits), threshold, isDistData=True)
    lab = np.empty(len(bits), dtype=object)
    for c, members in enumerate(clusters):
        for m in members:
            lab[m] = f"c{c}"
    return pd.Series(lab, index=ids)


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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--n-jobs", type=int, default=14)
    a = ap.parse_args()
    run, out = Path(a.run).expanduser(), Path(a.out).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    F = pd.read_csv(run / "features.csv", index_col="id", float_precision="round_trip")
    Y = pd.read_csv(run / "targets.csv", index_col="id", float_precision="round_trip")
    M = pd.read_csv(run / "molecules.csv", index_col="id")
    ok = M[M.status == "ok"]
    conn = ok.connectivity_group.fillna(pd.Series(ok.index, index=ok.index))

    E, bits = ecfp(ok.smiles)
    clus = butina(bits, ok.index)
    first_cluster = clus.loc[conn.values].set_axis(ok.index)  # members -> first member's cluster
    scaf = pd.Series(ev.scaffold_groups([Chem.MolFromSmiles(s) for s in ok.smiles]), index=ok.index)
    first_scaf = scaf.loc[conn.values].set_axis(ok.index)
    labels = {"random": conn, "scaffold": first_scaf, "cluster": first_cluster}
    pd.DataFrame(dict(connectivity=conn, scaffold=first_scaf, cluster=first_cluster)) \
        .to_csv(out / "cluster_labels.csv")
    (out / "ecfp_sha256.txt").write_text(hashlib.sha256(E.to_numpy().tobytes()).hexdigest() + "\n")

    def block(p):
        return F[[c for c in F.columns if c.startswith(p)]]
    comp, ti, rd = block("composition."), block("ti."), block("rdkit2d.")
    size = F[["size.n_heavy", "size.n_bonds"]]
    sets_all = {"size": size, "composition": comp,
                "TIa_composition": pd.concat([np.arcsinh(ti), comp], axis=1),
                "TI_composition": pd.concat([ti, comp], axis=1),
                "ECFP": E, "rdkit2d_asinh": np.arcsinh(rd)}

    e1_path = Path("~/topoqspr_runs/audit_2026-10-12_sensitivity/summary.csv").expanduser()
    e1 = pd.read_csv(e1_path) if e1_path.exists() and not a.smoke else None
    targets, splits, seeds = TARGETS, SPLITS, (0, 1, 2, 3, 4)
    if a.smoke:
        targets, splits, seeds = ["logd74_az"], ["cluster"], (0,)

    summ, pairs, folds, preds, gstats = [], [], [], [], []
    for split in splits:
        for t in targets:
            ids = Y.index[Y[t].notna() & (M.loc[Y.index, "status"] == "ok")]
            y = Y.loc[ids, t].astype(float)
            if a.smoke:
                y = pd.Series(np.random.default_rng(12345).permutation(y.to_numpy()), index=ids)
            g = labels[split].loc[ids].to_numpy(dtype=object)
            vc = pd.Series(g).value_counts()
            gstats.append(dict(split=split, target=t, n=len(ids), n_groups=int(vc.size),
                               largest_group=int(vc.iloc[0]),
                               largest_share=float(vc.iloc[0] / len(ids)),
                               n_singletons=int((vc == 1).sum())))
            t0 = time.time()
            res = {k: ev.nested_cv(X.loc[ids], y, FAMILIES, groups=g, seeds=seeds,
                                   n_jobs=a.n_jobs) for k, X in sets_all.items()}
            print(f"{split} {t} n={len(y)} groups={vc.size} {time.time() - t0:.0f}s", flush=True)
            for k, r in res.items():
                s = r.summary()
                summ.append(dict(split=split, target=t, set=k, n=len(y), p=r.meta["p"],
                                 r2_fold_mean=s["r2_fold_mean"], r2_pooled_mean=s["r2_pooled_mean"],
                                 mae_fold_mean=s["mae_fold_mean"], min_fold_r2=float(r.folds.r2.min()),
                                 n_failed_folds=s["n_failed_folds"],
                                 family_counts=json.dumps(s["family_counts"]), reason=s["reason"]))
                folds.append(r.folds.assign(split=split, target=t, set=k))
                preds.append(r.predictions.assign(split=split, target=t, set=k))
            for x, ref in CONTRASTS:
                row = dict(split=split, target=t, **paired(res, x, ref))
                if (x, ref) == ("TIa_composition", "composition"):
                    lo, hi = row["d_r2_ci90_low"], row["d_r2_ci90_high"]
                    row["equivalent_within_margin"] = bool(lo >= -MARGIN and hi <= MARGIN)
                pairs.append(row)
            if e1 is not None and split in ("random", "scaffold"):
                for k, k1 in E1.items():
                    ref = e1[(e1.split == split) & (e1.target == t) & (e1.set == k1)].r2_fold_mean.item()
                    mine = [r for r in summ if r["split"] == split and r["target"] == t
                            and r["set"] == k][0]["r2_fold_mean"]
                    assert np.isclose(ref, mine, rtol=1e-9, atol=1e-9), (split, t, k, ref, mine)
            pd.DataFrame(summ).to_csv(out / "summary.csv", index=False)
            pd.DataFrame(pairs).to_csv(out / "paired_tests.csv", index=False)
            pd.DataFrame(gstats).to_csv(out / "group_stats.csv", index=False)

    P = pd.DataFrame(pairs)
    if not a.smoke:
        h1 = P[(P.split == "cluster") & (P.set == "TIa_composition") & (P.reference == "composition")]
        q, _ = stats.benjamini_hochberg(h1.d_r2_p.to_numpy(), q=0.05)
        P.loc[h1.index, "h1_q_bh"] = q
        P.to_csv(out / "paired_tests.csv", index=False)
    pd.concat(folds).to_csv(out / "folds.csv", index=False)
    pd.concat(preds).to_csv(out / "predictions.csv.gz", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                            cwd=Path(__file__).parent).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"],
                                capture_output=True, text=True, cwd=Path(__file__).parent).stdout.strip())
    info = dict(amendment="AMD-2026-10-12-E2", smoke=a.smoke, commit=commit, dirty=dirty,
                inputs={f: sha(run / f) for f in ("features.csv", "molecules.csv", "targets.csv")},
                python=platform.python_version(), rdkit=rdkit.__version__, sklearn=sklearn.__version__,
                n_clusters_all=int(clus.nunique()), wall_s=round(time.time() - t_start, 1),
                e1_reproduction="asserted" if e1 is not None else "not applicable")
    (out / "run_info.json").write_text(json.dumps(info, indent=1) + "\n")
    print("done", json.dumps(info))


if __name__ == "__main__":
    main()
