"""Independent recomputation of E3 from archived predictions + 06b R6 failure analysis.

No topoqspr imports: fold metrics, paired corrected-t intervals, BH over the secondary family,
and the external-test statistics (point estimates exactly; intervals with an independent
formula-cluster bootstrap using a different seed, so agreement is a robustness check, not an
identity) are all recomputed from predictions.

    python -I docs/research/scripts/e3_recompute_and_failures.py ~/topoqspr_runs/research_2026-10-13_E3
"""
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from rdkit import Chem
from scipy import stats as st


def r2(y, p):
    s = ((y - y.mean()) ** 2).sum()
    return np.nan if s == 0 else 1 - ((y - p) ** 2).sum() / s


def corrected_t(d, n1, n2, alpha=0.05):
    d = np.asarray(d, float); J = d.size; s2 = d.var(ddof=1)
    if s2 == 0:
        return d.mean(), np.nan, np.nan, np.nan
    se = np.sqrt((1 / J + np.mean(n2) / np.mean(n1)) * s2); c = st.t.ppf(1 - alpha / 2, J - 1)
    return d.mean(), d.mean() - c * se, d.mean() + c * se, 2 * st.t.sf(abs(d.mean() / se), J - 1)


def bh(p):
    p = np.asarray(p, float); m = p.size; o = np.argsort(p); q = np.empty(m); run = 1.0
    for j in range(m - 1, -1, -1):
        run = min(run, p[o[j]] * m / (j + 1)); q[o[j]] = run
    return q


def main(run: Path) -> None:
    out = Path(__file__).resolve().parents[1] / "out" / "E3"
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for cdir in sorted((run / "cells").glob("*_E3*")):
        P = pd.read_csv(cdir / "predictions.csv.gz", float_precision="round_trip")
        F = pd.read_csv(cdir / "folds.csv", float_precision="round_trip")
        arch = pd.read_csv(cdir / "paired.csv", float_precision="round_trip")
        fm = P.groupby(["set", "repeat", "fold"]).apply(
            lambda d: pd.Series(dict(r2=r2(d.y_true.to_numpy(), d.y_pred.to_numpy()),
                                     mae=float(np.abs(d.y_true - d.y_pred).mean()), n_test=len(d))),
            include_groups=False).reset_index()
        fm = fm.merge(F[["set", "repeat", "fold", "n_train", "r2"]], on=["set", "repeat", "fold"],
                      suffixes=("", "_arch"))
        max_fold_diff = float(np.nanmax(np.abs(fm.r2 - fm.r2_arch)))
        for _, a in arch.iterrows():
            A = fm[fm.set == a.set].sort_values(["repeat", "fold"])
            B = fm[fm.set == a.reference].sort_values(["repeat", "fold"])
            mean, lo, hi, p = corrected_t(A.r2.to_numpy() - B.r2.to_numpy(), A.n_train, A.n_test)
            dm, dlo, dhi, _ = corrected_t(B.mae.to_numpy() - A.mae.to_numpy(), A.n_train, A.n_test)
            rows.append(dict(cell=cdir.name, set=a.set, reference=a.reference, d_r2=mean, lo=lo, hi=hi,
                             p=p, d_mae=dm, mae_lo=dlo, mae_hi=dhi,
                             max_abs_diff_vs_archived=float(np.nanmax(np.abs(
                                 [mean - a.d_r2_mean, lo - a.d_r2_ci_low, hi - a.d_r2_ci_high, p - a.d_r2_p]))),
                             max_fold_r2_diff=max_fold_diff))
    R = pd.DataFrame(rows)
    if len(R):
        primary = (R.cell == "mp_bradley_E3a") & (R.set == "TIa_composition") & (R.reference == "composition")
        R["role"] = np.where(primary, "primary", "secondary")
        sec = R.role == "secondary"
        R.loc[sec, "q_bh_internal_only"] = bh(R.loc[sec, "p"].to_numpy())
    # external (06b R3): point estimates recomputed exactly; CIs by an independent bootstrap
    ext = run / "external"
    erows = []
    if (ext / "predictions.csv.gz").exists():
        sys.path.insert(0, str(Path(__file__).parent))
        import e3_isomer as E3  # data-building only (no fitting): rebuilds the training formulas
        _, tmeta_train, _, _ = E3.build("mp_bradley", restrict=E3.physprop_free_bradley_keys())
        pr = pd.read_csv(ext / "predictions.csv.gz", index_col=0, float_precision="round_trip")
        known = pr.formula.isin(set(tmeta_train.formula))
        arch_e = pd.read_csv(ext / "paired.csv", float_precision="round_trip")
        rng = np.random.default_rng(99)  # deliberately not the run's seed (2026)
        for test, mask in (("known_formula", known), ("unseen_formula", ~known)):
            d = pr[mask]; yy = d.y.to_numpy(); Fm = d.formula.to_numpy()
            uf = np.unique(Fm); mem = [np.flatnonzero(Fm == f) for f in uf]
            boots = [np.concatenate([mem[j] for j in rng.integers(0, uf.size, uf.size)]) for _ in range(2000)]
            for _, a_ in arch_e[arch_e.test == test].iterrows():
                pa, pb = d[a_.set].to_numpy(), d[a_.reference].to_numpy()
                dr = r2(yy, pa) - r2(yy, pb)
                bs = np.array([r2(yy[i], pa[i]) - r2(yy[i], pb[i]) for i in boots])
                erows.append(dict(test=test, n=len(d), set=a_.set, reference=a_.reference, d_r2=dr,
                                  lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                                  p_boot=float(min(1, 2 * min((bs <= 0).mean(), (bs >= 0).mean()))),
                                  point_diff_vs_archived=abs(dr - a_.d_r2),
                                  ci_shift_vs_archived=max(abs(np.percentile(bs, 2.5) - a_.d_r2_ci_low),
                                                           abs(np.percentile(bs, 97.5) - a_.d_r2_ci_high))))
    E = pd.DataFrame(erows)
    if len(E) and len(R):
        allp = np.r_[R.loc[R.role == "secondary", "p"].to_numpy(), E.p_boot.to_numpy()]
        q = bh(allp)
        R.loc[R.role == "secondary", "q_bh_family"] = q[: (R.role == "secondary").sum()]
        E["q_bh_family"] = q[(R.role == "secondary").sum():]
        E.to_csv(out / "recomputed_external.csv", index=False)
    # failure analysis on MP E3a (06b R6)
    fa = {}
    c = run / "cells" / "mp_bradley_E3a"
    if (c / "predictions.csv.gz").exists():
        P = pd.read_csv(c / "predictions.csv.gz", float_precision="round_trip")
        meta = pd.read_csv(run / "inclusion_mp_bradley.csv", index_col=0)
        sys.path.insert(0, str(Path(__file__).parent))
        import e3_isomer as E3
        recs, _ = E3.parents("mp_bradley")
        smi = recs.drop_duplicates("inchikey").set_index("inchikey").canonical_smiles
        err = P.assign(ae=(P.y_true - P.y_pred).abs()).groupby(["set", "row"]).ae.mean().unstack(0)
        fa["n_molecules"] = len(err)
        fa["mae_by_set"] = err.mean().round(3).to_dict()
        fa["n_abs_error_gt_100C_TIa_composition"] = int((err["TIa_composition"] > 100).sum())
        g = err.join(meta).groupby("formula")["TIa_composition"].mean().sort_values(ascending=False)
        fa["worst_10_formulas_TIa_composition_MAE"] = g.head(10).round(1).to_dict()
        if smi is not None:
            def graph(s):
                m = Chem.MolFromSmiles(s); G = nx.Graph(); G.add_nodes_from(range(m.GetNumAtoms()))
                G.add_edges_from((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()); return G
            iso_ids = set()
            for f, grp in meta.groupby("formula"):
                ids = list(grp.index); Gs = {i: graph(smi[i]) for i in ids}
                for x in range(len(ids)):
                    for y in range(x + 1, len(ids)):
                        if ids[x][:14] != ids[y][:14] and nx.is_isomorphic(Gs[ids[x]], Gs[ids[y]]):
                            iso_ids.update([ids[x], ids[y]])
            m_iso = err.index.isin(iso_ids)
            fa["n_molecules_in_isomorphic_pairs"] = int(m_iso.sum())
            d = err["TIa_composition"] - err["ECFP"]
            fa["mean_AE_TIa_comp_minus_ECFP_isomorphic"] = round(float(d[m_iso].mean()), 3)
            fa["mean_AE_TIa_comp_minus_ECFP_rest"] = round(float(d[~m_iso].mean()), 3)
            d2 = err["composition"] - err["TIa_composition"]
            fa["mean_AE_gain_TIa_over_comp_isomorphic"] = round(float(d2[m_iso].mean()), 3)
            fa["mean_AE_gain_TIa_over_comp_rest"] = round(float(d2[~m_iso].mean()), 3)
    R.to_csv(out / "recomputed_paired.csv", index=False)
    pd.Series(fa, dtype=object).to_json(out / "failure_analysis.json", indent=1)
    with pd.option_context("display.width", 220):
        print(R.round(4).to_string(index=False))
        print(E.round(4).to_string(index=False) if len(E) else "external: not run")
        print(fa)


if __name__ == "__main__":
    main(Path(sys.argv[1]).expanduser())
