"""Notebook 02 / doc 14 (revised): four-arm fixed-learner comparison on melting point.

Arms (built from the E3 feature code, e3_isomer.featurise, unchanged):
    C       = composition                    (element counts)
    C+E     = edgetypes_composition          (9 degree-pair edge counts m_ij + max degree + C)
    C+T     = TIa_composition                (asinh of the 104 registry indices + C)
    C+E+T   = edge-type block + asinh(TI) block + C   (new; the union, column order E, T, C)

Learners: each of the six E3 families *alone* (families=[F]); hyperparameters tuned by inner
grouped 3-fold CV inside the training data with the E3 grids. All preprocessing (median
imputation, constant/duplicate-column removal, scaling, target scaling) is inside the sklearn
pipeline, so it is fitted on training rows only (topoqspr.models).

Contexts:
    mp_E3a    Bradley isomer set, groups = constitution, seeds 0-2 x 5 outer folds (E3 splits)
    mp_E3b    same, groups = formula
    external  train on the PHYSPROP-free Bradley isomer set, test on OPERA/PHYSPROP MP parents
              (known/unseen formula, sizes asserted 1,014 / 7,327). Training seeds 0-4 for
              extra trees (seed changes the forest and the inner splits); the other families
              are deterministic given the inner splits and are fitted with seed 0 only.

Work is split into units (context, family, arm, seed, fold); each unit is written atomically
with a sha256-verified record and skipped on resume only if its record verifies. Units can be
sharded (--shard i/N) across Colab sessions or machines; results are identical whichever
shard computes a unit.

SMOKE mode permutes the MP labels before fitting (seeded), so a smoke run exercises the full
code path but cannot reveal the answer to the question. Its numbers are meaningless. BENCH mode
is the same with all six families on one internal fold (mp_E3b, seed 0, fold 0) and the external
seed-0 fit: it times every (family, arm) unit for the compute budget, again on permuted labels.

    python -I docs/research/scripts/fixed_learner.py build  OUT
    python -I docs/research/scripts/fixed_learner.py run    OUT --mode smoke
    python -I docs/research/scripts/fixed_learner.py run    OUT --mode full --shard 0/4 --n-jobs 2
    python -I docs/research/scripts/fixed_learner.py status OUT --mode full
    python -I docs/research/scripts/fixed_learner.py analyse OUT --mode full
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import repro_common as RC  # noqa: E402

FAMILIES = ["ridge", "elasticnet", "pls", "knn", "svr", "extratrees"]
CO_PRIMARY_LEARNERS = ["extratrees", "svr"]
ARMS = ["C", "C+E", "C+T", "C+E+T"]
#: (set, reference, role). Co-primary estimands per the revised doc 14.
CONTRASTS = [("C+T", "C", "co-primary: topology beyond composition"),
             ("C+E+T", "C+E", "co-primary: registry beyond edge types (nested)"),
             ("C+T", "C+E", "secondary: non-nested compact-vs-full comparison"),
             ("C+E+T", "C+T", "secondary: redundancy diagnostic"),
             ("C+E", "C", "secondary: edge types beyond composition")]
CONTEXTS = ["external", "mp_E3a", "mp_E3b"]
SEEDS_INTERNAL = (0, 1, 2)
TRAIN_SEEDS_EXTERNAL = {"extratrees": (0, 1, 2, 3, 4)}  # others: (0,)
DELTA = 0.02  # candidate smallest effect of interest (R2); see doc 14 §6 for its status
SMOKE = dict(contexts={"mp_E3b": dict(seeds=(0,), folds=(0,)), "external": {}},
             families=["ridge", "extratrees"], permute_seed=20261017)
PROTOCOL_VERSION = "doc14-r2"


# ----------------------------------------------------------------------------- data
def arm_frames(sets: dict) -> dict:
    E = sets["edgetypes_composition"]; C = sets["composition"]; T = sets["TIa"]
    ecols = [c for c in E.columns if c not in C.columns]
    assert list(sets["TIa_composition"].columns) == list(T.columns) + list(C.columns)
    return {"C": C, "C+E": E, "C+T": sets["TIa_composition"],
            "C+E+T": pd.concat([E[ecols], T, C], axis=1)}


def feature_manifest(arms: dict, raw_ti: pd.DataFrame | None = None) -> dict:
    """Names, dimensions, numeric rank, constant and duplicate columns per arm."""
    out = {}
    for k, X in arms.items():
        A = X.to_numpy(float)
        const = [c for c in X.columns if X[c].nunique(dropna=False) <= 1]
        dup = int(X.T.duplicated().sum())
        Z = A - np.nanmean(A, 0)
        Z = Z / np.where(np.nanstd(Z, 0) > 0, np.nanstd(Z, 0), 1)
        out[k] = dict(p=int(A.shape[1]), n=int(A.shape[0]), n_missing=int(np.isnan(A).sum()),
                      rank_standardised=int(np.linalg.matrix_rank(np.nan_to_num(Z), tol=1e-8 * max(A.shape))),
                      constant_columns=const, n_duplicate_columns=dup,
                      columns=list(X.columns), sha256=RC_hash(X))
    return out


def RC_hash(df: pd.DataFrame) -> str:  # identical to e3_isomer.hsh
    return hashlib.sha256(np.ascontiguousarray(df.to_numpy(dtype=float)).tobytes()
                          + "|".join(map(str, df.columns)).encode()).hexdigest()


def build(out: Path) -> dict:
    """Rebuild every context's arms from public sources with the E3 code; fail closed on any
    mismatch with the archived E3 inclusion lists, labels, splits or feature hashes."""
    import e3_isomer as E3
    from topoqspr import expdata
    ck = RC.Checks("02_build")
    ck.info("sources", json.dumps(RC.prefetch_sources()))
    ddir = out / "data"; ddir.mkdir(parents=True, exist_ok=True)
    ref_int = json.loads((RC.E3_DATA / "feature_sha256_mp_internal.json").read_text())
    ext_info = json.loads((RC.E3_DATA / "external/info.json").read_text())
    manifest = {}
    # internal (mp_E3a / mp_E3b share the data)
    sets, meta, y, info = E3.build("mp_bradley")
    arch_meta = pd.read_csv(RC.E3_DATA / "inclusion_mp_bradley.csv", index_col=0)
    RC.check_ids(ck, "MP internal inclusion", meta.index, arch_meta.index)
    ck.require("MP internal: inclusion order and formulas identical to archive",
               meta.index.equals(arch_meta.index) and meta.formula.equals(arch_meta.formula), "", "")
    Pa = pd.read_csv(RC.E3_DATA / "cells/mp_bradley_E3a/predictions.csv.gz", float_precision="round_trip")
    ya = Pa.drop_duplicates("row").set_index("row").y_true
    ck.close("MP internal: labels equal archived y_true (max |diff|, degC)",
             float(np.abs(y.loc[ya.index] - ya).max()), 0.0, 1e-9)
    for k in ("composition", "edgetypes_composition", "TIa_composition", "TIa"):
        ck.require(f"MP internal feature hash {k}", info["feature_sha256"][k] == ref_int["feature_sha256"][k],
                   info["feature_sha256"][k][:16], ref_int["feature_sha256"][k][:16])
    arms = arm_frames(sets)
    manifest["internal"] = feature_manifest(arms)
    np.savez_compressed(ddir / "internal.npz", **{a: arms[a].to_numpy(float) for a in ARMS},
                        y=y.to_numpy(float), ids=np.array(meta.index, dtype=str), formula=meta.formula.to_numpy(str),
                        constitution=meta.constitution.to_numpy(str))
    # external
    tr_sets, tr_meta, tr_y, tr_info = E3.build("mp_bradley", restrict=E3.physprop_free_bradley_keys())
    for k in ("composition", "edgetypes_composition", "TIa_composition", "TIa"):
        ck.require(f"external train feature hash {k} = archived external/info.json",
                   tr_info["feature_sha256"][k] == ext_info["train"]["feature_sha256"][k],
                   tr_info["feature_sha256"][k][:16], ext_info["train"]["feature_sha256"][k][:16])
    recs, _ = E3.parents("mp_opera")
    rep = recs.drop_duplicates("inchikey")
    rep = rep[~rep.inchikey.str[:14].isin(set(tr_meta.constitution))]
    te_sets, te_meta, _ = E3.featurise(rep.inchikey, rep.canonical_smiles)
    split = pd.read_csv(RC.E3_DATA / "external/test_split.csv").set_index("id")
    RC.check_ids(ck, "external test", te_meta.index, split.index)
    ck.require("external test: order identical to archived split", list(te_meta.index) == list(split.index), "", "")
    known = te_meta.formula.isin(set(tr_meta.formula))
    ck.require("external: known/unseen split identical to archived", bool(
        (np.where(known, "known_formula", "unseen_formula") == split.test.to_numpy()).all()), "", "")
    ck.require("external: sizes 1,014 / 7,327", (int(known.sum()), int((~known).sum())) == (1014, 7327),
               (int(known.sum()), int((~known).sum())), (1014, 7327))
    RC.check_group_disjoint(ck, "external constitutions", tr_meta.constitution, te_meta.index.str[:14])
    ck.require("external: no test formula marked unseen occurs in training",
               not (set(te_meta.formula[~known]) & set(tr_meta.formula)), "", "")
    tr_arms, te_arms = arm_frames(tr_sets), arm_frames(te_sets)
    for a in ARMS:
        ck.require(f"external {a}: test columns = train columns", list(te_arms[a].columns) == list(tr_arms[a].columns), "", "")
    manifest["external_train"] = feature_manifest(tr_arms)
    np.savez_compressed(ddir / "external.npz",
                        **{f"train_{a}": tr_arms[a].to_numpy(float) for a in ARMS},
                        **{f"test_{a}": te_arms[a].to_numpy(float) for a in ARMS},
                        y_train=tr_y.to_numpy(float), train_ids=np.array(tr_meta.index, dtype=str),
                        train_constitution=tr_meta.constitution.to_numpy(str),
                        train_formula=tr_meta.formula.to_numpy(str),
                        test_ids=np.array(te_meta.index, dtype=str), test_formula=te_meta.formula.to_numpy(str),
                        test_split=split.test.to_numpy(str))
    # the external *test labels* are deliberately not stored here: units never see them;
    # analyse() reads them from the restored public dataset.
    RC.atomic_write_text(json.dumps(manifest, indent=1), out / "feature_manifest.json")
    (ddir / "SHA256").write_text("".join(f"{RC.sha256_file(ddir / f)}  {f}\n" for f in ("internal.npz", "external.npz")))
    ck.finish(out)
    return manifest


def load(out: Path):
    ddir = out / "data"
    for line in (ddir / "SHA256").read_text().splitlines():
        h, f = line.split("  ")
        if RC.sha256_file(ddir / f) != h:
            raise SystemExit(f"{ddir / f} does not match its recorded sha256; run `build` again")
    with np.load(ddir / "internal.npz", allow_pickle=False) as a, np.load(ddir / "external.npz", allow_pickle=False) as b:
        return dict(a), dict(b)  # plain dicts: picklable for joblib workers


# ----------------------------------------------------------------------------- units
def units(mode: str) -> list[dict]:
    fams = SMOKE["families"] if mode == "smoke" else FAMILIES
    out = []
    for ctx in CONTEXTS:
        if mode in ("smoke", "bench") and ctx not in SMOKE["contexts"]:
            continue
        for F in fams:
            for a in ARMS:
                if ctx == "external":
                    for s in TRAIN_SEEDS_EXTERNAL.get(F, (0,)):
                        if mode in ("smoke", "bench") and s > 0:
                            continue
                        out.append(dict(context=ctx, family=F, arm=a, seed=s, fold=-1))
                else:
                    seeds = SMOKE["contexts"][ctx]["seeds"] if mode in ("smoke", "bench") else SEEDS_INTERNAL
                    folds = SMOKE["contexts"][ctx]["folds"] if mode in ("smoke", "bench") else range(5)
                    for s in seeds:
                        for k in folds:
                            out.append(dict(context=ctx, family=F, arm=a, seed=s, fold=k))
    return out


def unit_id(u: dict) -> str:
    return f"{u['context']}/{u['family']}/{u['arm'].replace('+', 'p')}/s{u['seed']}_f{u['fold']}"


def config_hash(mode: str, data_sha: str) -> str:
    cfg = dict(PROTOCOL_VERSION=PROTOCOL_VERSION, mode=mode, FAMILIES=FAMILIES, ARMS=ARMS,
               SEEDS_INTERNAL=SEEDS_INTERNAL, TRAIN_SEEDS_EXTERNAL=TRAIN_SEEDS_EXTERNAL,
               SMOKE=SMOKE if mode != "full" else None, data=data_sha)
    return hashlib.sha256(json.dumps(cfg, sort_keys=True, default=str).encode()).hexdigest()[:16]


def _done(rec_path: Path, cfg: str) -> bool:
    if not rec_path.exists():
        return False
    try:
        rec = json.loads(rec_path.read_text())
        p = rec_path.with_name(rec_path.name.replace(".json", ".csv.gz"))
        return rec.get("config") == cfg and p.exists() and RC.sha256_file(p) == rec["predictions_sha256"]
    except Exception:  # noqa: BLE001 - an unreadable record is simply recomputed
        return False


def run_unit(u: dict, I, X, mode: str, cfg: str, root: Path) -> dict:
    from topoqspr import evaluate as ev
    t0 = time.time()
    rng_y = (lambda y: np.random.default_rng(SMOKE["permute_seed"]).permutation(y)) if mode != "full" else (lambda y: y)
    if u["context"] == "external":
        Xtr, Xte = X[f"train_{u['arm']}"], X[f"test_{u['arm']}"]
        ytr = rng_y(X["y_train"])
        sel = ev.select_and_refit(Xtr, ytr, [u["family"]], inner=3, random_state=u["seed"],
                                  groups=X["train_constitution"].astype(object))
        yp = sel.estimator.predict(Xte) if sel.estimator is not None else np.full(len(Xte), np.nan)
        pred = pd.DataFrame(dict(id=X["test_ids"], y_pred=np.asarray(yp, float).ravel()))
        rec = dict(n_train=int(len(ytr)), n_test=int(len(Xte)), params=sel.params, inner_mse=sel.inner_mse,
                   reason=sel.reason, inner_table=sel.table.to_dict("records"),
                   train_groups_sha256=hashlib.sha256("|".join(sorted(set(X["train_constitution"]))).encode()).hexdigest())
    else:
        A = I[u["arm"]]; y = rng_y(I["y"])
        groups = (I["constitution"] if u["context"] == "mp_E3a" else I["formula"]).astype(object)
        splits = ev.make_outer_splits(len(y), 5, SEEDS_INTERNAL, groups)
        split = next(s for s in splits if s[2] == u["seed"] and s[1] == u["fold"])
        tr, te = split[3], split[4]
        inter = set(groups[tr]) & set(groups[te])
        if inter:
            raise SystemExit(f"LEAKAGE: {len(inter)} groups in both train and test for {unit_id(u)}")
        fold_row, p, tab, _ = ev._run_fold(A, y, split, [u["family"]], 3, groups, 0, None, False)
        pred = pd.DataFrame(dict(id=I["ids"][p.pos.to_numpy()], y_true=p.y_true.to_numpy(), y_pred=p.y_pred.to_numpy()))
        rec = dict(fold_row, inner_table=tab.to_dict("records"),
                   train_groups_sha256=hashlib.sha256("|".join(sorted(set(groups[tr]))).encode()).hexdigest(),
                   test_groups_sha256=hashlib.sha256("|".join(sorted(set(groups[te]))).encode()).hexdigest(),
                   n_train_groups=len(set(groups[tr])), n_test_groups=len(set(groups[te])))
        rec["params"] = json.loads(rec["params"]) if isinstance(rec.get("params"), str) and rec["params"] else rec.get("params")
    base = root / unit_id(u)
    base.parent.mkdir(parents=True, exist_ok=True)
    buf = io.BytesIO(); pred.to_csv(buf, index=False, compression={"method": "gzip", "mtime": 0}); data = buf.getvalue()
    tmp = base.with_name(base.name + ".csv.gz.tmp"); tmp.write_bytes(data); tmp.replace(base.with_name(base.name + ".csv.gz"))
    rec.update(u, config=cfg, mode=mode, seconds=round(time.time() - t0, 2), finite=bool(np.isfinite(pred.y_pred).all()),
               predictions_sha256=hashlib.sha256(data).hexdigest(), git_sha=RC.git_sha())
    RC.atomic_write_text(json.dumps(rec, default=str, indent=1), base.with_name(base.name + ".json"))
    return rec


def run(out: Path, mode: str, shard: str = "0/1", n_jobs: int = 1, max_units: int | None = None) -> pd.DataFrame:
    I, X = load(out)
    data_sha = (out / "data/SHA256").read_text()
    cfg = config_hash(mode, data_sha)
    root = out / f"units_{mode}"
    i, n = map(int, shard.split("/"))
    todo = [u for j, u in enumerate(units(mode)) if j % n == i]
    pending = [u for u in todo if not _done(root / (unit_id(u) + ".json"), cfg)]
    print(f"{mode}: shard {shard}: {len(todo)} units, {len(todo) - len(pending)} verified done, "
          f"{len(pending)} to run (n_jobs={n_jobs})", flush=True)
    if max_units is not None:
        pending = pending[:max_units]
    from joblib import Parallel, delayed
    recs = Parallel(n_jobs=n_jobs, verbose=5)(delayed(run_unit)(u, I, X, mode, cfg, root) for u in pending)
    bad = [r for r in recs if not r["finite"]]
    if bad:
        print(f"WARNING: {len(bad)} unit(s) produced non-finite predictions: {[unit_id(r) for r in bad]}")
    return status(out, mode)


def status(out: Path, mode: str) -> pd.DataFrame:
    data_sha = (out / "data/SHA256").read_text()
    cfg = config_hash(mode, data_sha)
    root = out / f"units_{mode}"
    rows = []
    for u in units(mode):
        p = root / (unit_id(u) + ".json")
        ok = _done(p, cfg)
        sec = json.loads(p.read_text())["seconds"] if ok else np.nan
        rows.append(dict(u, done=ok, seconds=sec))
    S = pd.DataFrame(rows)
    print(S.groupby(["context", "family"]).done.agg(["sum", "size"]).rename(columns={"sum": "done", "size": "units"}).to_string())
    print(f"TOTAL {int(S.done.sum())}/{len(S)} units verified; CPU-seconds so far {np.nansum(S.seconds):.0f}", flush=True)
    return S


def budget(out: Path, cpus: int = 2) -> pd.DataFrame:
    """Project FULL CPU time from BENCH unit timings (one internal fold, external seed 0)."""
    B = status(out, "bench")
    if not B.done.all():
        raise SystemExit("run --mode bench first")
    full = pd.DataFrame(units("full"))
    t = B.set_index(["context", "family", "arm"]).seconds
    full["seconds"] = [t.get(("mp_E3b" if u.context != "external" else "external", u.family, u.arm), np.nan)
                       for u in full.itertuples()]
    S = full.groupby(["context", "family"]).seconds.sum().unstack(0) / 3600
    tot = full.seconds.sum() / 3600
    print(S.round(2).to_string())
    print(f"projected FULL: {tot:.1f} CPU-hours (from {len(B)} bench units, E3b fold as proxy for E3a); "
          f"about {tot / cpus:.1f} wall-hours at {cpus} parallel workers if memory allows")
    RC.atomic_write_csv(full, out / "budget_projection.csv", index=False)
    return S


# ----------------------------------------------------------------------------- analysis
def _r2(y, p):
    s = ((y - y.mean()) ** 2).sum()
    return np.nan if s == 0 else 1 - ((y - p) ** 2).sum() / s


def holm(p) -> np.ndarray:
    """Holm step-down adjusted p-values."""
    p = np.asarray(p, float); o = np.argsort(p); k = p.size; out = np.empty(k); run_ = 0.0
    for j, i in enumerate(o):
        run_ = max(run_, min(1.0, (k - j) * p[i])); out[i] = run_
    return out


def classify(lo95, hi95, lo90, hi90, delta=DELTA) -> str:
    """Decision categories of revised doc 14 §6 (one-sided bounds are the 90% two-sided ends)."""
    if lo90 > delta:
        return "superior_by_more_than_delta"     # one-sided 95% LCB > delta
    if lo95 > 0:
        return "superior_not_shown_beyond_delta"  # 2-sided 95% CI excludes 0, LCB <= delta
    if hi90 < -delta:
        return "inferior_by_more_than_delta"
    if hi95 < 0:
        return "inferior_not_shown_beyond_delta"
    if lo90 > -delta and hi90 < delta:
        return "equivalent_within_delta"         # TOST at alpha .05
    return "inconclusive"


def analyse(out: Path, mode: str, restored: Path | None = None) -> dict:
    """Assemble all verified units into tables. External test labels come from the restored
    public dataset (restore_opera_values.py output) and are read only here."""
    from scipy import stats as st
    import reproduce_e3 as R3
    S = status(out, mode)
    if not S.done.all():
        print(f"NOTE: {int((~S.done).sum())} units missing; analysing completed contexts only")
    root = out / f"units_{mode}"
    res = out / f"results_{mode}"; res.mkdir(parents=True, exist_ok=True)
    tables = {}
    # internal
    fold_rows, ctr_rows = [], []
    for ctx in ("mp_E3a", "mp_E3b"):
        d = S[(S.context == ctx)]
        for F in d.family.unique():
            ok = d[(d.family == F)]
            if not ok.done.all():
                continue
            for a in ARMS:
                for _, u in ok[ok.arm == a].iterrows():
                    rec = json.loads((root / (unit_id(u.to_dict()) + ".json")).read_text())
                    fold_rows.append(dict(context=ctx, family=F, arm=a, seed=u.seed, fold=u.fold, r2=rec["r2"],
                                          mae=rec["mae"], rmse=rec["rmse"], n_train=rec["n_train"],
                                          n_test=rec["n_test"], params=json.dumps(rec["params"]), status=rec["status"]))
    FR = pd.DataFrame(fold_rows)
    if len(FR):
        for (ctx, F), g in FR.groupby(["context", "family"]):
            for a, b, role in CONTRASTS:
                A = g[g.arm == a].sort_values(["seed", "fold"]); B = g[g.arm == b].sort_values(["seed", "fold"])
                dd = A.r2.to_numpy() - B.r2.to_numpy()
                m, lo, hi, p = R3.corrected_t(dd, A.n_train, A.n_test)
                _, lo90, hi90, _ = R3.corrected_t(dd, A.n_train, A.n_test, 0.10)
                # repeat-level sensitivity: mean over folds per seed, then t over 3 seeds (no overlap correction needed
                # across seeds beyond shared data; reported as a descriptive check only)
                rep = pd.Series(dd).groupby(A.seed.to_numpy()).mean()
                se = (hi - m) / st.t.ppf(0.975, len(dd) - 1)
                p_gt_delta = float(st.t.sf((m - DELTA) / se, len(dd) - 1)) if se > 0 else np.nan
                ctr_rows.append(dict(context=ctx, family=F, set=a, reference=b, role=role, d_r2=m, ci95_low=lo, ci95_high=hi,
                                     ci90_low=lo90, ci90_high=hi90, p_two_sided=p,
                                     p_one_sided_gt0=(p / 2 if m > 0 else 1 - p / 2),
                                     p_one_sided_gt_delta=p_gt_delta,
                                     r2_set=A.r2.mean(), r2_reference=B.r2.mean(), mae_set=A.mae.mean(),
                                     mae_reference=B.mae.mean(), rmse_set=A.rmse.mean(), rmse_reference=B.rmse.mean(),
                                     n_folds=len(dd), folds_positive=int((dd > 0).sum()),
                                     seed_level_means=json.dumps(rep.round(5).tolist()),
                                     inference="corrected resampled t (Nadeau-Bengio), 15 seed x fold scores; "
                                               "coverage under grouped CV checked in notebook 04"))
    # external
    ext = S[S.context == "external"]
    erows = []
    if len(ext) and restored is not None:
        pr = pd.read_csv(restored / "external/predictions.csv.gz", index_col=0, float_precision="round_trip")
        split = pd.read_csv(restored / "external/test_split.csv").set_index("id")
        ytest = pr.y
        for F in ext.family.unique():
            if not ext[ext.family == F].done.all():
                continue
            for seed in sorted(ext[ext.family == F].seed.unique()):
                P = {}
                for a in ARMS:
                    u = dict(context="external", family=F, arm=a, seed=seed, fold=-1)
                    p = pd.read_csv(root / (unit_id(u) + ".csv.gz"), float_precision="round_trip").set_index("id").y_pred
                    if not p.index.equals(split.index):
                        raise SystemExit("external predictions are not in the archived split order")
                    P[a] = p.to_numpy()
                rng = np.random.default_rng(2026)
                for test in ("known_formula", "unseen_formula"):
                    m = (split.test == test).to_numpy()
                    yy = ytest.to_numpy()[m]; Fm = split.formula.to_numpy()[m]
                    boots = R3.cluster_boot(Fm, rng)
                    for a, b, role in CONTRASTS:
                        pa, pb = P[a][m], P[b][m]
                        dr = _r2(yy, pa) - _r2(yy, pb)
                        bs = np.array([_r2(yy[i], pa[i]) - _r2(yy[i], pb[i]) for i in boots])
                        erows.append(dict(context=f"external_{test}", family=F, train_seed=seed, set=a, reference=b,
                                          role=role, d_r2=dr, ci95_low=np.percentile(bs, 2.5),
                                          ci95_high=np.percentile(bs, 97.5), ci90_low=np.percentile(bs, 5),
                                          ci90_high=np.percentile(bs, 95), boot_se=float(bs.std(ddof=1)),
                                          p_one_sided_gt0=float((bs <= 0).mean()),
                                          p_one_sided_gt_delta=float((bs <= DELTA).mean()),
                                          p_two_sided=float(min(1, 2 * min((bs <= 0).mean(), (bs >= 0).mean()))),
                                          r2_set=_r2(yy, pa), r2_reference=_r2(yy, pb),
                                          mae_set=float(np.abs(yy - pa).mean()), mae_reference=float(np.abs(yy - pb).mean()),
                                          rmse_set=float(np.sqrt(((yy - pa) ** 2).mean())),
                                          rmse_reference=float(np.sqrt(((yy - pb) ** 2).mean())),
                                          n=int(m.sum()), n_formulas=int(np.unique(Fm).size),
                                          inference="paired formula-cluster bootstrap (2,000, seed 2026); training held "
                                                    "fixed: covers test-set sampling only"))
    C = pd.concat([pd.DataFrame(ctr_rows), pd.DataFrame(erows)], ignore_index=True)
    if len(C):
        C["category"] = [classify(r.ci95_low, r.ci95_high, r.ci90_low, r.ci90_high) for r in C.itertuples()]
        # multiplicity: Holm over the co-primary family (2 estimands x 2 learners x 3 contexts, seed 0);
        # BH over everything else; both one-sided H0: delta <= 0
        main = C[(C.role.str.startswith("co-primary")) & C.family.isin(CO_PRIMARY_LEARNERS)
                 & C.context.isin(["mp_E3a", "mp_E3b", "external_unseen_formula"])
                 & (C.get("train_seed", pd.Series(0, index=C.index)).fillna(0) == 0)]
        if len(main):
            for col, name in (("p_one_sided_gt0", "p_holm_gt0"), ("p_one_sided_gt_delta", "p_holm_gt_delta")):
                C.loc[main.index, name] = holm(main[col].to_numpy())
            C.loc[main.index, "family_size_coprimary"] = len(main)
        rest = C.index.difference(main.index)
        C.loc[rest, "q_bh_secondary"] = RC.bh(C.loc[rest, "p_one_sided_gt0"].to_numpy())
        if mode != "full":
            C.insert(0, "WARNING", f"{mode.upper()}: permuted training labels, numbers meaningless")
        RC.atomic_write_csv(C, res / "contrasts.csv", index=False)
    if len(FR):
        RC.atomic_write_csv(FR, res / "folds.csv", index=False)
    # training-seed variability (extra trees, external)
    if len(erows):
        E = pd.DataFrame(erows)
        sv = E[E.family == "extratrees"].groupby(["context", "set", "reference"]).agg(
            n_seeds=("d_r2", "size"), d_r2_mean=("d_r2", "mean"), d_r2_seed_sd=("d_r2", "std"),
            boot_se_mean=("boot_se", "mean")).reset_index()
        RC.atomic_write_csv(sv, res / "external_training_seed_variability.csv", index=False)
        tables["seed_var"] = sv
    tables["contrasts"] = C
    return tables


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "run", "status", "analyse", "budget"])
    ap.add_argument("out")
    ap.add_argument("--mode", choices=["smoke", "bench", "full"], default="smoke")
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--max-units", type=int)
    ap.add_argument("--restored", help="restored public E3 dir (external test labels), for analyse")
    a = ap.parse_args()
    out = Path(a.out).expanduser()
    if a.cmd == "build":
        build(out)
    elif a.cmd == "run":
        if a.mode == "full" and not (out / "APPROVED_FULL_RUN").exists():
            raise SystemExit("full mode needs owner approval: create OUT/APPROVED_FULL_RUN containing the "
                             "approving commit SHA of the frozen doc 14 (see doc 14 §9)")
        run(out, a.mode, a.shard, a.n_jobs, a.max_units)
    elif a.cmd == "status":
        status(out, a.mode)
    elif a.cmd == "budget":
        budget(out, a.n_jobs)
    else:
        analyse(out, a.mode, Path(a.restored).expanduser() if a.restored else None)


if __name__ == "__main__":
    main()
