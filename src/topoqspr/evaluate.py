"""Leakage-free model evaluation: repeated nested CV and the baseline ladder.

Procedure (see docs/EVALUATION_PROTOCOL.md for the full specification):

For repeat r = 0..R-1 (seed s_r) the rows are split into K outer folds by
``KFold(K, shuffle=True, random_state=s_r)`` (or, with ``groups``,
``GroupKFold(K, shuffle=True, random_state=s_r)``). For each outer fold k:

1. Only the outer-training part D_train^(r,k) is used for selection. It is
   split by ``KFold(L, shuffle=True, random_state=1000*s_r + k)`` (grouped if
   groups are given). For every family and every grid point the inner mean
   squared error is computed; each family's best grid point is recorded and the
   family with the lowest inner MSE (ties: first in the given order) is chosen.
2. A fresh ``clone`` of the chosen family/parameters is fitted on all of
   D_train^(r,k) (preprocessing included, so imputation, column removal and
   scaling are fitted on training rows only).
3. It predicts the untouched outer-test part once.

Outputs: out-of-fold predictions per repeat, per-fold R^2/MAE/RMSE, pooled
per-repeat R^2/MAE/RMSE, and the inner selection table.

``legacy_protocol`` reproduces the v0.x procedure (tuning on all rows, then
``cross_val_score``, best family by the same scores) for Experiment C.

References
----------
Varma S, Simon R (2006) BMC Bioinformatics 7:91. doi:10.1186/1471-2105-7-91
Cawley GC, Talbot NLC (2010) J Mach Learn Res 11:2079-2107.
Krstajic D et al. (2014) J Cheminform 6:10. doi:10.1186/1758-2946-6-10
Bemis GW, Murcko MA (1996) J Med Chem 39:2887-2893. doi:10.1021/jm9602928
"""

from __future__ import annotations

import json
import math
import warnings
from dataclasses import dataclass, field
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.exceptions import ConvergenceWarning, FitFailedWarning
from sklearn.model_selection import (GridSearchCV, GroupKFold, KFold,
                                     cross_val_score)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from . import stats
from .models import (DEFAULT_FAMILIES, LEGACY_MODELS, LEGACY_RANDOM_STATE,
                     ModelFamily, get_family)

__all__ = [
    "DEFAULT_SEEDS",
    "NestedCVResult",
    "ComparisonResult",
    "LegacyResult",
    "make_outer_splits",
    "nested_cv",
    "select_and_refit",
    "compare_feature_sets",
    "legacy_protocol",
    "optimism",
    "scaffold_groups",
    "nested_r2_statistic",
]

DEFAULT_SEEDS: tuple[int, ...] = (0, 1, 2, 3, 4)


# --------------------------------------------------------------------- utils
def _metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    """R^2 = 1 - SSE/SST (SST about mean(y_true)), MAE, RMSE. NaN if undefined."""
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    n = y_true.size
    if n == 0 or not np.all(np.isfinite(y_pred)):
        return dict(r2=math.nan, mae=math.nan, rmse=math.nan)
    err = y_true - y_pred
    sse = float(np.sum(err ** 2))
    sst = float(np.sum((y_true - y_true.mean()) ** 2))
    return dict(r2=1.0 - sse / sst if n >= 2 and sst > 0 else math.nan,
                mae=float(np.mean(np.abs(err))),
                rmse=math.sqrt(sse / n))


def _params_json(params: dict) -> str:
    return json.dumps(params, sort_keys=True, default=repr)


def _as_xy(X, y):
    """Validate shapes; return (X ndarray, y ndarray, row labels, columns)."""
    if isinstance(X, pd.DataFrame):
        labels, columns = X.index, list(X.columns)
        Xa = X.to_numpy(dtype=float)
    else:
        Xa = np.asarray(X, dtype=float)
        if Xa.ndim == 1:
            Xa = Xa[:, None]
        labels, columns = pd.RangeIndex(Xa.shape[0]), list(range(Xa.shape[1]))
    if isinstance(y, pd.Series) and isinstance(X, pd.DataFrame):
        if not y.index.equals(X.index):
            raise ValueError("X and y indices differ; align them first")
    ya = np.asarray(y, dtype=float).ravel()
    if Xa.shape[0] != ya.size:
        raise ValueError(f"X has {Xa.shape[0]} rows, y has {ya.size}")
    if np.isinf(Xa).any():
        raise ValueError("X contains +/-inf; convert to NaN with a recorded reason")
    if not np.all(np.isfinite(ya)):
        raise ValueError("y contains NaN/inf; drop those rows explicitly before evaluation")
    return Xa, ya, labels, columns


def make_outer_splits(n: int, n_splits: int = 5, seeds: Sequence[int] = DEFAULT_SEEDS,
                 groups=None) -> list[tuple[int, int, int, np.ndarray, np.ndarray]]:
    """Outer splits as (repeat, fold, seed, train_idx, test_idx), indices sorted.

    Repeat r uses ``KFold(n_splits, shuffle=True, random_state=seeds[r])`` or,
    with groups, ``GroupKFold(n_splits, shuffle=True, random_state=seeds[r])``
    (scikit-learn >= 1.6: groups are permuted and split into n_splits blocks of
    nearly equal *group* count, so fold sizes in rows can differ).
    """
    out = []
    dummy = np.zeros((n, 1))
    for r, seed in enumerate(seeds):
        if groups is None:
            cv = KFold(n_splits, shuffle=True, random_state=seed)
            it = cv.split(dummy)
        else:
            cv = GroupKFold(n_splits, shuffle=True, random_state=seed)
            it = cv.split(dummy, groups=np.asarray(groups))
        for k, (tr, te) in enumerate(it):
            out.append((r, k, int(seed), np.sort(tr), np.sort(te)))
    return out


def _inner_splits(n_train: int, n_splits: int, seed: int, groups=None):
    dummy = np.zeros((n_train, 1))
    if groups is None:
        return list(KFold(n_splits, shuffle=True, random_state=seed).split(dummy))
    return list(GroupKFold(n_splits, shuffle=True, random_state=seed)
                .split(dummy, groups=np.asarray(groups)))


# ------------------------------------------------------------ inner selection
@dataclass
class Selection:
    """Result of inner-loop selection plus the refitted estimator."""

    estimator: object | None
    family: str | None
    params: dict
    inner_mse: float
    table: pd.DataFrame  # one row per family
    reason: str = ""


def select_and_refit(X: np.ndarray, y: np.ndarray, families: Sequence,
                     inner, random_state: int = 0, groups=None) -> Selection:
    """Choose family + hyperparameters by inner CV on (X, y); refit on all of it.

    ``inner`` is an int (number of shuffled KFold splits, seeded by
    ``random_state``) or an explicit list of (train, test) index pairs.
    Selection criterion: mean inner-fold MSE (equal-size folds: equivalent to
    pooled inner R^2). Grid points whose fit fails in any inner fold score
    NaN and are not selectable. Every fit is on a ``clone``.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    if isinstance(inner, int):
        splits = _inner_splits(len(y), inner, random_state, groups)
    else:
        splits = list(inner)
    n_inner_train = min(len(tr) for tr, _ in splits)
    p = X.shape[1]
    rows = []
    best = (math.inf, None, None, None)  # (mse, family, params, ModelFamily)
    for fam in map(get_family, families):
        ok, why = fam.feasible(n_inner_train, p)
        if not ok:
            rows.append(dict(family=fam.name, status="skipped", reason=why,
                             inner_mse=math.nan, params="", n_configs=0, n_failed=0))
            continue
        grid = fam.param_grid(n_inner_train, p)
        est = clone(fam.build(random_state))
        gs = GridSearchCV(est, [{k: [v] for k, v in g.items()} for g in grid],
                          cv=splits, scoring="neg_mean_squared_error", refit=False,
                          error_score=np.nan, n_jobs=1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            warnings.simplefilter("ignore", FitFailedWarning)
            warnings.filterwarnings("ignore", message=".*test scores are non-finite.*")
            warnings.filterwarnings("ignore", message=".*fits failed.*")
            try:
                gs.fit(X, y)
            except ValueError as exc:  # scikit-learn raises when *every* fit fails
                rows.append(dict(family=fam.name, status="failed",
                                 reason=f"every inner fit failed: {str(exc)[:200]}",
                                 inner_mse=math.nan, params="", n_configs=len(grid),
                                 n_failed=len(grid)))
                continue
        scores = np.asarray(gs.cv_results_["mean_test_score"], dtype=float)
        n_failed = int(np.sum(~np.isfinite(scores)))
        if n_failed == scores.size:
            rows.append(dict(family=fam.name, status="failed",
                             reason="every grid point failed in inner CV",
                             inner_mse=math.nan, params="", n_configs=scores.size,
                             n_failed=n_failed))
            continue
        i = int(np.nanargmax(scores))  # first maximum -> deterministic
        params = gs.cv_results_["params"][i]
        mse = float(-scores[i])
        rows.append(dict(family=fam.name, status="ok", reason="", inner_mse=mse,
                         params=_params_json(params), n_configs=scores.size,
                         n_failed=n_failed))
        if mse < best[0]:
            best = (mse, fam.name, params, fam)
    table = pd.DataFrame(rows, columns=["family", "status", "reason", "inner_mse",
                                        "params", "n_configs", "n_failed"])
    if best[1] is None:
        return Selection(None, None, {}, math.nan, table, "no family could be fitted")
    mse, name, params, fam = best
    est = clone(fam.build(random_state)).set_params(**params)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        try:
            est.fit(X, y)
        except Exception as exc:  # noqa: BLE001 - reported, never imputed
            return Selection(None, name, dict(params), mse, table,
                             f"refit of {name} failed: {exc!r}"[:300])
    return Selection(est, name, dict(params), mse, table)


# ---------------------------------------------------------------- nested CV
@dataclass
class NestedCVResult:
    """Output of :func:`nested_cv`.

    predictions : one row per (repeat, row): repeat, fold, row (label), pos
        (0-based position), y_true, y_pred, in_domain (AD flag; NaN if no AD).
    folds : one row per outer fold: repeat, fold, seed, n_train, n_test,
        family, params, inner_mse, r2, mae, rmse, ad_coverage, status, reason.
    repeats : one row per repeat: repeat, seed, n, r2, mae, rmse (pooled over
        that repeat's out-of-fold predictions).
    inner : the per-family inner selection table of every outer fold.
    """

    predictions: pd.DataFrame
    folds: pd.DataFrame
    repeats: pd.DataFrame
    inner: pd.DataFrame
    meta: dict = field(default_factory=dict)
    models: list | None = None

    def summary(self) -> dict:
        f, r = self.folds, self.repeats
        return dict(
            r2_fold_mean=float(f["r2"].mean()), r2_fold_sd=float(f["r2"].std(ddof=1)),
            r2_fold_median=float(f["r2"].median()),
            mae_fold_mean=float(f["mae"].mean()), rmse_fold_mean=float(f["rmse"].mean()),
            r2_pooled_mean=float(r["r2"].mean()), r2_pooled_sd=float(r["r2"].std(ddof=1)),
            mae_pooled_mean=float(r["mae"].mean()),
            n_folds=int(len(f)), n_failed_folds=int((f["status"] != "ok").sum()),
            family_counts=f["family"].value_counts(dropna=False).sort_index().to_dict(),
        )


def _run_fold(X, y, split, families, inner_splits, groups, random_state, ad,
              keep_model):
    r, k, seed, tr, te = split
    g_tr = None if groups is None else np.asarray(groups)[tr]
    inner = _inner_splits(len(tr), inner_splits, 1000 * seed + k, g_tr)
    sel = select_and_refit(X[tr], y[tr], families, inner, random_state)
    tab = sel.table.copy()
    tab.insert(0, "fold", k)
    tab.insert(0, "repeat", r)
    in_dom = np.full(te.size, np.nan)
    ad_cov = math.nan
    if ad is not None and X.shape[1] > 0:
        try:
            dom = stats.ApplicabilityDomain(method=ad).fit(X[tr])
            in_dom = dom.in_domain(X[te]).astype(float)
            ad_cov = float(np.mean(in_dom))
        except ValueError:
            pass
    if sel.estimator is None:
        y_pred = np.full(te.size, np.nan)
        status, reason = "failed", sel.reason
    else:
        try:
            y_pred = np.asarray(sel.estimator.predict(X[te]), dtype=float).ravel()
            status, reason = ("ok", "") if np.all(np.isfinite(y_pred)) else \
                ("failed", "non-finite predictions")
        except Exception as exc:  # noqa: BLE001 - reported, never imputed
            y_pred = np.full(te.size, np.nan)
            status, reason = "failed", f"predict failed: {exc!r}"[:300]
    m = _metrics(y[te], y_pred)
    fold_row = dict(repeat=r, fold=k, seed=seed, n_train=int(tr.size), n_test=int(te.size),
                    family=sel.family, params=_params_json(sel.params),
                    inner_mse=sel.inner_mse, **m, ad_coverage=ad_cov,
                    status=status, reason=reason)
    pred = pd.DataFrame(dict(repeat=r, fold=k, pos=te, y_true=y[te], y_pred=y_pred,
                             in_domain=in_dom))
    return fold_row, pred, tab, (sel.estimator if keep_model else None)


def nested_cv(X, y, families: Sequence[str | ModelFamily] = DEFAULT_FAMILIES, *,
              outer_splits: int = 5, seeds: Sequence[int] = DEFAULT_SEEDS,
              inner_splits: int = 3, groups=None, random_state: int = 0,
              ad: str | None = "knn", keep_models: bool = False,
              n_jobs: int = 1) -> NestedCVResult:
    """Repeated nested cross-validation with family + hyperparameter selection.

    Parameters
    ----------
    X : DataFrame or array (n, p); NaN allowed (median-imputed inside folds).
        p = 0 is allowed and evaluates only the ``mean`` family (null model).
    y : array-like (n,), finite.
    families : family names or :class:`ModelFamily` objects, in tie-break order.
    outer_splits : K, outer folds per repeat (default 5).
    seeds : one seed per repeat (default 0..4 -> 5 x 5 CV).
    inner_splits : L, inner folds (default 3).
    groups : optional group labels (e.g. :func:`scaffold_groups`); outer and
        inner splits then keep each group within one fold.
    random_state : seed passed to every estimator (fixed across folds).
    ad : "knn", "leverage" or None - applicability domain fitted per outer
        training part; flags outer-test rows.
    keep_models : keep the 25 refitted estimators in ``result.models``.
    n_jobs : outer folds in parallel (joblib); results are identical.
    """
    Xa, ya, labels, columns = _as_xy(X, y)
    families = ("mean",) if Xa.shape[1] == 0 else tuple(families)
    splits = make_outer_splits(len(ya), outer_splits, seeds, groups)
    args = (Xa, ya)
    if n_jobs == 1:
        outs = [_run_fold(*args, s, families, inner_splits, groups, random_state, ad,
                          keep_models) for s in splits]
    else:
        from joblib import Parallel, delayed
        outs = Parallel(n_jobs=n_jobs)(
            delayed(_run_fold)(*args, s, families, inner_splits, groups, random_state,
                               ad, keep_models) for s in splits)
    folds = pd.DataFrame([o[0] for o in outs])
    preds = pd.concat([o[1] for o in outs], ignore_index=True)
    preds.insert(2, "row", np.asarray(labels)[preds["pos"].to_numpy()])
    preds = preds.sort_values(["repeat", "pos"], kind="mergesort").reset_index(drop=True)
    inner = pd.concat([o[2] for o in outs], ignore_index=True)
    reps = []
    for r, g in preds.groupby("repeat", sort=True):
        m = _metrics(g["y_true"].to_numpy(), g["y_pred"].to_numpy())
        reps.append(dict(repeat=int(r), seed=int(seeds[int(r)]), n=int(len(g)), **m))
    meta = dict(n=int(len(ya)), p=int(Xa.shape[1]), columns=columns,
                families=[get_family(f).name for f in families],
                outer_splits=outer_splits, seeds=list(seeds), inner_splits=inner_splits,
                grouped=groups is not None, random_state=random_state, ad=ad,
                n_missing_cells=int(np.isnan(Xa).sum()))
    return NestedCVResult(preds, folds, pd.DataFrame(reps), inner, meta,
                          [o[3] for o in outs] if keep_models else None)


def nested_r2_statistic(**nested_kwargs):
    """Return ``score_fn(X, y)`` = mean pooled R^2 over repeats of
    :func:`nested_cv` (for :func:`topoqspr.stats.y_scrambling`)."""
    def score_fn(X, y):
        return nested_cv(X, y, **nested_kwargs).summary()["r2_pooled_mean"]
    return score_fn


# ------------------------------------------------------- feature comparison
@dataclass
class ComparisonResult:
    """Output of :func:`compare_feature_sets`.

    results : {set name: NestedCVResult} (identical outer folds).
    paired  : per (set, repeat, fold): d_r2 = R2_set - R2_ref,
              d_mae = MAE_ref - MAE_set, d_rmse = RMSE_ref - RMSE_set
              (positive = the set beats the reference).
    summary : per set: fold-mean and pooled R2/MAE, and for d_r2 / d_mae the
              mean, corrected-t CI, p, plus BH q across the non-reference sets.
    """

    results: dict
    paired: pd.DataFrame
    summary: pd.DataFrame
    reference: str


def compare_feature_sets(feature_sets: Mapping[str, pd.DataFrame], y, *,
                         reference: str = "size",
                         families: Sequence[str | ModelFamily] = DEFAULT_FAMILIES,
                         alpha: float = 0.05, **nested_kwargs) -> ComparisonResult:
    """Run :func:`nested_cv` on every feature set with identical outer folds.

    All frames must share one index equal to ``y.index`` (when y is a Series).
    A zero-column set is evaluated with the null ``mean`` family. Paired
    per-fold differences against ``reference`` get a corrected repeated-CV
    t-test CI (:func:`topoqspr.stats.corrected_ttest`); p-values of d_r2 are
    BH-adjusted across the compared sets.
    """
    if reference not in feature_sets:
        raise KeyError(f"reference set {reference!r} not in feature_sets")
    names = list(feature_sets)
    idx = feature_sets[reference].index
    for nm in names:
        if not feature_sets[nm].index.equals(idx):
            raise ValueError(f"feature set {nm!r} index differs from {reference!r}")
    y = pd.Series(np.asarray(y, dtype=float).ravel(), index=idx) \
        if not isinstance(y, pd.Series) else y
    results = {nm: nested_cv(feature_sets[nm], y, families, **nested_kwargs) for nm in names}
    keys = ["repeat", "fold"]
    ref = results[reference].folds.set_index(keys)
    paired_rows, summ = [], []
    for nm in names:
        f = results[nm].folds.set_index(keys)
        if not (f[["n_train", "n_test"]].equals(ref[["n_train", "n_test"]])):
            raise RuntimeError("outer folds differ between feature sets")  # defensive
        d = pd.DataFrame(dict(set=nm, d_r2=f["r2"] - ref["r2"], d_mae=ref["mae"] - f["mae"],
                              d_rmse=ref["rmse"] - f["rmse"], n_train=f["n_train"],
                              n_test=f["n_test"])).reset_index()
        paired_rows.append(d)
        s = results[nm].summary()
        row = dict(set=nm, p=results[nm].meta["p"], r2_fold_mean=s["r2_fold_mean"],
                   r2_pooled_mean=s["r2_pooled_mean"], mae_fold_mean=s["mae_fold_mean"],
                   n_failed_folds=s["n_failed_folds"])
        if nm != reference:
            for met in ("d_r2", "d_mae"):
                t = stats.corrected_ttest(d[met], d["n_train"], d["n_test"], alpha=alpha)
                row.update({f"{met}_mean": t["mean"], f"{met}_ci_low": t["ci_low"],
                            f"{met}_ci_high": t["ci_high"], f"{met}_p": t["p"]})
        summ.append(row)
    summary = pd.DataFrame(summ)
    if "d_r2_p" in summary:
        q, _ = stats.benjamini_hochberg(summary["d_r2_p"].to_numpy(), q=alpha)
        summary["d_r2_q"] = q
    paired = pd.concat(paired_rows, ignore_index=True)
    return ComparisonResult(results, paired, summary, reference)


# ------------------------------------------------------------ v0.x protocol
@dataclass
class LegacyResult:
    """Output of :func:`legacy_protocol` (v0.x ``run_ml_qspr`` reproduction).

    table : per family R2_mean, R2_std (population SD, as v0.x), MAE_mean,
            params (tuned on all rows).
    best_family, cv_r2 : the family with the largest R2_mean and that value.
    train_r2 : R^2 of the best pipeline refitted on, and predicting, all rows.
    n_used, dropped : rows used / row labels dropped (non-finite X or y).
    """

    table: pd.DataFrame
    best_family: str | None
    cv_r2: float
    train_r2: float
    n_used: int
    dropped: list


def legacy_protocol(X, y, families: Mapping | None = None, cv_folds: int = 5,
                    random_state: int = LEGACY_RANDOM_STATE) -> LegacyResult:
    """Reproduce the v0.x evaluation exactly (for Experiment C).

    v0.x (ml.py ``run_ml_qspr``): rows with any non-finite feature or target
    are dropped; for each family, ``GridSearchCV(Pipeline([scaler, model]),
    grid, cv=KFold(3, shuffle, 42), scoring="r2")`` is fitted on **all**
    rows; the tuned pipeline is scored with ``cross_val_score`` on
    ``KFold(5, shuffle, 42)`` (R^2 and MAE); the family with the largest mean
    fold R^2 is reported. This is *not* nested: hyperparameters and the family
    are chosen with the rows later used as test folds. Families with no grid
    (OLS) skip the grid search. The best pipeline is refitted on all rows and
    its training R^2 is returned. Unlike v0.x, every fit uses ``clone``.
    """
    families = LEGACY_MODELS if families is None else families
    Xa = np.asarray(X, dtype=float)
    if Xa.ndim == 1:
        Xa = Xa[:, None]
    ya = np.asarray(y, dtype=float).ravel()
    labels = X.index if isinstance(X, pd.DataFrame) else pd.RangeIndex(len(ya))
    ok = np.all(np.isfinite(Xa), axis=1) & np.isfinite(ya)
    dropped = list(np.asarray(labels)[~ok])
    Xa, ya = Xa[ok], ya[ok]
    if len(ya) < 10:
        raise ValueError(f"v0.x skips targets with < 10 rows (got {len(ya)})")
    outer = KFold(max(2, min(cv_folds, len(ya))), shuffle=True, random_state=random_state)
    inner = KFold(max(2, min(3, len(ya))), shuffle=True, random_state=random_state)
    rows, best_r2, best = [], -np.inf, None
    for name, (factory, grid) in families.items():
        pipe = Pipeline([("scaler", StandardScaler()), ("model", factory())])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                if grid:
                    gs = GridSearchCV(clone(pipe), grid, cv=inner, scoring="r2",
                                      error_score=np.nan, n_jobs=1).fit(Xa, ya)
                    tuned, params = clone(gs.best_estimator_), gs.best_params_
                else:
                    tuned, params = clone(pipe), {}
                cv_r2 = cross_val_score(tuned, Xa, ya, cv=outer, scoring="r2",
                                        error_score=np.nan)
                cv_mae = cross_val_score(tuned, Xa, ya, cv=outer,
                                         scoring="neg_mean_absolute_error",
                                         error_score=np.nan)
            except Exception as exc:  # noqa: BLE001 - v0.x skipped failing families
                rows.append(dict(family=name, R2_mean=math.nan, R2_std=math.nan,
                                 MAE_mean=math.nan, params="", reason=repr(exc)))
                continue
        mean_r2 = float(np.nanmean(cv_r2))
        rows.append(dict(family=name, R2_mean=mean_r2, R2_std=float(np.nanstd(cv_r2)),
                         MAE_mean=float(-np.nanmean(cv_mae)),
                         params=_params_json(params), reason=""))
        if mean_r2 > best_r2:
            best_r2, best = mean_r2, (name, tuned)
    train_r2 = math.nan
    if best is not None:
        final = clone(best[1])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            final.fit(Xa, ya)
        train_r2 = _metrics(ya, final.predict(Xa))["r2"]
    return LegacyResult(pd.DataFrame(rows), best[0] if best else None,
                        float(best_r2) if best else math.nan, train_r2,
                        int(len(ya)), dropped)


def optimism(legacy: LegacyResult, nested: NestedCVResult) -> dict:
    """Experiment C optimism terms.

    B_train = R2_train(legacy best) - R2_nested
    B_CV    = R2_CV(legacy best)    - R2_nested
    with R2_nested = mean outer-fold R^2 (same estimator type as the v0.x fold
    mean); the pooled variant (mean over repeats of pooled OOF R^2) is also
    returned as *_pooled.
    """
    s = nested.summary()
    return dict(r2_train=legacy.train_r2, r2_cv_legacy=legacy.cv_r2,
                r2_nested=s["r2_fold_mean"], r2_nested_pooled=s["r2_pooled_mean"],
                B_train=legacy.train_r2 - s["r2_fold_mean"],
                B_CV=legacy.cv_r2 - s["r2_fold_mean"],
                B_train_pooled=legacy.train_r2 - s["r2_pooled_mean"],
                B_CV_pooled=legacy.cv_r2 - s["r2_pooled_mean"])


# ------------------------------------------------------------------- groups
def scaffold_groups(mols, acyclic: str = "singleton") -> np.ndarray:
    """Bemis-Murcko scaffold SMILES per molecule, for grouped outer folds.

    The scaffold is RDKit's ``MurckoScaffoldSmiles`` (ring systems + linkers,
    atom and bond types kept; Bemis & Murcko 1996). Acyclic molecules have an
    empty scaffold: ``acyclic="singleton"`` (default) gives each its own group
    (they share no ring framework, so no scaffold leakage), ``"shared"`` puts
    them all in one group. ``None`` molecules raise ``ValueError``.
    """
    from rdkit.Chem.Scaffolds import MurckoScaffold

    if acyclic not in ("singleton", "shared"):
        raise ValueError("acyclic must be 'singleton' or 'shared'")
    out = []
    for i, mol in enumerate(mols):
        if mol is None:
            raise ValueError(f"molecule {i} is None")
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        if scaf == "":
            scaf = f"__acyclic_{i}" if acyclic == "singleton" else "__acyclic"
        out.append(scaf)
    return np.asarray(out, dtype=object)
