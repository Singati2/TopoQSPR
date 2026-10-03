"""Registry of regressor families with small hyperparameter grids.

Each :class:`ModelFamily` builds a complete, *unfitted* scikit-learn
estimator (preprocessing + model) and declares its grid. Callers must fit
``sklearn.base.clone(family.build(seed))``, never a shared instance; the
evaluation module does this for every fit.

Preprocessing (inside every pipeline, therefore fitted on the training part
of each split only):

1. ``SimpleImputer(strategy="median")`` - missing feature values are replaced
   by the training-fold median. Columns that are entirely missing in the
   training fold are dropped by the imputer.
2. :class:`DropConstantAndDuplicate` - removes columns that are constant in
   the training fold and exact duplicates of an earlier column (first kept).
3. ``StandardScaler`` - for the scale-sensitive families (ridge, elastic net,
   OLS, PLS, SVR, kNN). Tree ensembles are not scaled.

Elastic net and SVR also standardise the *target* (``TransformedTargetRegressor``)
because their penalty (alpha * l1_ratio) and epsilon-tube are in units of y;
without it the same grid means different things for MolWt (~400) and logP (~2).

The ``linear`` (OLS) family is feasible only when n_train > p + 1 (an intercept
plus p coefficients must be identifiable); otherwise it is skipped with a
reason. ``mean`` (``DummyRegressor``) is the null model M0: y_hat = mean(y_train).

``LEGACY_MODELS`` reproduces the v0.x ``ml.MODELS_AND_GRIDS`` (same
estimators, grids and seeds) for the Experiment C comparison only.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Mapping, Sequence

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import TransformedTargetRegressor
from sklearn.cross_decomposition import PLSRegression
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import (ExtraTreesRegressor, GradientBoostingRegressor,
                              RandomForestRegressor)
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet, Lasso, LinearRegression, Ridge
from sklearn.model_selection import ParameterGrid
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.utils.validation import check_is_fitted
from xgboost import XGBRegressor

__all__ = [
    "DropConstantAndDuplicate",
    "ModelFamily",
    "FAMILIES",
    "DEFAULT_FAMILIES",
    "FAST_FAMILIES",
    "get_family",
    "LEGACY_MODELS",
    "LEGACY_RANDOM_STATE",
]


class DropConstantAndDuplicate(TransformerMixin, BaseEstimator):
    """Drop training-constant columns and exact duplicate columns.

    A column j is constant if max(x_j) == min(x_j) over the training rows. A
    column is a duplicate if it is bitwise equal to an earlier kept column.
    Raises ``ValueError`` at fit time if no column remains.
    """

    def fit(self, X, y=None):
        X = np.asarray(X, dtype=float)
        keep: list[int] = []
        seen: dict[bytes, int] = {}
        for j in range(X.shape[1]):
            col = X[:, j]
            if X.shape[0] == 0 or np.nanmax(col) == np.nanmin(col):
                continue
            key = np.ascontiguousarray(col).tobytes()
            if key in seen:
                continue
            seen[key] = j
            keep.append(j)
        if not keep:
            raise ValueError("no non-constant feature in the training fold")
        self.keep_ = np.asarray(keep, dtype=int)
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X):
        check_is_fitted(self, "keep_")
        return np.asarray(X, dtype=float)[:, self.keep_]


def _prep(scale: bool) -> list[tuple[str, object]]:
    steps: list[tuple[str, object]] = [
        ("impute", SimpleImputer(strategy="median")),
        ("clean", DropConstantAndDuplicate()),
    ]
    if scale:
        steps.append(("scale", StandardScaler()))
    return steps


def _pipe(model, scale: bool) -> Pipeline:
    return Pipeline(_prep(scale) + [("model", model)])


def _scaled_target(pipe: Pipeline) -> TransformedTargetRegressor:
    return TransformedTargetRegressor(regressor=pipe, transformer=StandardScaler())


@dataclass(frozen=True)
class ModelFamily:
    """A regressor family: an estimator factory plus a hyperparameter grid.

    ``factory(random_state)`` returns a new unfitted estimator.
    ``grid`` maps parameter names (as accepted by ``set_params`` on the built
    estimator) to candidate values. ``constraint(params, n_train, p)`` returns
    False for grid points that cannot be fitted with ``n_train`` training rows
    and ``p`` input columns. ``min_rows(p)`` is the smallest training size at
    which the family is used at all.
    """

    name: str
    factory: Callable[[int], BaseEstimator]
    grid: Mapping[str, Sequence] = field(default_factory=dict)
    constraint: Callable[[dict, int, int], bool] | None = None
    min_rows: Callable[[int], int] | None = None
    needs_features: bool = True

    def build(self, random_state: int = 0) -> BaseEstimator:
        return self.factory(random_state)

    def feasible(self, n_train: int, p: int) -> tuple[bool, str]:
        if self.needs_features and p == 0:
            return False, "no features"
        if self.min_rows is not None and n_train < self.min_rows(p):
            return False, f"needs n_train >= {self.min_rows(p)} for p={p}, got {n_train}"
        if not self.param_grid(n_train, p):
            return False, "no feasible grid point"
        return True, ""

    def param_grid(self, n_train: int, p: int) -> list[dict]:
        points = list(ParameterGrid(dict(self.grid))) if self.grid else [{}]
        if self.constraint is not None:
            points = [g for g in points if self.constraint(g, n_train, p)]
        return points


def _f(name, factory, grid=None, constraint=None, min_rows=None, needs_features=True):
    return ModelFamily(name, factory, grid or {}, constraint, min_rows, needs_features)


FAMILIES: dict[str, ModelFamily] = {
    f.name: f
    for f in [
        _f("mean", lambda rs: DummyRegressor(strategy="mean"), needs_features=False),
        _f("linear", lambda rs: _pipe(LinearRegression(), scale=True),
           min_rows=lambda p: p + 2),
        _f("ridge", lambda rs: _pipe(Ridge(), scale=True),
           {"model__alpha": [0.01, 0.1, 1.0, 10.0, 100.0]}),
        _f("elasticnet",
           lambda rs: _scaled_target(_pipe(ElasticNet(max_iter=20000, random_state=rs),
                                           scale=True)),
           {"regressor__model__alpha": [0.001, 0.01, 0.1, 1.0],
            "regressor__model__l1_ratio": [0.2, 0.5, 0.8]}),
        _f("pls", lambda rs: _pipe(PLSRegression(scale=False), scale=True),
           {"model__n_components": [1, 2, 3, 5, 10]},
           constraint=lambda g, n, p: g["model__n_components"] <= min(p, n - 1)),
        _f("knn", lambda rs: _pipe(KNeighborsRegressor(), scale=True),
           {"model__n_neighbors": [3, 5, 10], "model__weights": ["uniform", "distance"]},
           constraint=lambda g, n, p: g["model__n_neighbors"] <= n),
        _f("svr", lambda rs: _scaled_target(_pipe(SVR(kernel="rbf", gamma="scale"),
                                                  scale=True)),
           {"regressor__model__C": [0.3, 3.0, 30.0],
            "regressor__model__epsilon": [0.1]}),
        _f("rf", lambda rs: _pipe(RandomForestRegressor(
               n_estimators=200, random_state=rs, n_jobs=1), scale=False),
           {"model__max_features": [1.0, 0.33], "model__min_samples_leaf": [1, 5]}),
        _f("extratrees", lambda rs: _pipe(ExtraTreesRegressor(
               n_estimators=200, random_state=rs, n_jobs=1), scale=False),
           {"model__max_features": [1.0, 0.33], "model__min_samples_leaf": [1, 5]}),
        _f("gbr", lambda rs: _pipe(GradientBoostingRegressor(
               n_estimators=200, random_state=rs), scale=False),
           {"model__learning_rate": [0.05, 0.1], "model__max_depth": [2, 3]}),
        _f("xgb", lambda rs: _pipe(XGBRegressor(
               n_estimators=200, random_state=rs, n_jobs=1, tree_method="hist",
               verbosity=0), scale=False),
           {"model__learning_rate": [0.05, 0.1], "model__max_depth": [3, 6]}),
    ]
}

#: Families used by default in nested CV (``mean`` is used for the null set).
DEFAULT_FAMILIES: tuple[str, ...] = (
    "ridge", "elasticnet", "pls", "linear", "knn", "svr",
    "rf", "extratrees", "gbr", "xgb",
)
#: A cheap subset for smoke tests and quick looks.
FAST_FAMILIES: tuple[str, ...] = ("ridge", "pls", "knn")


def get_family(f: str | ModelFamily) -> ModelFamily:
    if isinstance(f, ModelFamily):
        return f
    try:
        return FAMILIES[f]
    except KeyError:
        raise KeyError(f"unknown model family {f!r}; known: {sorted(FAMILIES)}") from None


# --- v0.x reproduction (Experiment C only) ----------------------------------
LEGACY_RANDOM_STATE = 42

#: name -> (factory, grid) exactly as ml.MODELS_AND_GRIDS in v0.x, wrapped
#: in Pipeline([("scaler", StandardScaler()), ("model", model)]). n_jobs is 1
#: (v0.x used -1; seeded results are identical, and runs are reproducible).
LEGACY_MODELS: dict[str, tuple[Callable[[], BaseEstimator], dict]] = {
    "LinearReg": (lambda: LinearRegression(), {}),
    "Ridge": (lambda: Ridge(random_state=LEGACY_RANDOM_STATE),
              {"model__alpha": [0.1, 1.0, 10.0]}),
    "Lasso": (lambda: Lasso(random_state=LEGACY_RANDOM_STATE, max_iter=3000),
              {"model__alpha": [0.01, 0.1, 1.0]}),
    "ElasticNet": (lambda: ElasticNet(random_state=LEGACY_RANDOM_STATE, max_iter=3000),
                   {"model__alpha": [0.1, 1.0], "model__l1_ratio": [0.3, 0.7]}),
    "RandomForest": (lambda: RandomForestRegressor(random_state=LEGACY_RANDOM_STATE, n_jobs=1),
                     {"model__n_estimators": [100], "model__max_depth": [None, 10]}),
    "ExtraTrees": (lambda: ExtraTreesRegressor(random_state=LEGACY_RANDOM_STATE, n_jobs=1),
                   {"model__n_estimators": [100], "model__max_depth": [None, 10]}),
    "GradientBoosting": (lambda: GradientBoostingRegressor(random_state=LEGACY_RANDOM_STATE),
                         {"model__n_estimators": [100], "model__learning_rate": [0.1],
                          "model__max_depth": [3, 5]}),
    "XGBoost": (lambda: XGBRegressor(random_state=LEGACY_RANDOM_STATE, verbosity=0, n_jobs=1),
                {"model__n_estimators": [100], "model__learning_rate": [0.1],
                 "model__max_depth": [3, 6]}),
    "SVR": (lambda: SVR(),
            {"model__C": [1.0, 10.0], "model__kernel": ["rbf"], "model__gamma": ["scale"]}),
    "NeuralNet": (lambda: MLPRegressor(max_iter=500, random_state=LEGACY_RANDOM_STATE),
                  {"model__hidden_layer_sizes": [(64, 32)], "model__alpha": [0.001]}),
}
