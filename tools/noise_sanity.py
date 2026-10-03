"""Pure-noise sanity check for the evaluation protocol (docs/EVALUATION_PROTOCOL.md).

With targets independent of the features, nested CV must report R^2 <= ~0; the v0.x protocol's
training R^2 shows how misleading in-sample fit is.

Usage: PYTHONPATH=src python tools/noise_sanity.py [n_data_seeds]
"""
import sys
import warnings

import numpy as np
import pandas as pd

from topoqspr.evaluate import legacy_protocol, nested_cv
from topoqspr.models import LEGACY_MODELS

warnings.filterwarnings("ignore")


def main(n_seeds: int = 3, n: int = 80, p: int = 60) -> None:
    fam = {k: LEGACY_MODELS[k] for k in ("Ridge", "RandomForest")}
    print("data_seed\tnested_r2_pooled\tnested_r2_fold_mean\tlegacy_cv_r2\tlegacy_train_r2")
    for seed in range(n_seeds):
        rng = np.random.default_rng(seed)
        X = pd.DataFrame(rng.normal(size=(n, p)), columns=[f"f{i}" for i in range(p)])
        y = pd.Series(rng.normal(size=n))
        s = nested_cv(X, y, families=["ridge", "rf"], seeds=[0], ad=None, min_n=10).summary()
        lg = legacy_protocol(X, y, families=fam)
        print(f"{seed}\t{s['r2_pooled_mean']:.3f}\t{s['r2_fold_mean']:.3f}\t{lg.cv_r2:.3f}\t{lg.train_r2:.3f}", flush=True)


if __name__ == "__main__":
    main(*(int(a) for a in sys.argv[1:2]))
