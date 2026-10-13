"""Counterexample to blanket claims of permutation validity under clustered outcomes.

Run: python tools/audit_permutation_exchangeability.py

This is an *adversarial null simulation*, not a reanalysis of the four archived
drug endpoints. The null has y independent of the predictor, but within-group
outcomes are correlated. A unique-size group cannot be swapped as a block by
the current permutation rule and is consequently shuffled row-wise.

Rejecting far above nominal alpha demonstrates that 'every draw is a
permutation of y' is insufficient to establish valid p-values for *arbitrary*
grouped null distributions. Assess exchangeability for the actual datasets.
"""
from __future__ import annotations

import numpy as np

from topoqspr.pipeline import permutation_indices


def main() -> None:
    groups = ["large"] * 8 + [f"singleton{i}" for i in range(8)]
    rng = np.random.default_rng(2212)
    n_sim, n_perm = 200, 99
    rejections = 0

    for _ in range(n_sim):
        # True null: cluster-specific random intercept is independent of x.
        # x_i = 1 for the large group, 0 for singleton observations.
        y = np.r_[np.repeat(rng.normal(), 8), rng.normal(size=8)]
        observed = abs(float(y[:8].mean() - y[8:].mean()))
        more_extreme = 0
        for _ in range(n_perm):
            perm = permutation_indices(16, groups, rng)
            yp = y[perm]
            simulated = abs(float(yp[:8].mean() - yp[8:].mean()))
            more_extreme += simulated >= observed
        p = (1 + more_extreme) / (n_perm + 1)
        rejections += p <= 0.05

    rate = rejections / n_sim
    print(f"Seed 2212 | {n_sim} independent null data sets | {n_perm} permutations")
    print(f"Rejections: {rejections}/{n_sim}; empirical false-positive rate: {rate:.3f}")
    print("Nominal rate: 0.050. This is a null-exchangeability counterexample, "
          "NOT an observed false-positive rate for TopoQSPR's drug endpoints.")
    if rate <= 0.15:
        raise AssertionError("Diagnostic did not reproduce expected anti-conservatism; investigate")


if __name__ == "__main__":
    main()
