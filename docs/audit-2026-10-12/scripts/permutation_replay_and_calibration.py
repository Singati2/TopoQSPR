"""Audit 2026-10-12, F17 follow-up (recommendations P0 'grouped permutation').

Part 1, replay. Regenerates the *actual* permutations of the archived runs. The pipeline
draws permutation b of a target as permutation_indices(n, groups,
numpy.random.default_rng([seed, b])), with the run record's seed and K = 200, using the CV
groups of that target (connectivity groups restricted to its rows). For every replayed
permutation it checks that each multi-member group's rows were moved as one block onto
another group's rows (no group split). Groups are reconstructed from molecules.csv because
the run does not save them; the reconstruction is checked against the archived fold sizes
(every n_train/n_test of folds.csv is reproduced by nested_cv's splitter).

Part 2, calibration. Empirical type-I error of the permutation test of a simple statistic
(|mean(y | x = 1) - mean(y | x = 0)|, x independent of y) under clustered nulls with
intra-group correlation rho in {0, 0.3, 0.7, 1.0}, for two group structures:
  (a) PR #2's worst case: one unique-size group of 8 plus 8 singletons, x = group indicator;
  (b) the archived structure of mp_bradley: 6 stereo pairs + 699 singletons, x a random
      half split of *groups* (pairs share x, since partners share all 2D features).
Reports the rejection rate at alpha = 0.05 with a Wilson 95% interval.

    python -I docs/audit-2026-10-12/scripts/permutation_replay_and_calibration.py \
        ~/topoqspr_runs/main_confirmatory ~/topoqspr_runs/main_confirmatory_protocol_svr
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from topoqspr.pipeline import permutation_indices


def groups_for(run: Path, target: str) -> np.ndarray:
    M = pd.read_csv(run / "molecules.csv", index_col="id")
    Y = pd.read_csv(run / "targets.csv", index_col="id", float_precision="round_trip")
    ids = Y.index[Y[target].notna() & (M.loc[Y.index, "status"] == "ok")]
    return M.loc[ids, "connectivity_group"].fillna(pd.Series(ids, index=ids)).to_numpy(dtype=object)


def replay(run: Path) -> list[dict]:
    rec = json.loads((run / "run_record.json").read_text())
    perm_cfg = rec["config"]["protocol"]["permutation"]
    seed, K = int(perm_cfg["seed"]), int(perm_cfg["n_permutations"])
    out = []
    for t in rec["config"]["targets"][0]["names"]:
        g = groups_for(run, t).astype(str)
        uniq, inv = np.unique(g, return_inverse=True)
        members = [np.flatnonzero(inv == i) for i in range(uniq.size)]
        multi = [m for m in members if m.size > 1]
        split_groups = moved_blocks = 0
        for b in range(K):
            perm = permutation_indices(g.size, g, np.random.default_rng([seed, b]))
            assert np.array_equal(np.sort(perm), np.arange(g.size))
            for m in multi:
                src = perm[m]  # source rows of the y values placed on this group's rows
                donor = members[inv[src[0]]]
                # a block move: all sources come from one group, and they are that whole group
                if np.unique(inv[src]).size == 1 and np.array_equal(np.sort(src), np.sort(donor)):
                    moved_blocks += 1
                else:
                    split_groups += 1
        out.append(dict(run=run.name, target=t, seed=seed, K=K, n=g.size,
                        multi_groups=len(multi), group_draws=len(multi) * K,
                        moved_as_block=moved_blocks, split=split_groups))
    return out


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def calibrate(structure: str, rho: float, n_sim: int, n_perm: int, rng) -> dict:
    if structure == "worst_case":
        groups = np.array(["big"] * 8 + [f"s{i}" for i in range(8)], dtype=object)
        x = np.r_[np.ones(8), np.zeros(8)]
    else:  # archived mp_bradley structure: 6 pairs + 699 singletons
        groups = np.array([f"p{i // 2}" for i in range(12)] + [f"s{i}" for i in range(699)],
                          dtype=object)
    uniq, inv = np.unique(groups.astype(str), return_inverse=True)
    rej = 0
    for _ in range(n_sim):
        u = rng.normal(size=uniq.size)[inv]  # shared group effect
        y = math.sqrt(rho) * u + math.sqrt(1 - rho) * rng.normal(size=groups.size)
        if structure != "worst_case":
            xg = rng.permutation(np.r_[np.ones(uniq.size // 2), np.zeros(uniq.size - uniq.size // 2)])
            x = xg[inv]
        obs = abs(y[x == 1].mean() - y[x == 0].mean())
        ge = 0
        for _ in range(n_perm):
            yp = y[permutation_indices(groups.size, groups, rng)]
            ge += abs(yp[x == 1].mean() - yp[x == 0].mean()) >= obs
        rej += (1 + ge) / (1 + n_perm) <= 0.05
    lo, hi = wilson(rej, n_sim)
    return dict(structure=structure, rho=rho, n_sim=n_sim, n_perm=n_perm, rejections=rej,
                rate=rej / n_sim, ci_low=lo, ci_high=hi)


def main() -> None:
    rows = [r for a in sys.argv[1:] for r in replay(Path(a).expanduser())]
    rep = pd.DataFrame(rows)
    print(rep.to_string(index=False))
    rng = np.random.default_rng(20261012)
    cal = pd.DataFrame([calibrate(s, rho, 400 if s == "worst_case" else 200, 99, rng)
                        for s in ("worst_case", "archived_mp") for rho in (0.0, 0.3, 0.7, 1.0)])
    print(cal.round(3).to_string(index=False))
    out = Path(__file__).resolve().parents[1] / "out"
    rep.to_csv(out / "permutation_replay.csv", index=False)
    cal.to_csv(out / "permutation_calibration.csv", index=False)


if __name__ == "__main__":
    main()
