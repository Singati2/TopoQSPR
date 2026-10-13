"""Audit 2026-10-12, F17 (PR #2 item A3): is the row-wise fallback of
pipeline.permutation_indices ever triggered in an archived run?

permutation_indices swaps equal-size CV groups as whole blocks, but a group of
size >= 2 whose size no other group shares is shuffled row-wise with the
singletons. On a clustered null this is anti-conservative (PR #2's simulation:
0.48 false-positive rate at nominal 0.05). This script reconstructs, per
target, the CV groups that the permutation stage received (connectivity
groups restricted to the rows with that target, as pipeline.py passes
ctx[t]["groups"]) and lists every group that would be shuffled row-wise.

    python -I docs/audit-2026-10-12/scripts/permutation_block_check.py \
        ~/topoqspr_runs/main_confirmatory ~/topoqspr_runs/main_confirmatory_protocol_svr
"""
import sys
from pathlib import Path

import pandas as pd


def check(run: Path) -> pd.DataFrame:
    M = pd.read_csv(run / "molecules.csv", index_col="id")
    Y = pd.read_csv(run / "targets.csv", index_col="id", float_precision="round_trip")
    rows = []
    for t in Y.columns:
        ids = Y.index[Y[t].notna() & (M.loc[Y.index, "status"] == "ok")]
        g = M.loc[ids, "connectivity_group"].fillna(pd.Series(ids, index=ids))
        sizes = g.value_counts()
        multi = sizes[sizes > 1]
        per_size = multi.value_counts()
        rowwise = [k for k, s in multi.items() if per_size[s] == 1]
        rows.append(dict(run=run.name, target=t, n=len(ids), n_multi_groups=len(multi),
                         rows_in_multi_groups=int(multi.sum()),
                         sizes=dict(sorted(per_size.items())),
                         rowwise_groups=";".join(rowwise),
                         rowwise_rows=int(sum(multi[k] for k in rowwise))))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    out = pd.concat([check(Path(a).expanduser()) for a in sys.argv[1:]])
    print(out.to_string(index=False))
    if out.rowwise_rows.sum():
        sys.exit("row-wise fallback was triggered for at least one target")
