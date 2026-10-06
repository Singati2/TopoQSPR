"""Experiment C (evaluation optimism) table, figure and macros from `run_experiment.py exp-c` outputs.

Usage: python tools/make_expC.py OUT_DIR RUN_DIR [RUN_DIR ...]

Each RUN_DIR holds summary.json, chain.csv (step, r2, ci_low, ci_high, detail) and optimism.csv
(term, value, ci_low, ci_high, pooled).  Writes tab_expC.tex, fig_expC.pdf/.png, expC_macros.tex.
Nothing is typed by hand.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

SERIES = "#2a78d6"   # reference palette slot 1 (validated)
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
TARGET = {"mp_bradley": "melting point (Bradley)", "logp_opera": "log P (OPERA)",
          "logd74_az": "log D$_{7.4}$ (AstraZeneca)", "logs_aqsoldb": "log S (AqSolDB)"}
SHORT = {"mp_bradley": "melting point", "logp_opera": "log P", "logd74_az": "log D", "logs_aqsoldb": "log S"}
STEPS = [("train", "training fit"), ("cv_legacy", "v0.x CV"), ("nested_random", "nested CV"),
         ("nested_scaffold", "scaffold CV")]

plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK_2, "axes.labelcolor": INK, "xtick.color": INK_2,
                     "ytick.color": INK_2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})


def _f(x, nd=2) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "--"
    return "--" if not math.isfinite(v) else f"{v:.{nd}f}"


def _ci(lo, hi) -> str:
    return "" if not (math.isfinite(float(lo)) and math.isfinite(float(hi))) else f" [{_f(lo)}, {_f(hi)}]"


def load(run: Path) -> dict:
    s = json.loads((run / "summary.json").read_text())
    chain = pd.read_csv(run / "chain.csv", keep_default_na=False, na_values=[""]).set_index("step")
    opt = pd.read_csv(run / "optimism.csv", keep_default_na=False, na_values=[""]).set_index("term")
    target = run.name
    info = json.loads((run / "run_info.json").read_text()) if (run / "run_info.json").exists() else {}
    tpath = (info.get("target") or {}).get("path", "")
    for k in TARGET:
        if k in tpath or k == target:
            target = k
    # B_train and B_CV are referenced to the nested CV on the v0.x protocol's complete-case rows
    # when that re-run is a usable R² (|R²| <= 1); when it exploded (extrapolation on the rows it
    # dropped), the all-rows nested CV is the reference and the row is flagged.
    r2_legacy_rows = float(s.get("r2_nested_legacy_rows", math.nan))
    exploded = math.isfinite(r2_legacy_rows) and abs(r2_legacy_rows) > 1 and "nested_random_legacy_rows" in chain.index
    if exploded:
        ref = chain.loc["nested_random"]
        for term, fixed in (("B_train", float(chain.loc["train", "r2"])), ("B_CV", float(chain.loc["cv_legacy", "r2"]))):
            opt.loc[term, ["value", "ci_low", "ci_high", "rows"]] = [
                fixed - float(ref["r2"]), fixed - float(ref["ci_high"]), fixed - float(ref["ci_low"]),
                "all rows (complete-case nested CV exploded)"]
    return dict(target=target, n=int(s["n"]), p=int(s["p"]), chain=chain, opt=opt,
                best=s.get("legacy_best_family", ""), n_legacy=s.get("legacy_n_used"),
                scaffold=s.get("scaffold") or {}, underpowered=bool(s.get("underpowered")),
                ref_all_rows=exploded)


def table(runs: list[dict], out: Path) -> None:
    lines = []
    for r in runs:
        c, o = r["chain"], r["opt"]
        cells = [TARGET.get(r["target"], r["target"]), str(r["n"]), r["best"].replace("_", "\\_")]
        for step, _ in STEPS:
            if step in c.index:
                row = c.loc[step]
                cells.append(_f(row["r2"]) + _ci(row["ci_low"], row["ci_high"]))
            else:
                cells.append("--")
        for term in ("B_train", "B_CV", "B_scaffold"):
            row = o.loc[term]
            cells.append(_f(row["value"]) + _ci(row["ci_low"], row["ci_high"]) +
                         ("$^{\\ast}$" if term != "B_scaffold" and r["ref_all_rows"] else ""))
        lines.append(" & ".join(cells) + " \\\\")
    tex = r"""\begin{table}[t]
\centering\footnotesize
\caption{Experiment C, evaluation optimism on the pooled approved-drug set, topological-index feature set (102 indices of the frozen registry). Training fit: the v0.x protocol's best family refitted and scored on all rows; v0.x CV: that protocol's fold-mean $R^2$ (hyper-parameters tuned on all rows, best of ten families); nested CV: $5\times5$ repeated nested CV with random outer folds (stereoisomers kept together); scaffold CV: the same with Bemis--Murcko scaffold groups as outer folds. $B_{\mathrm{train}}=R^2_{\mathrm{train}}-R^2_{\mathrm{nested}}$, $B_{\mathrm{CV}}=R^2_{\mathrm{CV}}-R^2_{\mathrm{nested}}$, $B_{\mathrm{scaffold}}=R^2_{\mathrm{nested}}-R^2_{\mathrm{scaffold}}$; intervals are corrected repeated-CV $t$ intervals (the training and v0.x CV values are treated as fixed; the scaffold interval combines two standard errors and is approximate). $B_{\mathrm{train}}$ and $B_{\mathrm{CV}}$ are referenced to a nested CV on the v0.x protocol's complete-case rows; $^{\\ast}$: that re-run extrapolated on out-of-domain molecules (pooled $R^2<-1$), so the all-rows nested CV is the reference. Generated by \texttt{tools/make\_expC.py}.}
\label{tab:expC}
\resizebox{\linewidth}{!}{%%
\begin{tabular}{lrlllllllll}
\toprule
Target & $n$ & v0.x best & Training & v0.x CV & Nested CV & Scaffold CV & $B_{\mathrm{train}}$ & $B_{\mathrm{CV}}$ & $B_{\mathrm{scaffold}}$ \\
\midrule
%s
\bottomrule
\end{tabular}}
\end{table}
""" % "\n".join(lines)
    (out / "tab_expC.tex").write_text(tex)


def figure(runs: list[dict], out: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.6, 2.6))
    xs = list(range(len(STEPS)))
    ends = []
    for i, r in enumerate(runs):
        c = r["chain"]
        ys = [float(c.loc[s, "r2"]) if s in c.index else math.nan for s, _ in STEPS]
        ax.plot(xs, ys, color=SERIES, lw=1.4, alpha=0.9, zorder=2)
        ax.scatter(xs, ys, s=16, color=SERIES, zorder=3)
        for s_i, (s, _) in enumerate(STEPS):
            if s in c.index and math.isfinite(float(c.loc[s, "ci_low"])):
                ax.plot([s_i, s_i], [float(c.loc[s, "ci_low"]), float(c.loc[s, "ci_high"])],
                        color=SERIES, lw=1, alpha=0.6, zorder=1)
        ends.append((ys[-1], SHORT.get(r["target"], r["target"])))
    # end labels, pushed apart so that close lines stay readable
    ends.sort()
    placed = []
    for y, lab in ends:
        yy = y
        if placed and yy - placed[-1] < 0.05:
            yy = placed[-1] + 0.05
        placed.append(yy)
        ax.text(xs[-1] + 0.08, yy, lab, fontsize=6.5, color=INK, va="center")
    ax.set_xticks(xs, [lab for _, lab in STEPS])
    ax.set_ylabel("$R^2$ (fold mean)")
    ax.axhline(0, color=INK_2, lw=0.6)
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_xlim(-0.3, len(STEPS) - 0.3)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"fig_expC.{ext}")
    plt.close(fig)


def macros(runs: list[dict], out: Path) -> dict:
    bt = [float(r["opt"].loc["B_train", "value"]) for r in runs]
    bc = [float(r["opt"].loc["B_CV", "value"]) for r in runs]
    bs = [float(r["opt"].loc["B_scaffold", "value"]) for r in runs]
    tr = [float(r["chain"].loc["train", "r2"]) for r in runs if "train" in r["chain"].index]
    ne = [float(r["chain"].loc["nested_random", "r2"]) for r in runs]
    m = {
        "expCnTargets": len(runs),
        "expCminTrain": f"{min(tr):.2f}", "expCmaxTrain": f"{max(tr):.2f}",
        "expCminNested": f"{min(ne):.2f}", "expCmaxNested": f"{max(ne):.2f}",
        "expCminBtrain": f"{min(bt):.2f}", "expCmaxBtrain": f"{max(bt):.2f}",
        "expCminBcv": f"{min(bc):+.2f}", "expCmaxBcv": f"{max(bc):+.2f}",
        "expCminBscaffold": f"{min(bs):+.2f}", "expCmaxBscaffold": f"{max(bs):+.2f}",
        "expCnBcvPositive": sum(1 for b in bc if b > 0),
        "expCnBscaffoldCIexcludesZero": sum(1 for r in runs if float(r["opt"].loc["B_scaffold", "ci_low"]) > 0),
        "expCbestFamilies": ", ".join(sorted({r["best"] for r in runs})).replace("_", "\\_"),
        "expCnRefAllRows": sum(1 for r in runs if r["ref_all_rows"]),
    }
    (out / "expC_macros.tex").write_text("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in m.items()))
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    runs = [load(Path(r)) for r in a.runs]
    table(runs, out)
    figure(runs, out)
    m = macros(runs, out)
    for r in runs:
        c = r["chain"]
        print(f"{r['target']:13s} n={r['n']:3d} best={r['best']:16s} " +
              " ".join(f"{s}={_f(c.loc[s, 'r2'])}" for s, _ in STEPS if s in c.index))
    print(m)


if __name__ == "__main__":
    main()
