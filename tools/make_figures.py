"""Paper figures generated only from a run's result files (no hand-copied numbers).

Usage: python tools/make_figures.py RESULTS_DIR OUT_DIR [--learner selected]

Writes, for the run in RESULTS_DIR (summary.csv from `topoqspr run`):
  fig_ladder.pdf/.png  pooled R² per feature set (baseline ladder), one panel per target
  fig_delta.pdf/.png   paired ΔR² versus the size baseline with 95 % corrected-t CIs

Design (dataviz method, print figure, light theme): one series per panel, so a single colour
(reference palette slot 1) and no legend; text in ink colours; reference lines recessive;
underpowered results drawn hollow and labelled, so meaning never relies on colour alone.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

SERIES = "#2a78d6"        # palette slot 1 (blue), validated reference palette
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#d9d8d4"
LADDER = ["null", "size", "size_rings", "composition", "TI", "TI_composition", "rdkit2d"]
LABEL = {"null": "null (mean)", "size": "size", "size_rings": "size + rings",
         "composition": "composition", "TI": "topological indices",
         "TI_composition": "indices + composition", "rdkit2d": "RDKit 2D"}

plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK_2, "axes.labelcolor": INK,
                     "xtick.color": INK_2, "ytick.color": INK_2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150, "savefig.bbox": "tight"})


def _flag(x) -> bool:
    return str(x).strip().lower() in ("true", "1")


def ladder(summary: pd.DataFrame, learner: str, out: Path) -> None:
    targets = list(dict.fromkeys(summary["target"]))
    ncol = 2 if len(targets) > 1 else 1
    nrow = math.ceil(len(targets) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.9 * ncol, 2.4 * nrow), squeeze=False,
                             gridspec_kw=dict(wspace=0.75, hspace=0.75))
    for ax, t in zip(axes.flat, targets):
        sub = summary[(summary["target"] == t) & (summary["learner"] == learner)]
        sets = [s for s in LADDER if s in set(sub["set"])]
        ys = list(range(len(sets)))[::-1]
        for y, s in zip(ys, sets):
            r = sub[sub["set"] == s].iloc[0]
            v = r["r2_pooled_mean"]
            if pd.isna(v):
                ax.text(0, y, "undefined", va="center", fontsize=7, color=INK_2)
                continue
            hollow = _flag(r["underpowered"])
            ax.plot([0, v], [y, y], color=GRID, lw=1, zorder=1)
            ax.scatter([v], [y], s=36, zorder=3, color="white" if hollow else SERIES,
                       edgecolors=SERIES, linewidths=1.5)
            ax.text(v, y + 0.28, f"{v:.2f}", ha="center", fontsize=6.5, color=INK)
        ax.set_yticks(ys, [LABEL[s] for s in sets])
        ax.axvline(0, color=INK_2, lw=0.6)
        n = int(sub["n"].iloc[0]) if len(sub) else 0
        ax.set_title(f"{t}  (n = {n})", fontsize=8, color=INK, loc="left")
        ax.set_xlabel("pooled cross-validated R²")
        ax.grid(axis="x", color=GRID, lw=0.5)
    for ax in list(axes.flat)[len(targets):]:
        ax.axis("off")
    fig.text(0, -0.02, f"Learner: {learner}. Hollow marker: underpowered (n below the plan's minimum).",
             fontsize=6.5, color=INK_2)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"fig_ladder.{ext}")
    plt.close(fig)


def delta(summary: pd.DataFrame, learner: str, out: Path, threshold: float) -> None:
    sub = summary[(summary["learner"] == learner) & summary["d_r2_mean"].notna()]
    sub = sub[sub["set"] != sub["reference"]]
    keep = ("composition", "TI", "TI_composition", "rdkit2d")  # confirmatory/secondary comparisons only
    rows = [(r["target"], r["set"], r) for _, r in sub.iterrows() if r["set"] in keep]
    lim = 1.0  # R² differences beyond ±1 are clipped and annotated, so one wide CI cannot hide the rest
    fig, ax = plt.subplots(figsize=(4.6, 0.26 * len(rows) + 0.9))
    for i, (t, s, r) in enumerate(rows[::-1]):
        lo, hi, m = r["d_r2_ci_low"], r["d_r2_ci_high"], r["d_r2_mean"]
        hollow = _flag(r["underpowered"])
        if not (pd.isna(lo) or pd.isna(hi)):
            ax.plot([max(lo, -lim), min(hi, lim)], [i, i], color=SERIES, lw=2, solid_capstyle="round", zorder=2)
            for edge, side in ((lo, -1), (hi, 1)):
                if abs(edge) > lim and edge * side > 0:
                    ax.scatter([side * lim], [i], marker=">" if side > 0 else "<", s=28, color=SERIES, zorder=3)
                    ax.text(side * lim, i + 0.3, f"{edge:.1f}", fontsize=6, color=INK_2,
                            ha="right" if side > 0 else "left")
        mc = max(-lim, min(lim, m))
        if mc != m:  # the mean itself is off-scale: say so, never let the edge stand in for it
            ax.text(mc, i - 0.42, f"mean {m:.1f}", fontsize=6, color=INK_2,
                    ha="left" if mc < 0 else "right")
        ax.scatter([mc], [i], s=30, zorder=4, color="white" if hollow else SERIES,
                   edgecolors=SERIES, linewidths=1.5)
    ax.set_yticks(range(len(rows)), [f"{t}: {LABEL.get(s, s)}" for t, s, _ in rows[::-1]])
    ax.set_xlim(-lim * 1.08, lim * 1.08)
    ax.axvline(0, color=INK_2, lw=0.8)
    ax.axvline(threshold, color=INK_2, lw=0.8, ls=(0, (3, 3)))
    ax.text(threshold, len(rows) - 0.4, f" {threshold:+.2f} (practical threshold)", fontsize=6.5,
            color=INK_2, va="bottom")
    ax.set_xlabel("ΔR² versus size baseline (95 % corrected-t CI)")
    ax.grid(axis="x", color=GRID, lw=0.5)
    fig.text(0, -0.02, f"Learner: {learner}. Hollow marker: underpowered.", fontsize=6.5, color=INK_2)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"fig_delta.{ext}")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("out")
    ap.add_argument("--learner", default="selected")
    ap.add_argument("--threshold", type=float, default=0.05)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    # keep_default_na=False: the baseline set is literally named "null", which pandas would read as NaN
    summary = pd.read_csv(Path(a.results) / "summary.csv", keep_default_na=False, na_values=[""])
    ladder(summary, a.learner, out)
    delta(summary, a.learner, out, a.threshold)
    print(f"figures written to {out}")


if __name__ == "__main__":
    main()
