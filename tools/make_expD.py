"""Experiment D (replication of published studies) table, figure and macros.

Usage: python tools/make_expD.py OUT_DIR RUN_DIR [RUN_DIR ...]

Each RUN_DIR is one `run_experiment.py exp-d` output (summary.json, reported_fits.csv, fits.csv,
nested.csv, controls.csv, index_check.csv).  Writes tab_expD.tex, fig_expD.pdf/.png and
expD_macros.tex.  Only *compared* reported fits (a form the harness reproduces, first
occurrence of a duplicate) appear in the table; the number of listed-only fits is a macro.
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
STAGES = ["published", "in-sample", "LOO $Q^2$", "nested CV", "size baseline"]
CITE = {"almulla2026_antifungal": ("Almulla et al.\\ 2026", "almulla2026antifungal"),
        "kabeer2026_antidepressants": ("Kabeer et al.\\ 2026", "kabeer2026antidepressants"),
        "godlin2026_mdd": ("Godlin \\& Radha 2026", "godlin2026mdd")}
FORM = {"linear": "lin.", "quadratic": "quad.", "cubic": "cubic"}

plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK_2, "axes.labelcolor": INK, "xtick.color": INK_2,
                     "ytick.color": INK_2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, keep_default_na=False, na_values=[""])


def _f(x, nd=2) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "--"
    return "--" if not math.isfinite(v) else f"{v:.{nd}f}"


def _p(x) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "--"
    return "--" if not math.isfinite(v) else ("$<$0.01" if v < 0.01 else f"{v:.2f}")


def _tex(s) -> str:
    return str(s).replace("_", "\\_").replace("%", "\\%").replace("&", "\\&")


def load(run: Path) -> dict:
    s = json.loads((run / "summary.json").read_text())
    fits = _read(run / "fits.csv")
    nested = _read(run / "nested.csv")
    ctrl = _read(run / "controls.csv")
    chk = _read(run / "index_check.csv") if (run / "index_check.csv").exists() else pd.DataFrame()
    rep = _read(run / "reported_fits.csv")
    report, loo_pub = {}, {}
    info = json.loads((run / "run_info.json").read_text()) if (run / "run_info.json").exists() else {}
    spec_path = Path((info.get("spec") or {}).get("path", ""))
    for cand in (spec_path, Path(__file__).resolve().parents[1] / spec_path):
        rp = cand.with_suffix(".report.json")
        if rp.exists():
            report = json.loads(rp.read_text())
            break
    for r in report.get("reported_loocv", []):
        if r.get("q2_loo") is not None:
            loo_pub[(r["property"], r["index"], str(r["form"]).lower())] = float(r["q2_loo"])
    rows = []
    for _, r in fits[fits["reported_r2"].notna()].iterrows():
        p = r["property"]
        g = nested[(nested["property"] == p) & (nested["learner"] == "genre")]
        ti = g[(g["set"] == "TI_study") & (g["split"] == "random")]
        sz = g[(g["set"] == "size") & (g["split"] == "random")]
        sc = g[(g["set"] == "TI_study") & (g["split"] == "scaffold")]
        c = ctrl[ctrl["property"] == p]
        one = lambda df, k: float(df.iloc[0][k]) if len(df) and pd.notna(df.iloc[0][k]) else math.nan  # noqa: E731
        rows.append(dict(
            study=s["study"], property=p, index=r["index"], form=r["form"], n=int(r["n"]),
            published=float(r["reported_r2"]), insample=float(r["r2"]), loo=float(r["loo_q2"]),
            loo_published=loo_pub.get((p, r["index"], str(r["form"]).lower()), math.nan),
            nested_ti=one(ti, "r2_pooled_mean"), nested_size=one(sz, "r2_pooled_mean"),
            d_r2=one(ti, "d_r2_vs_size"), d_lo=one(ti, "d_r2_ci_low"), d_hi=one(ti, "d_r2_ci_high"),
            scaffold=one(sc, "r2_pooled_mean"), registry_p=one(c, "registry_p"), perm_p=one(c, "perm_p")))
    # printed index values: compounds whose recomputable printed values all agree
    agree_comp = checked_comp = 0
    if len(chk):
        rec = chk[chk["recomputed"].notna()]
        by = rec.groupby("compound")["agree"].apply(lambda x: all(str(v).lower() == "true" for v in x))
        checked_comp, agree_comp = int(len(by)), int(by.sum())
    return dict(study=s["study"], n_compounds=int(s["n_compounds_used"]), n_spec=int(s["n_compounds_spec"]),
                n_reported=int(len(rep)), n_compared=int((rep["compared"].astype(str).str.lower() == "true").sum()),
                n_index_labels=len(report.get("index_map", {})) + len(report.get("unmapped_indices", []))
                if report else int(s["n_index_labels"]),
                n_mapped=len(report.get("index_map", {})) if report else int(s["n_mapped"]),
                n_loocv_published=len(report.get("reported_loocv", [])),
                values_checked=int(len(chk[chk["recomputed"].notna()])) if len(chk) else 0,
                values_agree=int((chk["agree"].astype(str).str.lower() == "true").sum()) if len(chk) else 0,
                compounds_checked=checked_comp, compounds_agree=agree_comp, fits=rows)


def table(studies: list[dict], out: Path) -> None:
    lines = []
    for st in studies:
        name, key = CITE.get(st["study"], (_tex(st["study"]), ""))
        cite = f"~\\cite{{{key}}}" if key else ""
        lines.append(
            f"\\multicolumn{{12}}{{l}}{{\\textit{{{name}{cite}: $n={st['n_compounds']}$, "
            f"{st['n_mapped']}/{st['n_index_labels']} indices in the registry, {st['n_compared']}/{st['n_reported']} "
            f"reported fits compared; printed index values agree for {st['compounds_agree']}/{st['compounds_checked']} "
            f"compounds}}}} \\\\")
        for r in st["fits"]:
            lines.append(" & ".join([
                _tex(r["property"]), "$" + _tex(r["index"]).replace("S\\_", "S_{") + ("}$" if r["index"].startswith("S_") else "$"),
                FORM.get(r["form"], r["form"]), _f(r["published"], 3), _f(r["insample"], 3),
                _f(r["loo_published"], 2), _f(r["loo"], 2),
                _f(r["nested_ti"], 2), _f(r["nested_size"], 2),
                f"{_f(r['d_r2'], 2)} [{_f(r['d_lo'], 2)}, {_f(r['d_hi'], 2)}]", _f(r["scaffold"], 2),
                _p(r["registry_p"]), _p(r["perm_p"])]) + " \\\\")
        lines.append("\\addlinespace")
    tex = r"""\begin{table}[t]
\centering\scriptsize
\caption{Experiment D, replication of published single-index studies on their own compounds and property values. Published: the $R^2$ printed in the paper; in-sample: the same fit recomputed by \ml{} on PubChem structures; publ.\ LOO: the leave-one-out $Q^2$ printed in the paper, where given; LOO: leave-one-out $Q^2$ of the recomputed fit; nested CV: pooled $R^2$ of the genre's own procedure (index and form chosen by inner CV among the study's indices) under $5\times5$ repeated nested CV, and the same procedure on the size baseline; $\Delta R^2$: TI minus size on identical folds with the corrected repeated $k$-fold $t$ interval; scaffold: nested CV with Bemis--Murcko scaffold groups; $p_{\mathrm{reg}}$: share of random registry indices whose best in-sample fit reaches the study's best $R^2$; $p_{\mathrm{perm}}$: Y-randomisation. All studies have $n<50$ and are underpowered by the plan's rule. Generated by \texttt{tools/make\_expD.py}.}
\label{tab:expD}
\resizebox{\linewidth}{!}{%%
\begin{tabular}{lllrrrrrrlrll}
\toprule
Property & Index & Form & Published & In-sample & Publ.\ LOO & LOO $Q^2$ & Nested TI & Nested size & $\Delta R^2$ [95\%% CI] & Scaffold & $p_{\mathrm{reg}}$ & $p_{\mathrm{perm}}$ \\
\midrule
%s
\bottomrule
\end{tabular}}
\end{table}
""" % "\n".join(lines)
    (out / "tab_expD.tex").write_text(tex)


def figure(studies: list[dict], out: Path) -> None:
    fig, axes = plt.subplots(1, len(studies), figsize=(2.3 * len(studies) + 0.6, 2.6), sharey=True,
                             gridspec_kw=dict(wspace=0.12))
    axes = list(axes) if len(studies) > 1 else [axes]
    for ax, st in zip(axes, studies):
        for r in st["fits"]:
            ys = [r["published"], r["insample"], r["loo"], r["nested_ti"], r["nested_size"]]
            ax.plot(range(5), ys, color=SERIES, lw=1, alpha=0.55, zorder=2)
            ax.scatter(range(4), ys[:4], s=12, color=SERIES, zorder=3)
            ax.scatter([4], [ys[4]], s=14, facecolors="white", edgecolors=INK_2, linewidths=1, zorder=3)
        ax.set_title(CITE.get(st["study"], (st["study"], ""))[0].replace("\\ ", " ").replace("\\&", "&")
                     + f" (n={st['n_compounds']})", fontsize=7.5, color=INK)
        ax.set_xticks(range(5), STAGES, rotation=35, ha="right", fontsize=6.5)
        ax.axhline(0, color=INK_2, lw=0.6)
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_xlim(-0.3, 4.3)
    lo = min(min(r["published"], r["insample"], r["loo"], r["nested_ti"], r["nested_size"])
             for st in studies for r in st["fits"] if all(map(math.isfinite, [r["loo"], r["nested_ti"], r["nested_size"]])))
    axes[0].set_ylim(min(-0.05, lo - 0.05), 1.02)
    axes[0].set_ylabel("$R^2$ (or $Q^2$)")
    fig.text(0.01, -0.12, "One line per reported fit (property × index × form). Filled: the study's index; "
             "hollow: the size baseline under the same nested CV.", fontsize=6.5, color=INK_2, va="top")
    for ext in ("pdf", "png"):
        fig.savefig(out / f"fig_expD.{ext}")
    plt.close(fig)


def macros(studies: list[dict], out: Path) -> None:
    fits = [r for st in studies for r in st["fits"]]
    diffs = [r["insample"] - r["published"] for r in fits]
    loo_drop = [r["published"] - r["loo"] for r in fits if math.isfinite(r["loo"])]
    ti_le_size = sum(1 for r in fits if math.isfinite(r["nested_ti"]) and r["nested_ti"] <= r["nested_size"])
    d_sig = sum(1 for r in fits if math.isfinite(r["d_lo"]) and r["d_lo"] > 0)
    reg_sig = sum(1 for r in fits if math.isfinite(r["registry_p"]) and r["registry_p"] < 0.05)
    props = {(r["study"], r["property"]) for r in fits}
    m = {
        "expDnStudies": len(studies),
        "expDnFitsReported": sum(st["n_reported"] for st in studies),
        "expDnFitsCompared": len(fits),
        "expDnLoocvPublished": sum(st["n_loocv_published"] for st in studies),
        "expDnLoocvCompared": sum(1 for r in fits if math.isfinite(r["loo_published"])),
        "expDmaxLooPublishedDiff": f"{max((abs(r['loo'] - r['loo_published']) for r in fits if math.isfinite(r['loo_published'])), default=math.nan):.2f}",
        "expDnProperties": len(props),
        "expDminDiff": f"{min(diffs):+.3f}", "expDmaxDiff": f"{max(diffs):+.3f}",
        "expDnInsampleAbove": sum(1 for d in diffs if d > 0),
        "expDnExactFits": sum(1 for d in diffs if abs(d) <= 0.0015),
        "expDnCompoundsAllIndicesDiffer": sum(st["compounds_checked"] - st["compounds_agree"] for st in studies),
        "expDmaxLooDrop": f"{max(loo_drop):.2f}",
        "expDnLooDropOverTenth": sum(1 for d in loo_drop if d > 0.1),
        "expDnTIleSize": ti_le_size,
        "expDnDeltaSig": d_sig,
        "expDnRegistrySig": reg_sig,
        "expDcompoundsAgree": sum(st["compounds_agree"] for st in studies),
        "expDcompoundsChecked": sum(st["compounds_checked"] for st in studies),
        "expDvaluesAgree": sum(st["values_agree"] for st in studies),
        "expDvaluesChecked": sum(st["values_checked"] for st in studies),
    }
    (out / "expD_macros.tex").write_text("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in m.items()))
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    studies = [load(Path(r)) for r in a.runs]
    table(studies, out)
    figure(studies, out)
    m = macros(studies, out)
    for st in studies:
        print(f"{st['study']:28s} n={st['n_compounds']:2d} fits {st['n_compared']}/{st['n_reported']} "
              f"index values {st['values_agree']}/{st['values_checked']} (compounds {st['compounds_agree']}/{st['compounds_checked']})")
    print({k: v for k, v in m.items() if k.startswith("expDn") or "Diff" in k or "Drop" in k})


if __name__ == "__main__":
    main()
