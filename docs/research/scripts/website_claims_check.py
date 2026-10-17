"""Notebook 05: check the website's public-facing claims against its implementation.

Run against a checkout of the website code (branch fix/website-claims, exact SHA in the notebook).
Reads index.html / README.md and imports config.py and ml.py; trains nothing except via the
branch's own smoke test (colab/smoke_website_baseline.py), which the notebook runs separately.
Writes a PASS/FAIL report. Deploys and merges nothing.

    python -I docs/research/scripts/website_claims_check.py SITE_DIR OUT_DIR
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import repro_common as RC  # noqa: E402

FORBIDDEN = [r"pre-?registered", r"deployed v1\.0", r"state[- ]of[- ]the[- ]art", r"clinically validated",
             r"independent(ly)? validat"]


def main(site: Path, out: Path) -> RC.Checks:
    sys.path.insert(0, str(site))
    import config  # noqa: PLC0415  (website module)
    import ml      # noqa: PLC0415
    ck = RC.Checks("05_website_claims")
    html = (site / "index.html").read_text(errors="replace")
    readme = (site / "README.md").read_text(errors="replace")
    n_idx = len(config.TOPO_INDICES)
    claims = sorted(set(int(x) for x in re.findall(r"(\d+) (?:graph-theoretic )?topological indices", html)))
    claims += sorted(set(int(x) for x in re.findall(r"Computing (\d+) indices", html)))
    # README line 6 describes the separate TopoQSPR package (104 verified indices); the site is index.html
    ck.require("every 'N topological indices' claim on the site equals len(config.TOPO_INDICES)",
               bool(claims) and all(c == n_idx for c in claims), claims, n_idx)
    n_models = len(ml.MODELS_AND_GRIDS)
    mc = sorted(set(int(x) for x in re.findall(r"(\d+) (?:ML )?models", html)))
    ck.require("every 'N ML models' claim equals len(ml.MODELS_AND_GRIDS)", bool(mc) and all(c == n_models for c in mc),
               mc, n_models)
    cvc = sorted(set(int(x) for x in re.findall(r"(\d+)-fold CV", html)))
    ck.require("'K-fold CV' claim equals config.CV_FOLDS", bool(cvc) and all(c == config.CV_FOLDS for c in cvc),
               cvc, config.CV_FOLDS)
    dc = sorted(set(int(x) for x in re.findall(r"up to (\d+) drugs", html)))
    ck.require("'up to N drugs' claim equals config.MAX_DRUGS", all(c == config.MAX_DRUGS for c in dc), dc, config.MAX_DRUGS)
    src = (site / "ml.py").read_text()
    ck.require("composition-only baseline implemented (ml._composition_matrix, Comp_R2_mean, Delta_R2)",
               "_composition_matrix" in src and "Comp_R2_mean" in src and "Delta_R2" in src, "", "")
    ck.require("results table shows the composition baseline and dR2 columns",
               "Composition-only R²" in html and "ΔR²" in html, "", "")
    fet = (site / "fetcher.py").read_text()
    ck.require("no Caco-2 prediction column produced by fetcher.py",
               not re.search(r"['\"]ADMET_Caco2", fet), re.findall(r".{20}Caco2.{20}", fet)[:2], "none")
    for pat in FORBIDDEN:
        hits = re.findall(rf".{{0,40}}{pat}.{{0,40}}", html + readme, flags=re.I)
        ck.require(f"forbidden public wording absent: /{pat}/", not hits, hits[:2], "none")
    ck.info("site git SHA", RC.git_sha(site))
    ck.finish(out)
    return ck


if __name__ == "__main__":
    main(Path(sys.argv[1]).expanduser(), Path(sys.argv[2]).expanduser())
