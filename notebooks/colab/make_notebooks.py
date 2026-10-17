"""Generate the Colab notebooks in this folder from one template, so every notebook shares the same
configuration cell, exact-SHA checkout, pinned install and subprocess execution.

    python notebooks/colab/make_notebooks.py CODE_SHA WEBSITE_SHA

CODE_SHA is the commit whose scripts the notebooks check out (the notebooks themselves are
committed afterwards, in a later commit). All science runs in `python -I` subprocesses calling
docs/research/scripts/*.py, i.e. the same code path as a Linux/HPC run.
"""
import sys
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent
PINS = ('"rdkit==2026.3.6" "chembl_structure_pipeline==1.2.4" "numpy==2.2.6" "pandas==2.3.3" '
        '"scipy==1.15.3" "scikit-learn==1.7.2" "openpyxl==3.1.5" "joblib==1.6.0" "xgboost==3.2.0" '
        '"requests>=2.31" "networkx>=3.2"')

CONFIG = '''# ===== CONFIGURATION (edit here only) =====
import os, pathlib
REPO_URL  = os.environ.get("NB_REPO_URL", "https://github.com/Singati2/TopoQSPR.git")
GIT_SHA   = os.environ.get("NB_GIT_SHA", "{sha}")   # exact commit of the scripts; checked after checkout
MODE      = os.environ.get("NB_MODE", "{mode}")      # {modes}
WORK      = pathlib.Path(os.environ.get("NB_WORK", "/content"))
USE_DRIVE = os.environ.get("NB_USE_DRIVE", "0") == "1"   # persist cache + outputs to your Google Drive
DRIVE_DIR = "/content/drive/MyDrive/topoqspr_repro"      # used only if USE_DRIVE
N_JOBS    = int(os.environ.get("NB_N_JOBS", os.cpu_count() or 1))   # worker processes (Colab free: 2)
{extra}'''

SETUP = '''# ===== SETUP: storage, exact checkout, pinned install (about 2-4 min) =====
import json, subprocess, sys, time
if USE_DRIVE:
    from google.colab import drive
    drive.mount("/content/drive")
    BASE = pathlib.Path(DRIVE_DIR)
else:
    BASE = WORK / "topoqspr_repro"
BASE.mkdir(parents=True, exist_ok=True)
REPO = WORK / "TopoQSPR"
OUT  = BASE / "out"
os.environ["TOPOQSPR_CACHE"] = str(BASE / "cache")   # pinned public downloads + processed records

def sh(cmd, check=True):
    """Run a shell command, stream its output, stop the notebook on failure."""
    print(f"$ {{cmd}}", flush=True)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in p.stdout:
        print(line, end="")
    if p.wait() != 0 and check:
        raise RuntimeError(f"command failed (exit {{p.returncode}}): {{cmd}}")
    return p.returncode

if not (REPO / ".git").exists():
    sh(f"git init -q {{REPO}} && git -C {{REPO}} remote add origin {{REPO_URL}}")
if sh(f"git -C {{REPO}} fetch -q --depth 1 origin {{GIT_SHA}} && git -C {{REPO}} checkout -q FETCH_HEAD", check=False):
    sh(f"git -C {{REPO}} fetch -q origin '+refs/heads/*:refs/remotes/origin/*' && git -C {{REPO}} checkout -q {{GIT_SHA}}")
head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
assert head == GIT_SHA, f"checked out {{head}}, expected {{GIT_SHA}}"
print("checked out", head)
sh(f'{{sys.executable}} -m pip install -q {pins}')
sh(f"{{sys.executable}} -m pip install -q --no-deps -e {{REPO}}")
os.chdir(REPO)
PY = f"{{sys.executable}} -I"   # every analysis runs in a fresh interpreter with the pinned packages
sh(f"""{{PY}} -c "import sys; sys.path.insert(0, 'docs/research/scripts'); import repro_common as RC, json; m = RC.runtime_metadata(); print(json.dumps(m, indent=1)); mm = RC.version_mismatches(); print('VERSION MISMATCHES:', mm or 'none'); sys.exit(1 if mm else 0)" """)
'''

RESTORE = '''# ===== Restore the OPERA/PHYSPROP values (pinned public zip, sha256-checked; values verified by hash) =====
RESTORED = BASE / "e3_restored"
if not (RESTORED / "external/predictions.csv.gz").exists():
    sh(f"{PY} docs/research/scripts/restore_opera_values.py docs/research/data/E3 {RESTORED}")
else:
    print("restored dataset already present:", RESTORED)
'''

HEADER = '''# {title}

[Open in Colab]({colab}) · code commit `{sha}` · branch `repro/colab-first`

{body}

**How to run.** Runtime → Run all, in a fresh CPU runtime. Edit only the configuration cell.
Everything runs as `python -I docs/research/scripts/...` in subprocesses, the same commands as on
Linux/HPC (listed at the end). A `FAIL` line stops the notebook; its numbers must not be used.
**Status:** {status}
'''


def nb(name, title, body, mode, modes, extra, cells, status, restore=True, sha="", pins=PINS):
    n = nbf.v4.new_notebook()
    colab = f"https://colab.research.google.com/github/Singati2/TopoQSPR/blob/repro/colab-first/notebooks/colab/{name}"
    n.cells = [nbf.v4.new_markdown_cell(HEADER.format(title=title, body=body.strip(), sha=sha, colab=colab, status=status)),
               nbf.v4.new_code_cell(CONFIG.format(sha=sha, mode=mode, modes=modes, extra=extra)),
               nbf.v4.new_code_cell(SETUP.format(pins=pins))]
    if restore:
        n.cells.append(nbf.v4.new_code_cell(RESTORE))
    for kind, src in cells:
        n.cells.append(nbf.v4.new_markdown_cell(src) if kind == "md" else nbf.v4.new_code_cell(src))
    n.metadata = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                  "language_info": {"name": "python"}, "colab": {"provenance": []}}
    nbf.write(n, HERE / name)


def main(sha: str, site_sha: str) -> None:
    status = "local clean-environment execution tested (see notebooks/colab/README.md); **Colab unverified**"
    nb("00_environment_and_data_audit.ipynb", "00 · Environment and data audit", """
Checks the runtime, the exact checkout, the pinned versions, the public sources (URLs, licences,
sha256), every published E3 file against `SHA256SUMS`, and the OPERA restore (identifiers, row
counts, hashes of the restored values). Writes `out/00_audit/checks_00_environment_and_data_audit.csv`.
CPU only, about 5 min, about 100 MB download. Nothing is trained.""", "audit", "audit only", "", [
        ("code", 'sh(f"{PY} docs/research/scripts/audit_environment.py docs/research/data/E3 {RESTORED} {OUT}/00_audit")'),
        ("md", "**Linux/HPC equivalent:** `python -I docs/research/scripts/restore_opera_values.py docs/research/data/E3 R && "
               "python -I docs/research/scripts/audit_environment.py docs/research/data/E3 R OUT/00_audit`"),
    ], status, sha=sha)

    nb("01_E3_P0_reproduction.ipynb", "01 · Recompute the E3 and P0 tables from the public artifacts", """
Recomputes, with independent numpy/scipy code, every committed E3 number that can be recomputed
without training: fold R²/MAE, per-cell summaries, all corrected-t paired contrasts, the external
R²/MAE and formula-cluster bootstrap intervals (seed 2026, E3's own), the BH family (seed 99), both
within-formula metrics, and P0. Each value is compared with the committed CSV at |diff| ≤ 1e-9.
The negative control and S-noflt need re-training and are listed, not recomputed.
CPU only, about 3–6 min after the restore. Nothing is trained.""", "reproduce", "reproduce only", "", [
        ("code", 'sh(f"{PY} docs/research/scripts/reproduce_e3.py {RESTORED} {OUT}/01_reproduction")'),
        ("md", "**Linux/HPC equivalent:** `python -I docs/research/scripts/reproduce_e3.py R OUT/01_reproduction`"),
    ], status, sha=sha)

    nb("02_fixed_learner_four_arm.ipynb", "02 · Fixed-learner comparison: C, C+E, C+T, C+E+T", """
Doc 14 (revised, **draft, not approved**). Four representations within each of six fixed learner
families, on the E3 splits, with inner grouped tuning and all preprocessing inside the training
folds. `MODE="smoke"` (default) and `"bench"` train on **permuted labels**: they test the code path
and time it, and cannot answer the question. `MODE="full"` refuses to start unless
`APPROVAL_SHA` holds the commit that froze doc 14 (owner approval).

Long runs: set `USE_DRIVE=True`. Every unit (context, learner, arm, seed, fold) is written
atomically with a sha256 record; re-running the notebook skips only verified units. Split the full
run across sessions or machines with `SHARD="i/N"`; `MAX_UNITS` caps one session.
Smoke: about 10–20 min on 2 CPUs including the feature build. Full: see `budget` below and doc 14 §7.""",
       "smoke", "smoke | bench | full", 'SHARD     = os.environ.get("NB_SHARD", "0/1")\nMAX_UNITS = None          # e.g. 40 to stop after 40 units in this session\n'
       'APPROVAL_SHA = ""        # FULL mode only: commit that froze doc 14 (doc 14 §9)\n', [
        ("md", "### 1. Build the four arms from public sources (cached; fails closed on any hash/split mismatch)"),
        ("code", 'FL = OUT / "fixed_learner"\nif not (FL / "data/SHA256").exists():\n    sh(f"{PY} docs/research/scripts/fixed_learner.py build {FL}")\nelse:\n    print("feature data present; verified by sha256 when loaded")\n'
                 'sh(f"{PY} -c \\"import json; m = json.load(open(\'{FL}/feature_manifest.json\')); [print(c, a, v[\'p\'], \'rank\', v[\'rank_standardised\']) for c in m for a, v in m[c].items()]\\"")'),
        ("md", "### 2. Run (resumable)"),
        ("code", 'if MODE == "full":\n    assert APPROVAL_SHA, "FULL mode needs owner approval (doc 14 §9)"\n    (FL / "APPROVED_FULL_RUN").write_text(APPROVAL_SHA + "\\n")\n'
                 'cap = f"--max-units {MAX_UNITS}" if MAX_UNITS else ""\nsh(f"{PY} docs/research/scripts/fixed_learner.py run {FL} --mode {MODE} --shard {SHARD} --n-jobs {N_JOBS} {cap}")'),
        ("md", "### 3. Analyse whatever is complete (external test labels are read only here)"),
        ("code", 'sh(f"{PY} docs/research/scripts/fixed_learner.py analyse {FL} --mode {MODE} --restored {RESTORED}")\n'
                 'sh(f"{PY} -c \\"import pandas as pd; pd.set_option(\'display.width\', 250); C = pd.read_csv(\'{FL}/results_{MODE}/contrasts.csv\'); print(C.round(4).to_string(index=False))\\"", check=False)'),
        ("code", 'if MODE == "bench":\n    sh(f"{PY} docs/research/scripts/fixed_learner.py budget {FL} --n-jobs {N_JOBS}")'),
        ("md", "**Linux/HPC equivalent:** `python -I docs/research/scripts/fixed_learner.py build OUT; "
               "... run OUT --mode full --shard 0/4 --n-jobs 8; ... analyse OUT --mode full --restored R`"),
    ], status, sha=sha)

    nb("03_external_replication.ipynb", "03 · Patent melting-point set: availability and target-blind census", """
Doc 15. Default: checks the figshare record (CC0 licence, file size, md5) **without downloading**.
`APPROVE_DOWNLOAD=True` (owner approval, doc 15 §7) downloads the 638 MB SDF and runs the
target-blind census: molblocks only (every SD data field, including the melting point, is cut off
before parsing), package standardisation, removal of every constitution present in Bradley or
OPERA, and the frozen go/no-go counts. The confirmation itself is deliberately not implemented: it
is written only after GO and after its protocol is frozen, before any value is read.
Census: about 1–2 h on 2 CPUs; needs about 3 GB disk.""", "metadata", "metadata | census",
       'APPROVE_DOWNLOAD = False   # owner approval to download ~640 MB (doc 15 §7)\n', [
        ("code", 'flag = "--approve-download" if APPROVE_DOWNLOAD else ""\nsh(f"{PY} docs/research/scripts/patent_census.py {OUT}/03_patent {flag} --n-jobs {N_JOBS}")'),
        ("md", "**Linux/HPC equivalent:** `python -I docs/research/scripts/patent_census.py OUT/03_patent [--approve-download]`"),
    ], status, restore=False, sha=sha)

    nb("04_statistical_robustness.ipynb", "04 · Statistical robustness (exploratory)", """
(a) formula-cluster vs molecule-level bootstrap on the external test (training fixed: test-sampling
uncertainty only); (b) coverage of the corrected resampled t under formula-grouped repeated CV, by
simulation on the MP data; (c) training-seed variability of the external pipeline: seed 0 must
reproduce the archived predictions, seeds 1–4 run only in `full` mode (owner approval);
(d) the **pre-declared** within-formula R² next to the **post hoc** centred metric, labelled.
Needs the feature data from notebook 02 step 1 (built here if absent). Smoke: about 10–20 min on 2 CPUs.""",
       "smoke", "smoke | full", "", [
        ("code", 'FL = OUT / "fixed_learner"\nif not (FL / "data/SHA256").exists():\n    sh(f"{PY} docs/research/scripts/fixed_learner.py build {FL}")\n'
                 'sh(f"{PY} docs/research/scripts/robustness.py {RESTORED} {FL} {OUT}/04_robustness --mode {MODE} --n-jobs {N_JOBS}")'),
        ("md", "**Linux/HPC equivalent:** `python -I docs/research/scripts/robustness.py R FL_OUT OUT/04_robustness --mode smoke --n-jobs 8`"),
    ], status, sha=sha)

    nb("05_website_smoke_test.ipynb", "05 · Website branch: smoke test and public-claims check", """
Clones the website code at the exact commit below (branch `fix/website-claims`, **not merged, not
deployed**), runs its composition-baseline smoke test, and checks the public claims in `index.html`
against the code (index count, model count, CV folds, composition baseline, removed Caco-2 column,
forbidden wording). Deploys and merges nothing. CPU only, about 5 min.""", "smoke", "smoke",
       'SITE_SHA = os.environ.get("NB_SITE_SHA", "{site_sha}")   # fix/website-claims commit under test\n'.replace("{site_sha}", site_sha), [
        ("code", 'SITE = WORK / "site"\nif not (SITE / ".git").exists():\n    sh(f"git init -q {SITE} && git -C {SITE} remote add origin {REPO_URL}")\n'
                 'if sh(f"git -C {SITE} fetch -q --depth 1 origin {SITE_SHA} && git -C {SITE} checkout -q FETCH_HEAD", check=False):\n'
                 '    sh(f"git -C {SITE} fetch -q origin \'+refs/heads/*:refs/remotes/origin/*\' && git -C {SITE} checkout -q {SITE_SHA}")\n'
                 'sh(f"{sys.executable} -m pip install -q shap==0.49.1 chembl_webresource_client==0.10.9 rapidfuzz==3.14.5")\n'
                 'sh(f"cd {SITE} && {PY} colab/smoke_website_baseline.py 2>&1 | grep -v -i warn")\n'
                 'sh(f"{PY} docs/research/scripts/website_claims_check.py {SITE} {OUT}/05_website")'),
        ("md", "**Linux/HPC equivalent:** `git worktree add SITE <SITE_SHA>; (cd SITE && python colab/smoke_website_baseline.py); "
               "python -I docs/research/scripts/website_claims_check.py SITE OUT/05_website`"),
    ], status, restore=False, sha=sha)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
