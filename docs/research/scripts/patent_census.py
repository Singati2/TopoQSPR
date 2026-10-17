"""Notebook 03 / doc 15 step 1: target-blind census of the patent melting-point set.

Source: Tetko IV, Lowe DM, Williams AJ. J Cheminform 8:2 (2016), doi:10.1186/s13321-016-0113-y;
data figshare doi:10.6084/m9.figshare.2007426.v1 (CC0), file 3575841 `melting_point.sdf`
(637,981,904 bytes, md5 20618014e0fe69ce7dbd1de87413e1f1 per the figshare API, 2026-10-17).

Three stages, each gated:
  metadata  (default, no download)  figshare API: licence, file size, md5 still as pinned.
  census    (--approve-download)    download the SDF (md5 + sha256 recorded), read the MOLBLOCKS
            ONLY (every SD data field, including the melting point, is cut off before RDKit sees
            the record), standardise with the package pipeline, remove every parent whose
            constitution (InChIKey first block) occurs in Bradley (all rows, incl. donotuse) or
            OPERA MP/logP, and count formula groups with >= 5 constitutions. GO/NO-GO on the
            counts frozen in doc 15 §3.
  confirmation  not implemented here on purpose: it may be written only after GO and after the
            owner freezes doc 15 §5, before any melting-point value is read.

    python -I docs/research/scripts/patent_census.py OUT                       # metadata only
    python -I docs/research/scripts/patent_census.py OUT --approve-download    # census (~640 MB)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import repro_common as RC  # noqa: E402

ARTICLE = 2007426
FILE_ID = 3575841
FILE_NAME = "melting_point.sdf"
PINNED_MD5 = "20618014e0fe69ce7dbd1de87413e1f1"
PINNED_SIZE = 637981904
GO = dict(min_molecules_in_ge5_groups=3400, min_ge5_formula_groups=300)  # doc 15 §3, frozen there


def metadata(ck: RC.Checks) -> dict:
    import requests
    d = requests.get(f"https://api.figshare.com/v2/articles/{ARTICLE}", timeout=60).json()
    f = next((x for x in d["files"] if x["id"] == FILE_ID), None)
    ck.require("figshare licence is CC0", d["license"]["name"] == "CC0", d["license"]["name"], "CC0")
    ck.require("file present", f is not None, "", FILE_NAME)
    ck.require("file md5 as pinned", f is not None and f.get("computed_md5") == PINNED_MD5,
               f and f.get("computed_md5"), PINNED_MD5)
    ck.require("file size as pinned", f is not None and f["size"] == PINNED_SIZE, f and f["size"], PINNED_SIZE)
    return dict(title=d["title"], doi=d.get("doi"), licence=d["license"], file=f)


def molblocks(path: Path):
    """Yield (record_index, title, molblock) with every SD data field removed (target-blind)."""
    buf, i = [], 0
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("$$$$"):
                text = "".join(buf)
                cut = text.find("M  END")
                block = text[: cut + len("M  END")] + "\n" if cut >= 0 else ""
                yield i, (text.splitlines() or [""])[0].strip(), block
                buf, i = [], i + 1
            else:
                buf.append(line)


def census(out: Path, ck: RC.Checks, n_jobs: int | None, limit: int | None = None) -> dict:
    import requests
    from rdkit import Chem, RDLogger
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    from topoqspr import expdata
    RDLogger.DisableLog("rdApp.*")
    raw = out / "raw"; raw.mkdir(parents=True, exist_ok=True)
    sdf = raw / FILE_NAME
    if not sdf.exists():
        tmp = sdf.with_suffix(".part")
        with requests.get(f"https://ndownloader.figshare.com/files/{FILE_ID}", stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(tmp, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
        tmp.replace(sdf)
    md5 = hashlib.md5(sdf.read_bytes()).hexdigest()
    ck.require("downloaded md5 as pinned", md5 == PINNED_MD5, md5, PINNED_MD5)
    sha = RC.sha256_file(sdf)
    ck.info("downloaded sha256 (record this)", sha)
    smiles, n_parse_fail = {}, 0
    for i, title, block in molblocks(sdf):
        if limit is not None and i >= limit:
            break
        m = Chem.MolFromMolBlock(block) if block else None
        if m is None:
            n_parse_fail += 1
            continue
        smiles[i] = Chem.MolToSmiles(m)
    std = expdata._std_many(smiles.values(), n_jobs=n_jobs)
    rows = []
    for i, s in smiles.items():
        code, detail, ik, can, note, ncomp = std[s]
        rows.append(dict(record=i, ok=code == "", reject=code, inchikey=ik, smiles=can, n_input_components=ncomp))
    R = pd.DataFrame(rows)
    ok = R[R.ok].drop_duplicates("inchikey").copy()
    ok["formula"] = [CalcMolFormula(Chem.MolFromSmiles(s)) for s in ok.smiles]
    ok["constitution"] = ok.inchikey.str[:14]
    used = set()
    for name in ("mp_bradley", "mp_opera", "logp_opera"):
        recs = expdata.load_records(name)  # all rows, every status (Bradley donotuse included)
        used |= set(recs.inchikey.dropna().astype(str).str[:14])
    fresh = ok[~ok.constitution.isin(used)]
    ncons = fresh.drop_duplicates("constitution").groupby("formula").size()
    ge5 = fresh[fresh.formula.map(ncons) >= 5]
    res = dict(records_read=len(smiles) + n_parse_fail, molblock_parse_failures=n_parse_fail,
               standardised_ok=int(R.ok.sum()), rejected_by_reason=R[~R.ok].reject.value_counts().to_dict(),
               multi_component_inputs=int((R.n_input_components > 1).sum()),
               unique_parents=len(ok), parents_overlapping_bradley_or_opera=int(len(ok) - len(fresh)),
               parents_after_overlap_removal=len(fresh),
               ge5_formula_groups=int((ncons >= 5).sum()), molecules_in_ge5_groups=len(ge5),
               sdf_sha256=sha, limit=limit, thresholds=GO)
    res["GO"] = bool(res["molecules_in_ge5_groups"] >= GO["min_molecules_in_ge5_groups"]
                     and res["ge5_formula_groups"] >= GO["min_ge5_formula_groups"])
    if limit is not None:
        res["GO"] = None
        res["NOTE"] = "limited run (smoke): counts are not a census, no go/no-go"
    RC.atomic_write_text(json.dumps(res, indent=1, default=str), out / "census.json")
    ge5[["inchikey", "formula", "constitution"]].to_csv(out / "census_ge5_parents.csv", index=False)
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--approve-download", action="store_true",
                    help="owner approval to download ~640 MB and run the target-blind census")
    ap.add_argument("--limit", type=int, help="smoke: only the first N records (no go/no-go)")
    ap.add_argument("--n-jobs", type=int)
    a = ap.parse_args()
    out = Path(a.out).expanduser(); out.mkdir(parents=True, exist_ok=True)
    ck = RC.Checks("03_patent_census")
    md = metadata(ck)
    RC.atomic_write_text(json.dumps(md, indent=1), out / "figshare_metadata.json")
    if a.approve_download:
        print(json.dumps(census(out, ck, a.n_jobs, a.limit), indent=1, default=str))
    else:
        ck.info("census", "not run: needs --approve-download (owner approval, doc 15 §7)")
    ck.finish(out)
