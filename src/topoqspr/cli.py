"""Command-line interface (``topoqspr``; also ``python -m topoqspr``).

Subcommands::

    topoqspr headings TEXT                       suggest MeSH IDs for free text
    topoqspr cohort --mesh D008881 --out DIR     controlled cohort -> DIR/manifest.csv + .json
    topoqspr indices --smiles-csv F --out G      verified indices for a SMILES file
    topoqspr run --config study.json --out DIR   full study (docs/PIPELINE.md)
    topoqspr rerun RUN_RECORD --out DIR          re-execute and compare byte for byte

Time stamps written into outputs default to the current UTC time and can be
fixed with ``--timestamp``. Exit status: 0 success, 1 a rerun whose result
CSVs differ, 2 usage or input errors.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import logging
import sys
from pathlib import Path
from typing import Sequence

from . import __version__

__all__ = ["main", "build_parser"]


def _now_utc() -> str:
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z")


def _client(args):
    from .chembl import ChemblClient

    return ChemblClient(args.cache_dir, offline=args.offline, release=args.release)


def _add_chembl_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--cache-dir", default=None,
                   help="ChEMBL cache root (default ~/.cache/topoqspr/chembl)")
    p.add_argument("--offline", action="store_true", help="serve from the cache only")
    p.add_argument("--release", default=None,
                   help="pin the ChEMBL release (e.g. ChEMBL_37); online: fail if the API "
                        "serves another one")


def cmd_headings(args) -> int:
    from .cohort import suggest_headings

    df = suggest_headings(args.text, _client(args), limit=args.limit)
    df.to_csv(sys.stdout, index=False, lineterminator="\n")
    return 0


def cmd_cohort(args) -> int:
    from .cohort import build_cohort

    client = _client(args)
    cohort = build_cohort(args.mesh, client, min_phase=None if args.all_phases else args.min_phase,
                          exclude_sole_evidence=args.exclude_sole_evidence or None,
                          evidence_sources=args.evidence_sources or None,
                          exclude_diagnostics=args.exclude_diagnostics,
                          disease_label=args.label)
    out = Path(args.out)
    csv_path, json_path = cohort.write(out / "manifest.csv", timestamp_utc=args.timestamp or _now_utc())
    print(f"{cohort.disease_label} ({';'.join(cohort.parameters['mesh_ids'])}, {cohort.release}): "
          f"{len(cohort.included)} included of {len(cohort.manifest)} manifest rows -> {csv_path}")
    return 0


def cmd_indices(args) -> int:
    import pandas as pd
    from rdkit import rdBase

    from .graph import standardize
    from .pipeline import compute_index_table

    df = pd.read_csv(args.smiles_csv, dtype=str, keep_default_na=False, na_values=[])
    for c in (args.id_column, args.smiles_column):
        if c not in df.columns:
            print(f"error: column {c!r} not in {args.smiles_csv}", file=sys.stderr)
            return 2
    ids = df[args.id_column].tolist()
    if len(set(ids)) != len(ids):
        print(f"error: duplicate ids in {args.smiles_csv}", file=sys.stderr)
        return 2
    rows, graphs, ok_ids = [], [], []
    with rdBase.BlockLogs():
        for rid, smi in zip(ids, df[args.smiles_column]):
            res = standardize(smi)
            rows.append(dict(id=rid, input_smiles=smi, status="ok" if res.ok else "rejected",
                             reason="" if res.ok else res.rejection.code,
                             smiles=res.canonical_smiles or "", inchikey=res.inchikey or ""))
            if res.ok:
                graphs.append(res.graph)
                ok_ids.append(rid)
    meta = pd.DataFrame(rows).set_index("id")
    vals, failures = compute_index_table(graphs, ok_ids, exact=True)
    # object columns keep exact integers (Hosoya Z > 2**53) next to rejected (empty) rows
    out = meta.join(vals.astype(object).reindex(meta.index)).reset_index()
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False, lineterminator="\n")
    fail = pd.DataFrame(
        [dict(row=r["id"], feature="*", reason=f"rejected: {r['reason']}")
         for r in rows if r["status"] != "ok"] + failures,
        columns=["row", "feature", "reason"])
    fail_path = out_path.with_name(out_path.stem + ".failures.csv")
    fail.to_csv(fail_path, index=False, lineterminator="\n")
    print(f"{len(ok_ids)} of {len(ids)} molecules standardised; {vals.shape[1]} indices; "
          f"{len(failures)} undefined values -> {out_path} (reasons: {fail_path})")
    return 0


def cmd_run(args) -> int:
    from .pipeline import load_config, run_study

    cfg_path = Path(args.config)
    cfg = load_config(cfg_path)
    if args.n_jobs is not None or args.cache_dir is not None or args.offline:
        rt = dict(cfg.get("runtime") or {})
        if args.n_jobs is not None:
            rt["n_jobs"] = args.n_jobs
        if args.cache_dir is not None:
            rt["chembl_cache_dir"] = args.cache_dir
        if args.offline:
            rt["offline"] = True
        cfg["runtime"] = rt
    res = run_study(cfg, args.out, timestamp_utc=args.timestamp or _now_utc(),
                    base_dir=cfg_path.resolve().parent, overwrite=args.overwrite)
    rec = res.run_record
    print(f"{rec['counts']['standardised']} molecules; results in {res.out_dir}; "
          f"config sha256 {rec['config_sha256'][:12]}; "
          f"evaluate {rec['timings_s'].get('evaluate', 0):.1f} s")
    return 0


def cmd_rerun(args) -> int:
    from .pipeline import rerun

    rep = rerun(args.run_record, args.out, timestamp_utc=args.timestamp or _now_utc(),
                cache_dir=args.cache_dir, base_dir=args.base_dir, n_jobs=args.n_jobs,
                offline=not args.online)
    print(json.dumps({k: rep[k] for k in ("ok", "identical", "different", "missing")}, indent=1))
    if rep["package_differences"]:
        print(f"note: {len(rep['package_differences'])} package versions differ from the record",
              file=sys.stderr)
    return 0 if rep["ok"] else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="topoqspr", description=__doc__.split("\n\n")[0])
    p.add_argument("--version", action="version", version=f"topoqspr {__version__}")
    p.add_argument("-v", "--verbose", action="store_true", help="log progress to stderr")
    sub = p.add_subparsers(dest="command", required=True)

    h = sub.add_parser("headings", help="suggest MeSH headings (IDs) for free text")
    h.add_argument("text")
    h.add_argument("--limit", type=int, default=10)
    _add_chembl_args(h)
    h.set_defaults(func=cmd_headings)

    c = sub.add_parser("cohort", help="build a controlled cohort from MeSH descriptor IDs")
    c.add_argument("--mesh", action="append", required=True, metavar="D######",
                   help="confirmed MeSH descriptor ID (repeatable)")
    c.add_argument("--out", required=True, help="output directory (manifest.csv + manifest.json)")
    c.add_argument("--min-phase", type=float, default=4.0)
    c.add_argument("--all-phases", action="store_true", help="disable the phase rule")
    c.add_argument("--evidence-sources", action="append", default=None, metavar="REF_TYPE")
    c.add_argument("--exclude-sole-evidence", action="append", default=None, metavar="REF_TYPE")
    c.add_argument("--exclude-diagnostics", action="store_true")
    c.add_argument("--label", default=None, help="disease label (default: the MeSH headings)")
    c.add_argument("--timestamp", default=None, help="UTC time stamp for the sidecar")
    _add_chembl_args(c)
    c.set_defaults(func=cmd_cohort)

    i = sub.add_parser("indices", help="compute the verified indices for a SMILES CSV")
    i.add_argument("--smiles-csv", required=True)
    i.add_argument("--out", required=True, help="output CSV (reasons in <out>.failures.csv)")
    i.add_argument("--id-column", default="id")
    i.add_argument("--smiles-column", default="smiles")
    i.set_defaults(func=cmd_indices)

    r = sub.add_parser("run", help="run a study configuration")
    r.add_argument("--config", required=True, help="study config (.json; .yaml if PyYAML is installed)")
    r.add_argument("--out", required=True, help="results directory (must be empty)")
    r.add_argument("--timestamp", default=None)
    r.add_argument("--n-jobs", type=int, default=None, help="outer folds in parallel")
    r.add_argument("--cache-dir", default=None, help="ChEMBL cache root")
    r.add_argument("--offline", action="store_true")
    r.add_argument("--overwrite", action="store_true", help="replace a non-empty --out")
    r.set_defaults(func=cmd_run)

    rr = sub.add_parser("rerun", help="re-execute a run record and compare result CSVs")
    rr.add_argument("run_record")
    rr.add_argument("--out", required=True)
    rr.add_argument("--timestamp", default=None)
    rr.add_argument("--cache-dir", default=None, help="ChEMBL cache root holding the release")
    rr.add_argument("--base-dir", default=None, help="directory for relative input paths")
    rr.add_argument("--n-jobs", type=int, default=None)
    rr.add_argument("--online", action="store_true",
                    help="allow network access (default: cache only)")
    rr.set_defaults(func=cmd_rerun)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(asctime)s %(name)s %(message)s", stream=sys.stderr)
    try:
        return int(args.func(args) or 0)
    except (ValueError, KeyError, FileNotFoundError, FileExistsError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
