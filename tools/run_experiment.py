#!/usr/bin/env python
"""Run the experiments of EXPERIMENT_PLAN §8 (A, C, D) and the §7 size-dependence analysis.

    PYTHONPATH=src python tools/run_experiment.py exp-a --mesh D008881 --naive-text migraine --out DIR
    PYTHONPATH=src python tools/run_experiment.py exp-c --mesh D008881 --target MolLogP --out DIR
    PYTHONPATH=src python tools/run_experiment.py exp-d --spec study.json --out DIR
    PYTHONPATH=src python tools/run_experiment.py size-dependence --mesh D008881 --out DIR

Inputs for exp-c and size-dependence: ``--mesh`` (controlled ChEMBL cohort), ``--manifest``
(a cohort manifest CSV written by ``topoqspr cohort``) or ``--smiles-csv`` (``id``, ``smiles``).
Targets: ``--target NAME`` from a target provider (default ``rdkit``) or ``--target-csv F
--target-column C`` (joined on the molecule id; missing values are dropped with a count).

**Quarantine.** Until docs/EXPERIMENT_PLAN.md is frozen, every modelling step on a ChEMBL
disease cohort (exp-a targets/correlations/``--model``, exp-c, size-dependence) is refused for
any MeSH ID other than the pilot D008881 (``experiments.check_modelling_allowed``), before
any network access. A ``--manifest`` whose included rows carry no MeSH ID (an ATC / pooled
approved-drug manifest, or no ``mesh_ids`` column) is refused too. exp-a on a non-pilot
disease runs only with ``--counts-only`` (cohort sizes, overlaps, size distributions).
``--plan-frozen`` lifts the guard; pass it only after the plan's header records the frozen
commit.

Outputs (``--out``, must be empty unless ``--overwrite``): the experiment's CSV tables,
``summary.json`` and ``run_info.json`` (command line, time stamp, git state, package
version, inputs and their SHA-256, ChEMBL release, wall time). Files are written to a
staging directory and moved into ``--out`` only when the run succeeded. Docs:
docs/EXPERIMENTS.md. Exit status: 0 success, 2 input/configuration error or quarantine
refusal, 3 unexpected internal error.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import logging
import shutil
import sys
import time
import traceback
from pathlib import Path
from typing import Sequence

EXIT_OK, EXIT_INPUT, EXIT_INTERNAL = 0, 2, 3
log = logging.getLogger("run_experiment")


def _now_utc() -> str:
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z")


# ------------------------------------------------------------------- inputs
def _client(args):
    from topoqspr.chembl import ChemblClient

    return ChemblClient(args.cache_dir, offline=args.offline, release=args.release)


def _read_csv(path):
    import pandas as pd

    return pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])


def _manifest_mesh_ids(df) -> list[str]:
    """MeSH IDs of the included rows of a manifest ([] if none or no column)."""
    from topoqspr.experiments import cohort_mesh_ids

    return cohort_mesh_ids(df)


def _load_structures(args, info: dict):
    """(Structures, mesh_ids or None) from --mesh / --manifest / --smiles-csv."""
    from topoqspr import experiments as ex
    from topoqspr import provenance as pv

    if args.mesh:
        ex.check_modelling_allowed(args.mesh, args.plan_frozen)
        from topoqspr.cohort import build_cohort

        client = _client(args)
        coh = build_cohort(args.mesh, client)
        inc = coh.included
        info["cohort"] = dict(source="chembl", mesh_ids=coh.parameters["mesh_ids"],
                              release=coh.release, release_date=coh.release_date,
                              disease_label=coh.disease_label,
                              manifest_sha256=pv.sha256_bytes(coh.to_csv_bytes()),
                              included=int(len(inc)))
        st = ex.prepare_structures(inc["parent_chembl_id"].astype(str),
                                   inc["canonical_smiles"].astype(str))
        return st, coh
    if args.manifest:
        df = _read_csv(args.manifest)
        if "decision" not in df.columns:
            raise ValueError(f"{args.manifest}: no 'decision' column (not a cohort manifest)")
        mesh = _manifest_mesh_ids(df)
        # [] (ATC / pooled manifest, blank or absent mesh_ids) is refused as "disease unknown"
        ex.check_modelling_allowed(mesh, args.plan_frozen)
        inc = df[df["decision"] == "include"]
        info["cohort"] = dict(source="manifest", path=str(args.manifest), mesh_ids=mesh,
                              sha256=pv.sha256_file(args.manifest), included=int(len(inc)))
        return ex.prepare_structures(inc["parent_chembl_id"], inc["canonical_smiles"]), None
    df = _read_csv(args.smiles_csv)
    for c in (args.id_column, args.smiles_column):
        if c not in df.columns:
            raise ValueError(f"{args.smiles_csv}: column {c!r} not found")
    info["cohort"] = dict(source="smiles_csv", path=str(args.smiles_csv),
                          sha256=pv.sha256_file(args.smiles_csv), rows=int(len(df)))
    return ex.prepare_structures(df[args.id_column], df[args.smiles_column]), None


def _target(args, st, info):
    """(y Series, name, rdkit2d exclusion key or None, failures) for exp-c."""
    import numpy as np
    import pandas as pd

    from topoqspr import experiments as ex
    from topoqspr import provenance as pv

    if args.target_csv:
        df = _read_csv(args.target_csv)
        for c in (args.target_id_column, args.target_column):
            if c not in df.columns:
                raise ValueError(f"{args.target_csv}: column {c!r} not found")
        if df[args.target_id_column].duplicated().any():
            raise ValueError(f"{args.target_csv}: duplicate ids")
        vals = pd.to_numeric(df.set_index(args.target_id_column)[args.target_column]
                             .replace("", np.nan), errors="coerce")
        y = vals.reindex(st.ids).astype(float)
        fails = [dict(row=i, feature=args.target_column, reason="no value in the target CSV")
                 for i in y.index[y.isna()]]
        info["target"] = dict(source="csv", path=str(args.target_csv),
                              sha256=pv.sha256_file(args.target_csv), column=args.target_column)
        return y, args.target_column, args.rdkit2d_exclusion, fails
    vals, fails, specs = ex.compute_targets(st, args.provider, [args.target])
    info["target"] = dict(source="provider", provider=args.provider, name=args.target,
                          spec=specs[0].as_dict())
    return vals[args.target], args.target, specs[0].rdkit2d_exclusion, fails


def _families(names: Sequence[str]):
    from topoqspr import experiments as ex

    return tuple(ex.POLY2_FAMILY if n == ex.POLY2_FAMILY.name else n for n in names)


# --------------------------------------------------------------- commands
def cmd_exp_a(args, info):
    from topoqspr import experiments as ex
    from topoqspr.cohort import build_cohort, naive_cohort

    if args.model and args.counts_only:
        raise ValueError("--model contradicts --counts-only")
    if not args.counts_only:  # refuse before any network access; the naive cohort's
        try:                  # headings are checked by experiment_a once it is retrieved
            ex.check_modelling_allowed(args.mesh, args.plan_frozen)
        except ex.QuarantineError as exc:
            raise ex.QuarantineError(f"{exc}; exp-a on this disease: pass --counts-only") from None
    client = _client(args)
    ctl = build_cohort(args.mesh, client)
    nav = naive_cohort(args.naive_text, client)
    info["cohort"] = dict(source="chembl", mesh_ids=ctl.parameters["mesh_ids"],
                          naive_text=args.naive_text, release=ctl.release,
                          release_date=ctl.release_date, disease_label=ctl.disease_label)
    res = ex.experiment_a(nav, ctl, target_provider=args.provider, target_names=args.targets,
                          index_names=args.index, top_k=args.top_k, alpha=args.alpha,
                          model=args.model, families=_families(args.families),
                          seeds=args.seeds, outer_splits=args.outer_splits,
                          inner_splits=args.inner_splits, n_jobs=args.n_jobs,
                          min_n=args.min_n, min_n_test=args.min_n_test,
                          mesh_ids=args.mesh, counts_only=args.counts_only,
                          plan_frozen=args.plan_frozen)
    return res, {"manifest_controlled.csv": ctl.to_csv_bytes(),
                 "manifest_naive.csv": nav.to_csv_bytes()}


def cmd_exp_c(args, info):
    from topoqspr import experiments as ex
    from topoqspr.models import LEGACY_MODELS

    st, _ = _load_structures(args, info)
    y, name, excl, fails = _target(args, st, info)
    X = ex.feature_set(st, args.feature_set, index_names=args.index, rdkit2d_exclusion=excl)
    legacy = None
    if args.legacy_families:
        unknown = sorted(set(args.legacy_families) - set(LEGACY_MODELS))
        if unknown:
            raise ValueError(f"unknown v0.x families {unknown}; known: {sorted(LEGACY_MODELS)}")
        legacy = {k: LEGACY_MODELS[k] for k in args.legacy_families}
    res = ex.experiment_c(X, y, connectivity=st.connectivity_labels(), mols=st.mols,
                          families=_families(args.families), legacy_families=legacy,
                          seeds=args.seeds, outer_splits=args.outer_splits,
                          inner_splits=args.inner_splits, alpha=args.alpha, n_jobs=args.n_jobs,
                          min_n=args.min_n, min_n_test=args.min_n_test,
                          label=f"{name} ~ {args.feature_set}")
    res.tables["molecules"] = st.table
    import pandas as pd

    res.tables["target_failures"] = pd.DataFrame(fails, columns=["row", "feature", "reason"])
    return res, {}


def cmd_exp_d(args, info):
    from topoqspr import experiments as ex
    from topoqspr import provenance as pv

    spec = ex.load_study_spec(args.spec)
    info["spec"] = dict(path=str(args.spec), sha256=pv.sha256_file(args.spec))
    res = ex.experiment_d(spec, families=_families(args.families), seeds=args.seeds,
                          outer_splits=args.outer_splits, inner_splits=args.inner_splits,
                          alpha=args.alpha, n_random=args.n_random,
                          n_permutations=args.n_permutations, random_seed=args.random_seed,
                          n_jobs=args.n_jobs, scaffold=not args.no_scaffold,
                          min_n=args.min_n, min_n_test=args.min_n_test)
    return res, {}


def cmd_size_dependence(args, info):
    import pandas as pd

    from topoqspr import experiments as ex

    st, _ = _load_structures(args, info)
    ti, fails = ex.index_table(st, args.index)
    res = ex.size_dependence(ti, ex.size_table(st), families=_families(args.families),
                             groups=ex.cv_groups(st, "random"), seeds=args.seeds,
                             outer_splits=args.outer_splits, inner_splits=args.inner_splits,
                             alpha=args.alpha, n_jobs=args.n_jobs, min_n=args.min_n,
                             min_n_test=args.min_n_test)
    res.tables["molecules"] = st.table
    res.tables["index_failures"] = pd.DataFrame(fails, columns=["row", "feature", "reason"])
    return res, {}


COMMANDS = {"exp-a": cmd_exp_a, "exp-c": cmd_exp_c, "exp-d": cmd_exp_d,
            "size-dependence": cmd_size_dependence}


# ------------------------------------------------------------------ parser
def build_parser() -> argparse.ArgumentParser:
    from topoqspr import experiments as ex

    ap = argparse.ArgumentParser(prog="run_experiment", description=__doc__.split("\n\n")[0])
    ap.add_argument("-v", "--verbose", action="store_true")
    sub = ap.add_subparsers(dest="command", required=True)

    def common(p, families_default):
        p.add_argument("--out", required=True, help="results directory (empty or new)")
        p.add_argument("--overwrite", action="store_true")
        p.add_argument("--timestamp", default=None, help="ISO time stamp for run_info.json")
        p.add_argument("--families", nargs="+", default=list(families_default),
                       help="model families chosen in the inner loop (topoqspr.models; "
                            f"'{ex.POLY2_FAMILY.name}' = degree-2 ridge)")
        p.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2, 3, 4])
        p.add_argument("--outer-splits", type=int, default=5)
        p.add_argument("--inner-splits", type=int, default=3)
        p.add_argument("--alpha", type=float, default=0.05)
        p.add_argument("--n-jobs", type=int, default=1, help="outer folds in parallel")
        p.add_argument("--min-n", type=int, default=50,
                       help="nested rows with n below this are flagged underpowered")
        p.add_argument("--min-n-test", type=int, default=10,
                       help="nested rows whose smallest outer test fold is below this are "
                            "flagged underpowered")
        p.add_argument("--plan-frozen", action="store_true",
                       help="lift the quarantine guard (only after the plan is frozen)")
        p.add_argument("--cache-dir", default=None)
        p.add_argument("--offline", action="store_true")
        p.add_argument("--release", default=None)

    def inputs(p):
        g = p.add_mutually_exclusive_group(required=True)
        g.add_argument("--mesh", nargs="+", help="confirmed MeSH descriptor IDs (controlled cohort)")
        g.add_argument("--manifest", help="cohort manifest CSV (topoqspr cohort)")
        g.add_argument("--smiles-csv", help="CSV with id and smiles columns")
        p.add_argument("--id-column", default="id")
        p.add_argument("--smiles-column", default="smiles")
        p.add_argument("--index", nargs="+", default=None, help="registry index names (default all)")

    p = sub.add_parser("exp-a", help="cohort sensitivity: naive vs controlled cohort")
    common(p, ("ridge", "pls", "knn"))
    p.add_argument("--mesh", nargs="+", required=True)
    p.add_argument("--naive-text", required=True, help="v0.x free-text disease query")
    p.add_argument("--provider", default="rdkit")
    p.add_argument("--targets", nargs="+", default=None)
    p.add_argument("--index", nargs="+", default=None)
    p.add_argument("--top-k", type=int, default=10)
    p.add_argument("--model", action="store_true",
                   help="also run nested dR2(TI - size) for both cohorts (quarantine applies)")
    p.add_argument("--counts-only", action="store_true",
                   help="cohort sizes, overlaps and size distributions only (no targets, "
                        "correlations or models): the only mode for a non-pilot disease "
                        "before the plan is frozen")

    p = sub.add_parser("exp-c", help="optimism chain: train -> v0.x CV -> nested -> scaffold")
    common(p, ("ridge", "pls", "knn"))
    inputs(p)
    p.add_argument("--feature-set", default="TI", choices=list(ex.FEATURE_SET_NAMES))
    p.add_argument("--provider", default="rdkit")
    tg = p.add_mutually_exclusive_group(required=True)
    tg.add_argument("--target", help="target name from --provider")
    tg.add_argument("--target-csv", help="CSV with measured values")
    p.add_argument("--target-column", default="value")
    p.add_argument("--target-id-column", default="id")
    p.add_argument("--rdkit2d-exclusion", default=None,
                   help="features.TARGET_EXCLUSIONS key (needed for rdkit2d with --target-csv)")
    p.add_argument("--legacy-families", nargs="+", default=None,
                   help="subset of the v0.x families (default: all ten)")

    p = sub.add_parser("exp-d", help="replication harness for one study spec")
    common(p, ("ridge", "pls", "linear"))
    p.add_argument("--spec", required=True, help="study spec (JSON, or YAML with PyYAML)")
    p.add_argument("--n-random", type=int, default=1000)
    p.add_argument("--n-permutations", type=int, default=1000)
    p.add_argument("--random-seed", type=int, default=0)
    p.add_argument("--no-scaffold", action="store_true")

    p = sub.add_parser("size-dependence", help="CV R2 of every index on size features")
    common(p, ("linear", ex.POLY2_FAMILY.name, "knn"))
    inputs(p)
    return ap


def _prepare_out(out: Path, overwrite: bool) -> Path:
    if out.exists() and any(out.iterdir()) and not overwrite:
        raise FileExistsError(f"{out} is not empty (pass --overwrite)")
    stage = out.parent / f".{out.name}.partial"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    return stage


def _publish(stage: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for p in list(out.iterdir()):
        shutil.rmtree(p) if p.is_dir() else p.unlink()
    for p in sorted(stage.iterdir()):
        shutil.move(str(p), str(out / p.name))
    stage.rmdir()


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")
    from topoqspr import __version__
    from topoqspr import provenance as pv
    from topoqspr.chembl import CacheMiss
    from topoqspr.experiments import QuarantineError

    out = Path(args.out)
    stage = None
    info = dict(command=args.command, argv=list(sys.argv[1:] if argv is None else argv),
                created_utc=args.timestamp or _now_utc(), topoqspr_version=__version__,
                code=pv.git_info())
    try:
        stage = _prepare_out(out, args.overwrite)
        t0 = time.perf_counter()
        res, extra = COMMANDS[args.command](args, info)
        info["wall_time_s"] = round(time.perf_counter() - t0, 1)
        res.write(stage)
        for name, data in extra.items():
            (stage / name).write_bytes(data)
        if "cohort" in info and info["cohort"].get("source") in ("chembl", "manifest"):
            from topoqspr.chembl import CHEMBL_DATA_LICENCE

            info["data_licence"] = CHEMBL_DATA_LICENCE
        (stage / "run_info.json").write_text(json.dumps(info, sort_keys=True, indent=1,
                                                        default=str) + "\n", encoding="utf-8")
        _publish(stage, out)
        print(f"{args.command}: wrote {out} in {info['wall_time_s']} s")
        return EXIT_OK
    except (QuarantineError, ValueError, KeyError, FileNotFoundError, FileExistsError,
            CacheMiss) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_INPUT
    except Exception:  # noqa: BLE001
        traceback.print_exc()
        return EXIT_INTERNAL
    finally:
        if stage is not None and stage.exists():
            shutil.rmtree(stage, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
