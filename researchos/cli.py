"""`python -m researchos <command>`. Exit codes: 0 ok, 1 check failed / blocked, 2 usage or config error."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from researchos import __version__, config, inspectors
from researchos import adapters as A
from researchos import executor as X
from researchos.agents import ROLES
from researchos.agents.coordinator import parse_task_file
from researchos.evidence import verify_manifest
from researchos.report import task_report
from researchos.safety import redact
from researchos.schema import Claim, Evidence, ValidationError
from researchos.store import Store


def _out(obj) -> None:
    print(redact(json.dumps(obj, indent=1, default=str) if not isinstance(obj, str) else obj))


def cmd_doctor(cfg, store, a) -> int:
    """Read-only environment check. Fails closed on anything it cannot establish."""
    checks = []

    def chk(name, ok, detail=""):
        checks.append({"check": name, "ok": bool(ok), "detail": detail})

    chk("python >= 3.10", sys.version_info >= (3, 10), sys.version.split()[0])
    try:
        st = inspectors.repo_state(cfg.repo, remote=a.remote)
        chk("git repository readable", True, f"{st['branch']} @ {st['head'][:9]}")
        chk("working tree clean", st["clean"], f"{st['dirty_paths']} changed paths")
        if a.remote:
            chk("remote heads fetched", "remote_heads" in st, st.get("remote_heads_error", ""))
    except inspectors.GitError as e:
        chk("git repository readable", False, str(e))
    chk("state store opens (schema migrated)", True, str(cfg.db_path))
    ign = Path(cfg.repo) / ".gitignore"
    chk("state dir is git-ignored", ign.exists() and ".researchos/" in ign.read_text(), ".researchos/ in .gitignore")
    procs = inspectors.processes(cfg.protected_processes)
    chk("protected processes (observed, not touched)", True, f"{len(procs)} running")
    h = inspectors.health()
    load_ok = h.get("load1") is not None and h["load1"] / max(1, h["cpus"]) <= cfg.max_load_per_cpu
    chk("load within compute ceiling", load_ok, f"load1 {h.get('load1')} on {h['cpus']} cpus")
    chk("memory available", (h.get("mem_available_gb") or 0) >= cfg.min_free_mem_gb, f"{h.get('mem_available_gb')} GB")
    chk("disk free", h["disk_free_gb"] >= cfg.min_free_disk_gb, f"{h['disk_free_gb']} GB")
    t = inspectors.tools()
    chk("claude CLI (optional; never invoked by v0.1)", True, "present" if t["claude"] else "absent")
    chk("SLURM (optional; no adapter in v0.1)", True, "present" if t["sbatch"] else "absent on this host")
    _out({"researchos": __version__, "checks": checks, "protected_processes": procs, "health": h})
    required = [c for c in checks if c["check"] in ("python >= 3.10", "git repository readable",
                                                   "state store opens (schema migrated)", "state dir is git-ignored")]
    return 0 if all(c["ok"] for c in required) else 1


def cmd_status(cfg, store, a) -> int:
    st = inspectors.repo_state(cfg.repo)
    tasks = [{"task_id": t["task_id"], "state": t["state"], "attempt": t["attempt"],
              "title": json.loads(t["spec"])["title"]} for t in store.tasks()]
    _out({"repo": st, "claims": store.claims(), "tasks": tasks,
          "protected_processes": inspectors.processes(cfg.protected_processes)})
    return 0


def cmd_agents(cfg, store, a) -> int:
    _out({"roles": ROLES, "adapters": {k: {"scope": v.scope.value, "role": v.role, "doc": v.doc}
                                       for k, v in A.REGISTRY.items()}})
    return 0


def cmd_claims(cfg, store, a) -> int:
    if a.action == "list":
        _out(store.claims())
    elif a.action == "show":
        c = store.claim(a.claim_id)
        if c is None:
            print(f"unknown claim {a.claim_id}", file=sys.stderr)
            return 1
        _out(c)
    elif a.action == "import":
        data = json.loads(Path(a.path).read_text())
        for c in data["claims"]:
            store.upsert_claim(Claim(**{k: v for k, v in c.items() if k != "evidence"}), actor=a.actor)
            for e in c.get("evidence", []):
                store.add_evidence(Evidence(claim_id=c["claim_id"], **e), actor=a.actor)
        _out({"imported": [c["claim_id"] for c in data["claims"]]})
    return 0


def cmd_tasks(cfg, store, a) -> int:
    _out([{"task_id": t["task_id"], "state": t["state"], "attempt": t["attempt"],
           "title": json.loads(t["spec"])["title"]} for t in store.tasks()])
    return 0


def cmd_task(cfg, store, a) -> int:
    if a.action == "propose":
        spec, flags = parse_task_file(Path(a.from_file))
        tid, new = store.add_task(spec, actor=a.actor)
        _out({"task_id": tid, "created": new, "steps": [s.id for s in spec.steps],
              "untrusted_instruction_like_text": flags})
        return 0
    if a.action == "dry-run":
        d = X.dry_run(store, cfg, a.task_id)
        _out(d)
        return 0
    if a.action == "approve":
        X.approve(store, cfg, a.task_id, a.scope, a.confirm, interactive=a.confirm is None and X.is_tty())
        _out({"approved": a.task_id, "scope": a.scope})
        return 0
    if a.action == "run":
        r = X.run(store, cfg, a.task_id)
        _out({k: v for k, v in r.items() if k != "steps"} |
             {"steps": {k: v["status"] for k, v in r["steps"].items()}})
        return 0 if not r.get("failed") else 1
    return 2


def cmd_report(cfg, store, a) -> int:
    md = task_report(store, cfg, a.task_id)
    if a.output:
        Path(a.output).write_text(md)
        print(a.output)
    else:
        print(md)
    return 0


def cmd_evidence(cfg, store, a) -> int:
    r = verify_manifest(cfg, Path(a.manifest))
    _out(r)
    return 0 if r["ok"] else 1


def cmd_audit(cfg, store, a) -> int:
    _out(store.audit_log(a.limit))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="researchos", description="TopoQSPR ResearchOS v0.1 (read-only by default)")
    p.add_argument("--repo", help="repository root (default: git toplevel of cwd)")
    p.add_argument("--state-dir", help="store location (default: <repo>/.researchos or $RESEARCHOS_HOME)")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("doctor"); d.add_argument("--remote", action="store_true", help="also read remote heads")
    d.set_defaults(fn=cmd_doctor)
    sub.add_parser("status").set_defaults(fn=cmd_status)
    ag = sub.add_parser("agents"); ag.add_argument("action", choices=["list"]); ag.set_defaults(fn=cmd_agents)
    c = sub.add_parser("claims")
    c.add_argument("action", choices=["list", "show", "import"])
    c.add_argument("claim_id", nargs="?")
    c.add_argument("--path", help="claims JSON for import")
    c.add_argument("--actor", default="owner-import")
    c.set_defaults(fn=cmd_claims)
    sub.add_parser("tasks").add_argument("action", choices=["list"])
    sub.choices["tasks"].set_defaults(fn=cmd_tasks)
    t = sub.add_parser("task")
    t.add_argument("action", choices=["propose", "dry-run", "approve", "run"])
    t.add_argument("task_id", nargs="?")
    t.add_argument("--from-file")
    t.add_argument("--scope", default="analysis-only", choices=["read-only", "analysis-only", "compute"])
    t.add_argument("--confirm", help="non-interactive approval: repeat the task id")
    t.add_argument("--actor", default="cli")
    t.set_defaults(fn=cmd_task)
    r = sub.add_parser("report"); r.add_argument("task_id"); r.add_argument("--format", choices=["md"], default="md")
    r.add_argument("--output"); r.set_defaults(fn=cmd_report)
    e = sub.add_parser("evidence"); e.add_argument("action", choices=["verify"]); e.add_argument("manifest")
    e.set_defaults(fn=cmd_evidence)
    au = sub.add_parser("audit"); au.add_argument("--limit", type=int, default=30); au.set_defaults(fn=cmd_audit)
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    try:
        cfg = config.load(Path(a.repo) if a.repo else None, Path(a.state_dir) if a.state_dir else None)
        store = Store(cfg.db_path)
    except (config.ConfigError, Exception) as e:  # fail closed on unknown repo or store state
        print(f"researchos: cannot start: {e}", file=sys.stderr)
        return 2
    if getattr(a, "action", None) in ("show",) and not a.claim_id:
        print("claims show needs a claim id", file=sys.stderr)
        return 2
    if a.cmd == "task" and a.action != "propose" and not a.task_id:
        print("task id required", file=sys.stderr)
        return 2
    if a.cmd == "task" and a.action == "propose" and not a.from_file:
        print("task propose needs --from-file", file=sys.stderr)
        return 2
    try:
        return a.fn(cfg, store, a)
    except (X.GateError, X.Busy) as e:
        print(f"refused: {e}", file=sys.stderr)
        return 1
    except (ValidationError, FileNotFoundError, KeyError, ValueError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    except Exception as e:  # noqa: BLE001 - CLI boundary; report, do not hide
        print(f"error: {type(e).__name__}: {e}", file=sys.stderr)
        return 1
