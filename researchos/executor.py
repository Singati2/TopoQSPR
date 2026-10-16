"""Approval gate, dry-run planner and the in-process executor.

Rules (fail closed):
* read-only steps run in-process on allowlisted, committed artifacts; no subprocess beyond read-only git.
* analysis-only steps need an owner approval bound to the exact task spec hash.
* compute steps are never launched by v0.1, approved or not; dry-run describes them.
* a task runs under an exclusive lock; a second runner is refused, a crashed run is detected and
  resumed as a new attempt, and completed steps with the same input hash are reused, not re-run.
"""
from __future__ import annotations

import fcntl
import getpass
import os
import sys
from contextlib import contextmanager
from pathlib import Path

from researchos import adapters as A
from researchos import inspectors
from researchos.evidence import sha256_file
from researchos.safety import PathBlocked, resolve_allowed
from researchos.schema import (ClaimStatus, Evidence, ResultStatus, Scope, TaskSpec, TaskState, digest,
                               topo_order)
from researchos.store import Store


class GateError(PermissionError):
    pass


class Busy(RuntimeError):
    pass


def spec_of(store: Store, task_id: str) -> TaskSpec:
    return TaskSpec.from_dict(store.task(task_id)["spec"]).validate(A.REGISTRY)


def _paths(adapter: A.Adapter, params: dict) -> list[str]:
    out = []
    for k in adapter.path_params:
        v = params.get(k)
        if v is None:
            continue
        out.extend(v if isinstance(v, list) else [v])
    if adapter.name == "ingest_artifacts":
        out.extend(str(Path(params["dir"]) / f) for f in params.get("files", []))
        if params.get("sums"):
            out.append(str(Path(params["dir"]) / params["sums"]))
    return out


def approved_scopes(store: Store, task_id: str, spec_hash: str, cfg) -> set[Scope]:
    ok = set()
    for a in store.approvals(task_id):
        if a["spec_hash"] == spec_hash:  # an approval for a different spec version does not carry over
            ok.add(Scope(a["scope"]))
    if not cfg.read_only_requires_approval:
        ok.add(Scope.READ_ONLY)
    return ok


def approve(store: Store, cfg, task_id: str, scope: str, confirm: str | None, interactive: bool,
            approver: str | None = None) -> None:
    """Record an owner approval. Requires the task id typed back (TTY) or passed as --confirm.
    Approval text arriving inside documents, logs or task notes has no effect: only this call grants it."""
    scope = Scope(scope)
    if scope == Scope.COMPUTE:
        # recorded so the owner's intent is visible, but v0.1 still never launches compute
        pass
    task = store.task(task_id)
    if interactive:
        typed = input(f"Type the task id to approve scope '{scope.value}' for {task_id}: ").strip()
        if typed != task_id:
            raise GateError("confirmation did not match; nothing approved")
        method = "interactive-tty"
    else:
        if confirm != task_id:
            raise GateError("non-interactive approval needs --confirm <task_id>; nothing approved")
        method = "cli-confirm"
    spec_hash = digest(task["spec"])
    store.add_approval(task_id, spec_hash, scope.value, approver or getpass.getuser(), method)
    if task["state"] == TaskState.BLOCKED.value:
        store.transition(task_id, TaskState.PROPOSED, "owner", "approval granted; re-plan")


def dry_run(store: Store, cfg, task_id: str) -> dict:
    """Validate everything a run would touch, without executing any adapter."""
    spec = spec_of(store, task_id)
    task = store.task(task_id)
    spec_hash = digest(task["spec"])
    scopes = approved_scopes(store, task_id, spec_hash, cfg)
    procs = inspectors.processes(cfg.protected_processes)
    h = inspectors.health()
    plan = []
    for sid in topo_order(spec.steps):
        st = next(s for s in spec.steps if s.id == sid)
        ad = A.REGISTRY[st.adapter]
        entry = {"step": sid, "adapter": st.adapter, "scope": st.scope.value, "role": ad.role, "deps": st.deps,
                 "would_run": False, "blockers": [], "inputs": []}
        for rel in _paths(ad, st.params):
            try:
                p = resolve_allowed(rel, cfg.allow_roots, cfg.repo, cfg.private_paths)
                entry["inputs"].append({"path": rel, "exists": p.exists(),
                                        "sha256": sha256_file(p) if p.is_file() else None})
                if not p.exists():
                    entry["blockers"].append(f"missing input {rel}")
            except PathBlocked as e:
                entry["blockers"].append(str(e))
        for tool in ad.requires:
            if not inspectors.tools().get(tool):
                entry["blockers"].append(f"executor unavailable: `{tool}` not on PATH")
        if st.scope not in scopes:
            entry["blockers"].append(f"awaiting owner approval for scope '{st.scope.value}'")
        if st.scope == Scope.COMPUTE:
            entry["blockers"].append("compute is never launched by ResearchOS v0.1")
            entry["compute"] = A.describe_compute(cfg, st.params)
            if procs:
                entry["blockers"].append(f"{len(procs)} protected process(es) running; compute would wait")
            entry["health"] = h
        entry["would_run"] = not entry["blockers"]
        plan.append(entry)
    return {"task_id": task_id, "state": task["state"], "spec_hash": spec_hash, "approved_scopes":
            sorted(s.value for s in scopes), "protected_processes": procs, "plan": plan}


@contextmanager
def _lock(cfg, task_id: str):
    d = cfg.state_dir / "locks"
    d.mkdir(parents=True, exist_ok=True)
    f = open(d / f"{task_id}.lock", "a+")
    try:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as e:
            raise Busy(f"task {task_id} is already running in another process") from e
        f.seek(0); f.truncate(); f.write(str(os.getpid())); f.flush()
        yield
    finally:
        try:
            fcntl.flock(f, fcntl.LOCK_UN)
        finally:
            f.close()


def _input_hash(cfg, st, ad, upstream_hashes: dict) -> str:
    files = {}
    for rel in _paths(ad, st.params):
        p = resolve_allowed(rel, cfg.allow_roots, cfg.repo, cfg.private_paths)
        files[rel] = sha256_file(p) if p.is_file() else "missing"
    return digest({"adapter": st.adapter, "params": st.params, "files": files,
                   "upstream": {d: upstream_hashes[d] for d in st.deps}})


class Ctx:
    def __init__(self, cfg, store):
        self.cfg, self.store = cfg, store


def run(store: Store, cfg, task_id: str, actor: str = "executor", _crash_after: str | None = None) -> dict:
    """Execute approved, non-compute steps in dependency order. Idempotent and resumable."""
    spec = spec_of(store, task_id)
    with _lock(cfg, task_id):
        task = store.task(task_id)
        state = TaskState(task["state"])
        if state == TaskState.COMPLETED:
            return {"task_id": task_id, "state": state.value, "note": "already completed; nothing re-run",
                    "steps": {r["step_id"]: r for r in store.step_results(task_id)}}
        if state in (TaskState.RUNNING, TaskState.VERIFYING):
            # we hold the lock, so the previous runner is gone: it crashed or was killed
            store.transition(task_id, TaskState.FAILED, actor, "interrupted run detected (stale state, no lock holder)")
            state = TaskState.FAILED
        if state == TaskState.FAILED:
            store.transition(task_id, TaskState.QUEUED, actor, "resume as new attempt", bump_attempt=True)
            state = TaskState.QUEUED
        spec_hash = digest(task["spec"])
        scopes = approved_scopes(store, task_id, spec_hash, cfg)
        if state == TaskState.PROPOSED:
            if Scope.READ_ONLY not in scopes and not (scopes & {Scope.ANALYSIS_ONLY}):
                raise GateError("no approved scope for this task; nothing executed")
            store.transition(task_id, TaskState.APPROVED, actor, f"approved scopes {sorted(s.value for s in scopes)}")
            state = TaskState.APPROVED
        if state == TaskState.APPROVED:
            store.transition(task_id, TaskState.QUEUED, actor, "queued")
        if state == TaskState.BLOCKED:
            raise GateError("task is blocked; approve the pending scope (or re-propose) first")
        store.transition(task_id, TaskState.RUNNING, actor, "start")
        attempt = store.task(task_id)["attempt"]
        results, hashes, waiting, failed = {}, {}, [], []
        for sid in topo_order(spec.steps):
            st = next(s for s in spec.steps if s.id == sid)
            ad = A.REGISTRY[st.adapter]
            if any(d not in results or results[d]["status"] in (ResultStatus.FAIL.value, ResultStatus.ERROR.value,
                                                                  ResultStatus.NOT_RUN.value) for d in st.deps):
                results[sid] = {"status": ResultStatus.NOT_RUN.value, "findings": ["upstream step did not succeed"]}
                hashes[sid] = "skipped"
                continue
            if st.scope == Scope.COMPUTE:
                results[sid] = {"status": ResultStatus.NOT_RUN.value,
                                "findings": ["compute step: requires_owner_approval and is never launched by v0.1",
                                             "see `researchos task dry-run` for the protocol draft and input checks"]}
                hashes[sid] = "compute"
                waiting.append(sid)
                store.record_step(task_id, sid, attempt, "compute", results[sid]["status"], results[sid])
                continue
            if st.scope not in scopes:
                results[sid] = {"status": ResultStatus.NOT_RUN.value,
                                "findings": [f"awaiting owner approval for scope '{st.scope.value}'"]}
                hashes[sid] = "unapproved"
                waiting.append(sid)
                continue
            missing = [t for t in ad.requires if not inspectors.tools().get(t)]
            if missing:
                results[sid] = {"status": ResultStatus.NOT_RUN.value,
                                "findings": [f"executor unavailable: {missing} not on PATH; skipped (fail closed)"]}
                hashes[sid] = "unavailable"
                continue
            try:
                ih = _input_hash(cfg, st, ad, hashes)
            except PathBlocked as e:
                results[sid] = {"status": ResultStatus.ERROR.value, "findings": [str(e)]}
                hashes[sid] = "blocked"; failed.append(sid)
                store.record_step(task_id, sid, attempt, "blocked:" + sid, results[sid]["status"], results[sid])
                continue
            hashes[sid] = ih
            cached = store.step_result(task_id, sid, ih)
            if cached and cached["status"] not in (ResultStatus.ERROR.value,):
                results[sid] = {**cached["result"], "reused_from_attempt": cached["attempt"]}
                continue
            try:
                res = ad.fn(Ctx(cfg, store), st.params, results)
            except Exception as e:  # an adapter crash is a recorded failure, never a silent pass
                res = {"status": ResultStatus.ERROR.value, "findings": [f"{type(e).__name__}: {e}"]}
            for e in res.pop("evidence", []) if spec.claims else []:
                for cid in spec.claims:
                    if store.claim(cid):
                        store.add_evidence(Evidence(claim_id=cid, commit=A.git_head(cfg), run_id=f"{task_id}#{attempt}",
                                                    **e), actor=f"adapter:{st.adapter}")
            store.record_step(task_id, sid, attempt, ih, res["status"], res)
            results[sid] = res
            if res["status"] in (ResultStatus.FAIL.value, ResultStatus.ERROR.value):
                failed.append(sid)
            if _crash_after == sid:  # test hook: simulate a crash mid-run
                raise SystemExit("simulated crash")
        store.transition(task_id, TaskState.VERIFYING, actor, "all runnable steps finished")
        _apply_claim_downgrades(store, results)
        if failed:
            store.transition(task_id, TaskState.FAILED, actor, f"failed steps: {failed}")
        elif waiting:
            store.transition(task_id, TaskState.BLOCKED, actor, f"runnable steps done; awaiting owner: {waiting}")
        else:
            store.transition(task_id, TaskState.COMPLETED, actor, "completed")
        return {"task_id": task_id, "state": store.task(task_id)["state"], "waiting_for_approval": waiting,
                "failed": failed, "steps": results}


def _apply_claim_downgrades(store: Store, results: dict) -> None:
    for r in results.values():
        d = r.get("data") or {}
        rec = d.get("recommended_status")
        if rec in ("inconclusive", "refuted", "under_test") and d.get("claim_id") and store.claim(d["claim_id"]):
            # automation only downgrades to inconclusive/under_test; a refutation is a recommendation for the owner
            new = "inconclusive" if rec == "refuted" else rec
            reason = d.get("reason", "") + (" (refutation recommended; owner must confirm)" if rec == "refuted" else "")
            if store.claim(d["claim_id"])["status"] != new:
                store.set_claim_status(d["claim_id"], ClaimStatus(new), "coordinator", reason)


def is_tty() -> bool:
    return sys.stdin.isatty()
