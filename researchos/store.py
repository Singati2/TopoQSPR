"""SQLite store: claims, append-only evidence, tasks, step results, approvals and an append-only audit log.

One file per project (default `<repo>/.researchos/state.db`, git-ignored). Every write is one
transaction. Evidence and audit rows cannot be updated or deleted (enforced by triggers).
"""
from __future__ import annotations

import json
import sqlite3
import time
from contextlib import contextmanager
from pathlib import Path

from researchos import statemachine
from researchos.schema import Claim, ClaimStatus, Evidence, TaskSpec, TaskState

SCHEMA_VERSION = 1

MIGRATIONS = {
    1: """
    CREATE TABLE claims (claim_id TEXT PRIMARY KEY, text TEXT NOT NULL, scope TEXT NOT NULL,
        status TEXT NOT NULL, stage TEXT NOT NULL, limitations TEXT NOT NULL,
        manuscript_refs TEXT NOT NULL, created REAL NOT NULL, updated REAL NOT NULL);
    CREATE TABLE evidence (evidence_id TEXT PRIMARY KEY, claim_id TEXT NOT NULL REFERENCES claims,
        record TEXT NOT NULL, created REAL NOT NULL);
    CREATE TABLE tasks (task_id TEXT PRIMARY KEY, spec TEXT NOT NULL, untrusted_notes TEXT NOT NULL,
        state TEXT NOT NULL, attempt INTEGER NOT NULL DEFAULT 0, created REAL NOT NULL, updated REAL NOT NULL);
    CREATE TABLE step_results (task_id TEXT NOT NULL, step_id TEXT NOT NULL, attempt INTEGER NOT NULL,
        input_hash TEXT NOT NULL, status TEXT NOT NULL, result TEXT NOT NULL, finished REAL NOT NULL,
        PRIMARY KEY (task_id, step_id, input_hash));
    CREATE TABLE approvals (task_id TEXT NOT NULL, spec_hash TEXT NOT NULL, scope TEXT NOT NULL,
        approver TEXT NOT NULL, method TEXT NOT NULL, created REAL NOT NULL,
        PRIMARY KEY (task_id, scope));
    CREATE TABLE audit (seq INTEGER PRIMARY KEY AUTOINCREMENT, ts REAL NOT NULL, actor TEXT NOT NULL,
        action TEXT NOT NULL, target TEXT NOT NULL, detail TEXT NOT NULL);
    CREATE TRIGGER evidence_no_update BEFORE UPDATE ON evidence BEGIN SELECT RAISE(ABORT, 'evidence is append-only'); END;
    CREATE TRIGGER evidence_no_delete BEFORE DELETE ON evidence BEGIN SELECT RAISE(ABORT, 'evidence is append-only'); END;
    CREATE TRIGGER audit_no_update BEFORE UPDATE ON audit BEGIN SELECT RAISE(ABORT, 'audit is append-only'); END;
    CREATE TRIGGER audit_no_delete BEFORE DELETE ON audit BEGIN SELECT RAISE(ABORT, 'audit is append-only'); END;
    """,
}


class StoreError(RuntimeError):
    pass


class Store:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(self.path, isolation_level=None, timeout=30)
        self.db.row_factory = sqlite3.Row
        self.db.execute("PRAGMA foreign_keys = ON")
        self.db.execute("PRAGMA journal_mode = WAL")
        self._migrate()

    # -- infrastructure -------------------------------------------------------------------
    @contextmanager
    def tx(self):
        self.db.execute("BEGIN IMMEDIATE")
        try:
            yield self.db
        except BaseException:
            self.db.execute("ROLLBACK")
            raise
        else:
            self.db.execute("COMMIT")

    def _migrate(self) -> None:
        v = self.db.execute("PRAGMA user_version").fetchone()[0]
        if v > SCHEMA_VERSION:
            raise StoreError(f"store schema v{v} is newer than this ResearchOS (v{SCHEMA_VERSION}); refusing to open")
        for n in range(v + 1, SCHEMA_VERSION + 1):
            with self.tx() as c:  # statement by statement: executescript() would commit implicitly
                for stmt in _split(MIGRATIONS[n]):
                    c.execute(stmt)
                c.execute(f"PRAGMA user_version = {n}")

    def audit(self, c, actor: str, action: str, target: str, detail: dict | None = None) -> None:
        c.execute("INSERT INTO audit (ts, actor, action, target, detail) VALUES (?,?,?,?,?)",
                  (time.time(), actor, action, target, json.dumps(detail or {}, sort_keys=True)))

    def audit_log(self, limit: int = 50) -> list[dict]:
        rows = self.db.execute("SELECT * FROM audit ORDER BY seq DESC LIMIT ?", (limit,)).fetchall()
        return [dict(r) for r in rows]

    # -- claims and evidence --------------------------------------------------------------
    def upsert_claim(self, claim: Claim, actor: str) -> None:
        claim.validate()
        r = claim.to_row()
        now = time.time()
        with self.tx() as c:
            old = c.execute("SELECT status FROM claims WHERE claim_id=?", (claim.claim_id,)).fetchone()
            if old is None:
                c.execute("INSERT INTO claims VALUES (?,?,?,?,?,?,?,?,?)",
                          (r["claim_id"], r["text"], r["scope"], r["status"], r["stage"],
                           json.dumps(r["limitations"]), json.dumps(r["manuscript_refs"]), now, now))
                self.audit(c, actor, "claim.create", claim.claim_id, {"status": r["status"]})
            else:
                c.execute("UPDATE claims SET text=?, scope=?, stage=?, limitations=?, manuscript_refs=?, updated=? "
                          "WHERE claim_id=?", (r["text"], r["scope"], r["stage"], json.dumps(r["limitations"]),
                                               json.dumps(r["manuscript_refs"]), now, claim.claim_id))
                self.audit(c, actor, "claim.update", claim.claim_id, {})

    def set_claim_status(self, claim_id: str, status: ClaimStatus, actor: str, reason: str) -> None:
        status = ClaimStatus(status)
        with self.tx() as c:
            old = c.execute("SELECT status FROM claims WHERE claim_id=?", (claim_id,)).fetchone()
            if old is None:
                raise StoreError(f"unknown claim {claim_id}")
            c.execute("UPDATE claims SET status=?, updated=? WHERE claim_id=?", (status.value, time.time(), claim_id))
            self.audit(c, actor, "claim.status", claim_id, {"from": old["status"], "to": status.value, "reason": reason})

    def claim(self, claim_id: str) -> dict | None:
        r = self.db.execute("SELECT * FROM claims WHERE claim_id=?", (claim_id,)).fetchone()
        if r is None:
            return None
        d = dict(r)
        d["limitations"] = json.loads(d["limitations"])
        d["manuscript_refs"] = json.loads(d["manuscript_refs"])
        d["evidence"] = self.evidence_for(claim_id)
        return d

    def claims(self) -> list[dict]:
        return [dict(r) for r in self.db.execute("SELECT claim_id, status, stage, scope, text FROM claims ORDER BY claim_id")]

    def add_evidence(self, ev: Evidence, actor: str) -> str:
        ev.validate()
        eid = ev.evidence_id()
        with self.tx() as c:
            if c.execute("SELECT 1 FROM claims WHERE claim_id=?", (ev.claim_id,)).fetchone() is None:
                raise StoreError(f"evidence for unknown claim {ev.claim_id}")
            if c.execute("SELECT 1 FROM evidence WHERE evidence_id=?", (eid,)).fetchone() is None:
                c.execute("INSERT INTO evidence VALUES (?,?,?,?)",
                          (eid, ev.claim_id, json.dumps(ev.to_row(), sort_keys=True), time.time()))
                self.audit(c, actor, "evidence.add", eid, {"claim": ev.claim_id, "independence": ev.independence.value,
                                                           "supports": ev.supports})
        return eid

    def evidence_for(self, claim_id: str) -> list[dict]:
        rows = self.db.execute("SELECT evidence_id, record FROM evidence WHERE claim_id=? ORDER BY created", (claim_id,))
        return [{"evidence_id": r["evidence_id"], **json.loads(r["record"])} for r in rows]

    # -- tasks ----------------------------------------------------------------------------
    def add_task(self, spec: TaskSpec, actor: str) -> tuple[str, bool]:
        tid = spec.task_id()
        now = time.time()
        with self.tx() as c:
            if c.execute("SELECT 1 FROM tasks WHERE task_id=?", (tid,)).fetchone():
                return tid, False  # idempotent: same spec, same task
            c.execute("INSERT INTO tasks VALUES (?,?,?,?,?,?,?)",
                      (tid, json.dumps(spec.to_dict(), sort_keys=True), spec.untrusted_notes,
                       TaskState.PROPOSED.value, 0, now, now))
            self.audit(c, actor, "task.propose", tid, {"title": spec.title})
        return tid, True

    def task(self, task_id: str) -> dict:
        r = self.db.execute("SELECT * FROM tasks WHERE task_id=?", (task_id,)).fetchone()
        if r is None:
            raise StoreError(f"unknown task {task_id}")
        d = dict(r)
        d["spec"] = json.loads(d["spec"])
        return d

    def tasks(self) -> list[dict]:
        return [dict(r) for r in self.db.execute("SELECT task_id, state, attempt, spec FROM tasks ORDER BY created")]

    def transition(self, task_id: str, new: TaskState, actor: str, reason: str = "", bump_attempt: bool = False) -> None:
        with self.tx() as c:
            r = c.execute("SELECT state FROM tasks WHERE task_id=?", (task_id,)).fetchone()
            if r is None:
                raise StoreError(f"unknown task {task_id}")
            statemachine.check(r["state"], new)
            c.execute("UPDATE tasks SET state=?, updated=?, attempt=attempt+? WHERE task_id=?",
                      (TaskState(new).value, time.time(), 1 if bump_attempt else 0, task_id))
            self.audit(c, actor, "task.transition", task_id, {"from": r["state"], "to": TaskState(new).value,
                                                              "reason": reason})

    def record_step(self, task_id: str, step_id: str, attempt: int, input_hash: str, status: str, result: dict) -> None:
        with self.tx() as c:
            c.execute("INSERT OR REPLACE INTO step_results VALUES (?,?,?,?,?,?,?)",
                      (task_id, step_id, attempt, input_hash, status, json.dumps(result, sort_keys=True, default=str),
                       time.time()))
            self.audit(c, "executor", "step.result", f"{task_id}/{step_id}", {"status": status, "input_hash": input_hash})

    def step_result(self, task_id: str, step_id: str, input_hash: str) -> dict | None:
        r = self.db.execute("SELECT * FROM step_results WHERE task_id=? AND step_id=? AND input_hash=?",
                            (task_id, step_id, input_hash)).fetchone()
        if r is None:
            return None
        d = dict(r)
        d["result"] = json.loads(d["result"])
        return d

    def step_results(self, task_id: str) -> list[dict]:
        rows = self.db.execute("SELECT * FROM step_results WHERE task_id=? ORDER BY finished", (task_id,)).fetchall()
        out = []
        for r in rows:
            d = dict(r)
            d["result"] = json.loads(d["result"])
            out.append(d)
        return out

    # -- approvals ------------------------------------------------------------------------
    def add_approval(self, task_id: str, spec_hash: str, scope: str, approver: str, method: str) -> None:
        with self.tx() as c:
            c.execute("INSERT OR REPLACE INTO approvals VALUES (?,?,?,?,?,?)",
                      (task_id, spec_hash, scope, approver, method, time.time()))
            self.audit(c, approver, "approval.grant", task_id, {"scope": scope, "method": method, "spec_hash": spec_hash})

    def approvals(self, task_id: str) -> list[dict]:
        return [dict(r) for r in self.db.execute("SELECT * FROM approvals WHERE task_id=?", (task_id,))]


def _split(script: str) -> list[str]:
    """Split a migration into statements, keeping trigger bodies (BEGIN ... END;) intact."""
    out, buf = [], []
    for line in script.strip().splitlines():
        buf.append(line)
        joined = "\n".join(buf).strip()
        if joined.endswith(";") and (not joined.upper().startswith("CREATE TRIGGER") or joined.upper().endswith("END;")):
            out.append(joined)
            buf = []
    if "".join(buf).strip():
        raise StoreError("unterminated statement in migration")
    return out
