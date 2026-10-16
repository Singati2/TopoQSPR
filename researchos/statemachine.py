"""Task lifecycle. Result status and scientific-claim status are tracked separately (schema.py)."""
from __future__ import annotations

from researchos.schema import TaskState as S

ALLOWED: dict[S, set[S]] = {
    S.PROPOSED: {S.APPROVED, S.BLOCKED},
    S.APPROVED: {S.QUEUED, S.BLOCKED},
    S.QUEUED: {S.RUNNING, S.BLOCKED},
    S.RUNNING: {S.VERIFYING, S.FAILED, S.BLOCKED},
    S.VERIFYING: {S.COMPLETED, S.FAILED, S.BLOCKED},  # BLOCKED: some steps await owner approval
    S.FAILED: {S.QUEUED},            # retry = new attempt, never silent re-execution
    S.BLOCKED: {S.PROPOSED},         # unblocking sends the task back for approval
    S.COMPLETED: set(),
}


class TransitionError(RuntimeError):
    pass


def check(old: S | str, new: S | str) -> None:
    old, new = S(old), S(new)
    if new not in ALLOWED[old]:
        raise TransitionError(f"invalid task transition {old.value} -> {new.value}")
