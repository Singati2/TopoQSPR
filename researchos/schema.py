"""Typed records and vocabularies. Validation is explicit so a bad record fails before it is stored."""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from enum import Enum


class ClaimStatus(str, Enum):
    PROPOSED = "proposed"
    UNDER_TEST = "under_test"
    SUPPORTED = "supported"
    REFUTED = "refuted"
    INCONCLUSIVE = "inconclusive"
    BLOCKED = "blocked"
    SUPERSEDED = "superseded"


class Independence(str, Enum):
    """How independent a piece of evidence is from the result it checks. Ordered weakest first."""
    QUALITATIVE_LLM_REVIEW = "qualitative_llm_review"      # a model read something and judged it
    OWNER_ASSERTION = "owner_assertion"
    SAME_CODE_RERUN = "same_code_rerun"                    # identical code re-executed
    CONSISTENCY_CHECK = "consistency_check"                # archived files agree with each other
    INDEPENDENT_RECOMPUTATION = "independent_recomputation"  # different code, same inputs
    ALTERNATIVE_IMPLEMENTATION = "alternative_implementation"
    SOURCE_DISJOINT_DATA = "source_disjoint_data"          # no shared records; independence not shown
    INDEPENDENT_DATA = "independent_data"                  # independence shown at reference level
    EXTERNAL_SOURCE = "external_source"                    # a published primary source

    @property
    def is_numerical(self) -> bool:
        return self not in (Independence.QUALITATIVE_LLM_REVIEW, Independence.OWNER_ASSERTION)


class Stage(str, Enum):
    EXPLORATORY = "exploratory"
    CONFIRMATORY = "confirmatory"
    POST_HOC = "post_hoc"
    NOT_EVALUATED = "not_evaluated"


class TaskState(str, Enum):
    PROPOSED = "proposed"
    APPROVED = "approved"
    QUEUED = "queued"
    RUNNING = "running"
    VERIFYING = "verifying"
    COMPLETED = "completed"
    FAILED = "failed"
    BLOCKED = "blocked"


class Scope(str, Enum):
    READ_ONLY = "read-only"          # in-process, allowlisted committed artifacts; runs without approval
    ANALYSIS_ONLY = "analysis-only"  # new statistics from existing outputs or paths outside the allowlist
    COMPUTE = "compute"              # training / refits; never launched by v0.1


class ResultStatus(str, Enum):
    PASS = "pass"
    FAIL = "fail"
    FLAG = "flag"          # ran, and found something a human must look at
    NOT_RUN = "not_run"
    ERROR = "error"


def canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def digest(obj) -> str:
    return hashlib.sha256(canonical(obj).encode()).hexdigest()


class ValidationError(ValueError):
    pass


def _req(cond: bool, msg: str) -> None:
    if not cond:
        raise ValidationError(msg)


@dataclass
class Claim:
    claim_id: str
    text: str
    scope: str
    status: ClaimStatus = ClaimStatus.PROPOSED
    stage: Stage = Stage.EXPLORATORY
    limitations: list[str] = field(default_factory=list)
    manuscript_refs: list[str] = field(default_factory=list)

    def validate(self) -> "Claim":
        _req(bool(self.claim_id) and self.claim_id.replace("-", "").replace("_", "").isalnum(),
             f"bad claim_id {self.claim_id!r}")
        _req(len(self.text.strip()) > 10, "claim text too short")
        _req(bool(self.scope.strip()), "claim scope required")
        self.status = ClaimStatus(self.status)
        self.stage = Stage(self.stage)
        return self

    def to_row(self) -> dict:
        d = asdict(self)
        d["status"] = self.status.value
        d["stage"] = self.stage.value
        return d


@dataclass
class Evidence:
    claim_id: str
    kind: str                      # e.g. "metric", "recomputation", "control", "review"
    independence: Independence
    summary: str
    supports: str                  # "supports" | "contradicts" | "neutral"
    reviewer: str                  # who or what produced it (code path, person, model)
    method: str                    # how: script, rule, read
    artifacts: list[dict] = field(default_factory=list)   # [{path, sha256, commit}]
    metrics: dict = field(default_factory=dict)
    commit: str | None = None
    command: str | None = None
    seed: str | None = None
    protocol_hash: str | None = None
    run_id: str | None = None
    dataset: str | None = None
    license_note: str | None = None
    limitations: list[str] = field(default_factory=list)

    def validate(self) -> "Evidence":
        self.independence = Independence(self.independence)
        _req(self.supports in ("supports", "contradicts", "neutral"), f"bad supports {self.supports!r}")
        _req(bool(self.reviewer) and bool(self.method), "reviewer and method required")
        for a in self.artifacts:
            _req("path" in a and "sha256" in a, f"artifact needs path and sha256: {a}")
            _req(len(a["sha256"]) == 64, f"sha256 must be 64 hex chars: {a}")
        if self.supports == "supports" and self.independence == Independence.QUALITATIVE_LLM_REVIEW:
            # A model's opinion is recorded, but it can never be the numerical basis for "supported".
            self.limitations = list(self.limitations) + ["qualitative review: not numerical verification"]
        return self

    def evidence_id(self) -> str:
        d = asdict(self)
        d["independence"] = self.independence.value
        return "E-" + digest(d)[:12]

    def to_row(self) -> dict:
        d = asdict(self)
        d["independence"] = self.independence.value
        return d


@dataclass
class Step:
    id: str
    adapter: str
    scope: Scope
    params: dict = field(default_factory=dict)
    deps: list[str] = field(default_factory=list)
    description: str = ""


@dataclass
class TaskSpec:
    title: str
    question: str
    steps: list[Step]
    claims: list[str] = field(default_factory=list)
    untrusted_notes: str = ""

    def validate(self, known_adapters: dict) -> "TaskSpec":
        _req(bool(self.title.strip()) and bool(self.question.strip()), "title and question required")
        _req(len(self.steps) > 0, "a task needs at least one step")
        ids = [s.id for s in self.steps]
        _req(len(ids) == len(set(ids)), "duplicate step ids")
        for s in self.steps:
            s.scope = Scope(s.scope)
            _req(s.adapter in known_adapters, f"unknown adapter {s.adapter!r} (only registered adapters run)")
            want = known_adapters[s.adapter].scope
            _req(s.scope == want, f"step {s.id}: adapter {s.adapter} has scope {want.value}, not {s.scope.value}")
            for d in s.deps:
                _req(d in ids and d != s.id, f"step {s.id}: bad dependency {d!r}")
        _topo(self.steps)  # raises on cycles
        return self

    def to_dict(self) -> dict:
        return {"title": self.title, "question": self.question, "claims": self.claims,
                "steps": [{**asdict(s), "scope": Scope(s.scope).value} for s in self.steps]}

    def task_id(self) -> str:
        # untrusted_notes are excluded: changing prose must not change what runs or its identity
        return "T-" + digest(self.to_dict())[:12]

    @staticmethod
    def from_dict(d: dict) -> "TaskSpec":
        _req(isinstance(d, dict), "task spec must be a JSON object")
        steps = [Step(id=s["id"], adapter=s["adapter"], scope=s.get("scope", "read-only"),
                      params=s.get("params", {}), deps=s.get("deps", []),
                      description=s.get("description", "")) for s in d.get("steps", [])]
        return TaskSpec(title=d.get("title", ""), question=d.get("question", ""), steps=steps,
                        claims=d.get("claims", []), untrusted_notes=d.get("untrusted_notes", ""))


def _topo(steps: list[Step]) -> list[str]:
    deps = {s.id: set(s.deps) for s in steps}
    order: list[str] = []
    while deps:
        ready = sorted(k for k, v in deps.items() if not v)
        _req(bool(ready), "dependency cycle in task steps")
        for k in ready:
            order.append(k)
            deps.pop(k)
        for v in deps.values():
            v.difference_update(ready)
    return order


def topo_order(steps: list[Step]) -> list[str]:
    return _topo(steps)
