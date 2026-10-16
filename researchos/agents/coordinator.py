"""Coordinator: turns a task file into a validated step DAG.

A task file is Markdown. Exactly one fenced block tagged `researchos-task` holds the JSON spec.
Everything else in the file is *untrusted prose*: it is stored, screened and shown, never parsed
for instructions. Instruction-like text (shell commands, "owner approved", "ignore previous
instructions") is flagged in reports and has no effect.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from researchos import adapters as A
from researchos.safety import scan_untrusted
from researchos.schema import TaskSpec, ValidationError

_BLOCK = re.compile(r"```researchos-task\s*\n(.*?)\n```", re.S)


def parse_task_file(path: Path) -> tuple[TaskSpec, list[dict]]:
    text = Path(path).read_text()
    blocks = _BLOCK.findall(text)
    if len(blocks) != 1:
        raise ValidationError(f"expected exactly one ```researchos-task block, found {len(blocks)}")
    try:
        d = json.loads(blocks[0])
    except json.JSONDecodeError as e:
        raise ValidationError(f"task block is not valid JSON: {e}") from e
    prose = _BLOCK.sub("[task block]", text)
    spec = TaskSpec.from_dict(d)
    spec.untrusted_notes = prose
    spec.validate(A.REGISTRY)
    flags = scan_untrusted(prose) + scan_untrusted(json.dumps(d.get("question", "")))
    return spec, flags
