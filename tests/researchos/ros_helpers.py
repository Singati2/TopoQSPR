"""Shared helpers for ResearchOS tests (imported by conftest and tests)."""


def task_md(steps_json: str, prose: str = "Context only.", claims: str = '["C-X"]',
            question: str = "Why does A beat B?") -> str:
    return (f"# t\n\n{prose}\n\n```researchos-task\n"
            f'{{"title": "t", "question": "{question}", "claims": {claims}, "steps": {steps_json}}}\n```\n')


STEPS = """[
 {"id": "ingest", "adapter": "ingest_artifacts", "scope": "read-only",
  "params": {"dir": "out/ext", "files": ["summary.csv", "paired.csv"], "sums": "SHA256SUMS"}},
 {"id": "points", "adapter": "verify_external_points", "scope": "read-only", "deps": ["ingest"],
  "params": {"summary": "out/ext/summary.csv", "paired": "out/ext/paired.csv", "recomputed": "out/recomputed.csv"}},
 {"id": "learners", "adapter": "learner_family_audit", "scope": "read-only", "deps": ["ingest"],
  "params": {"set_a": "A_comp", "set_b": "B_comp", "external_summary": "out/ext/summary.csv",
             "internal_summaries": ["out/internal.csv"]}},
 {"id": "classify", "adapter": "classify_claim", "scope": "read-only", "deps": ["learners"],
  "params": {"claim_id": "C-X", "from_step": "learners", "margin": 0.01}},
 {"id": "folds", "adapter": "fold_level_contrast", "scope": "analysis-only", "deps": ["learners"],
  "params": {"folds": "out/folds.csv", "set_a": "A_comp", "set_b": "B_comp"}},
 {"id": "fixed", "adapter": "fixed_learner_experiment", "scope": "compute", "deps": ["classify"],
  "params": {"protocol_inputs": ["out/internal.csv"], "protocol_draft": "draft", "resource_estimate": "small"}}
]"""
