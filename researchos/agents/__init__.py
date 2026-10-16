"""Agent roles. v0.1 implements three deterministic roles; the rest are declared interfaces only."""
from __future__ import annotations

ROLES = [
    {"role": "coordinator", "implemented": True, "kind": "deterministic",
     "does": "parses a task file's structured block into a step DAG, records it, classifies claims; "
             "may downgrade a claim, never upgrades one to supported or approves its own work",
     "module": "researchos.agents.coordinator, researchos.adapters.classify_claim"},
    {"role": "verifier/reviewer", "implemented": True, "kind": "deterministic",
     "does": "hash and provenance checks, re-derivation of archived statistics, learner audit, "
             "rule-based adversarial objections",
     "module": "researchos.adapters (verify_external_points, learner_family_audit, adversarial_review, "
               "fold_level_contrast)"},
    {"role": "repository steward", "implemented": True, "kind": "deterministic",
     "does": "read-only Git state, artifact ingest, protected-process discovery, machine health; cannot push",
     "module": "researchos.inspectors, researchos.adapters.ingest_artifacts"},
    {"role": "mathematics", "implemented": False, "kind": "interface only",
     "does": "planned: SymPy/Lean cross-checks of index identities", "module": "-"},
    {"role": "statistics", "implemented": False, "kind": "interface only",
     "does": "planned: protocol and leakage review beyond the verifier's fixed rules", "module": "-"},
    {"role": "experiment", "implemented": False, "kind": "interface only (dry-run description only)",
     "does": "compute steps are described by dry-run and never launched in v0.1",
     "module": "researchos.adapters.describe_compute"},
    {"role": "publication", "implemented": False, "kind": "interface only",
     "does": "planned: claim-to-manuscript mapping", "module": "-"},
    {"role": "llm reviewer (Claude Code adapter)", "implemented": False, "kind": "stub",
     "does": "detects the `claude` CLI; v0.1 never sends repository content to a model",
     "module": "researchos.adapters.llm_review"},
]
