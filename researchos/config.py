"""Project configuration. Defaults are conservative; an optional `researchos.json` at the repo root
can narrow or extend them. Unknown keys are rejected (fail closed)."""
from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_PROTECTED_PROCESSES = [
    # the owner's and the study's long-running jobs; ResearchOS never touches them and refuses
    # compute while any is alive
    "test_baselines.py", "baselines.py", "bp_tests.py", "box_sensitivity.py", "e3_isomer.py",
    "bin/topoqspr ", "-m topoqspr", "joblib.externals.loky",
]
DEFAULT_PRIVATE_PATHS = ["~/topoqspr_private"]  # never read into reports, never allowlisted


@dataclass
class Config:
    repo: Path
    state_dir: Path
    allow_roots: list[Path]
    private_paths: list[Path]
    protected_processes: list[str] = field(default_factory=lambda: list(DEFAULT_PROTECTED_PROCESSES))
    max_compute_concurrency: int = 1
    max_load_per_cpu: float = 0.5
    min_free_mem_gb: float = 8.0
    min_free_disk_gb: float = 20.0
    read_only_requires_approval: bool = False

    @property
    def db_path(self) -> Path:
        return self.state_dir / "state.db"


class ConfigError(RuntimeError):
    pass


def find_repo(start: Path | None = None) -> Path:
    start = Path(start or Path.cwd())
    try:
        out = subprocess.run(["git", "-C", str(start), "rev-parse", "--show-toplevel"],
                             capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired) as e:
        raise ConfigError(f"git not available: {e}") from e
    if out.returncode != 0:
        raise ConfigError(f"not inside a git repository: {start}")
    return Path(out.stdout.strip())


_KEYS = {"allow_roots", "private_paths", "protected_processes", "max_compute_concurrency", "max_load_per_cpu",
         "min_free_mem_gb", "min_free_disk_gb", "read_only_requires_approval", "state_dir"}


def load(repo: Path | None = None, state_dir: Path | None = None) -> Config:
    repo = Path(repo) if repo else find_repo()
    exp = lambda p: Path(os.path.expanduser(str(p)))  # noqa: E731
    cfg = Config(repo=repo,
                 state_dir=exp(state_dir or os.environ.get("RESEARCHOS_HOME") or repo / ".researchos"),
                 allow_roots=[repo], private_paths=[exp(p) for p in DEFAULT_PRIVATE_PATHS])
    f = repo / "researchos.json"
    if f.exists():
        data = json.loads(f.read_text())
        bad = set(data) - _KEYS
        if bad:
            raise ConfigError(f"unknown keys in researchos.json: {sorted(bad)}")
        for k, v in data.items():
            if k in ("allow_roots", "private_paths"):
                v = [exp(p) if Path(os.path.expanduser(p)).is_absolute() else repo / p for p in v]
            if k == "state_dir" and state_dir is None and not os.environ.get("RESEARCHOS_HOME"):
                v = exp(v)
            elif k == "state_dir":
                continue
            setattr(cfg, k, v)
    # private paths win over any allowlist entry
    for p in cfg.private_paths:
        for r in cfg.allow_roots:
            rp, pp = Path(os.path.realpath(r)), Path(os.path.realpath(p))
            if rp == pp or pp in rp.parents:
                raise ConfigError(f"allow root {r} contains or equals private path {p}")
    return cfg
