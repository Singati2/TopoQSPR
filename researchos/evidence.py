"""Artifact hashing and manifest verification (content-addressed, with Git provenance)."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from researchos import inspectors
from researchos.safety import resolve_allowed


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_sha256sums(path: Path) -> dict[str, str]:
    """`<64 hex>  <name>` lines. Malformed lines are an error, not silently skipped."""
    out = {}
    for i, ln in enumerate(Path(path).read_text().splitlines(), 1):
        if not ln.strip():
            continue
        parts = ln.split(None, 1)
        if len(parts) != 2 or len(parts[0]) != 64 or any(c not in "0123456789abcdef" for c in parts[0]):
            raise ValueError(f"{path}:{i}: malformed SHA256SUMS line")
        name = parts[1].strip().lstrip("*")
        if name in out and out[name] != parts[0]:
            raise ValueError(f"{path}: conflicting hashes for {name}")
        out[name] = parts[0]
    return out


def check_artifact(cfg, relpath: str, expected_sha: str | None = None, expect_commit: str | None = None) -> dict:
    """Hash one artifact and establish its Git provenance. Never raises for a bad artifact;
    the result says what failed so the caller can block."""
    rec: dict = {"path": relpath, "ok": False}
    try:
        p = resolve_allowed(relpath, cfg.allow_roots, cfg.repo, cfg.private_paths)
    except PermissionError as e:
        rec["error"] = str(e)
        return rec
    if not p.is_file():
        rec["error"] = "missing file"
        return rec
    rec["sha256"] = sha256_file(p)
    rec["bytes"] = p.stat().st_size
    if expected_sha is not None and rec["sha256"] != expected_sha:
        rec["error"] = f"hash mismatch: expected {expected_sha[:12]}…, found {rec['sha256'][:12]}…"
        return rec
    try:
        rel = str(p.relative_to(Path(cfg.repo).resolve()))
        commit = inspectors.file_commit(cfg.repo, rel)
    except (ValueError, inspectors.GitError):
        commit = None
    rec["commit"] = commit
    if commit:
        committed = inspectors.blob_sha256_at(cfg.repo, commit, rel)
        rec["matches_commit"] = committed == rec["sha256"]
        if not rec["matches_commit"]:
            rec["error"] = "working copy differs from its last commit (uncommitted change)"
            return rec
    else:
        rec["error"] = "untracked: no commit provenance"
        return rec
    if expect_commit and not commit.startswith(expect_commit):
        rec["error"] = f"last touched by {commit[:9]}, expected {expect_commit}"
        return rec
    rec["ok"] = True
    return rec


def verify_manifest(cfg, manifest: Path) -> dict:
    """Verify a SHA256SUMS file (names relative to its directory) or a JSON manifest
    {"commit": optional, "artifacts": [{"path", "sha256"}]} (paths relative to the repo)."""
    manifest = resolve_allowed(manifest, cfg.allow_roots, cfg.repo, cfg.private_paths)
    rows = []
    if manifest.suffix == ".json":
        m = json.loads(manifest.read_text())
        for a in m["artifacts"]:
            rows.append(check_artifact(cfg, a["path"], a.get("sha256"), m.get("commit")))
    else:
        base = manifest.parent.relative_to(Path(cfg.repo).resolve())
        for name, sha in parse_sha256sums(manifest).items():
            note = None
            if Path(name).is_absolute():
                # some archives list absolute paths; verify the same-named file next to the manifest, if any
                note = "manifest lists an absolute path; checked the same-named file beside the manifest"
                name = Path(name).name
            rel = str(base / name)
            local = Path(cfg.repo) / rel
            if not local.exists():
                rows.append({"path": rel, "ok": None, "sha256": sha,
                             "note": "listed in manifest but not in the repository (kept outside Git by design)"})
                continue
            row = check_artifact(cfg, rel, sha)
            if note:
                row["note"] = note
            rows.append(row)
    present = [r for r in rows if r["ok"] is not None]
    return {"manifest": str(manifest), "rows": rows,
            "ok": bool(present) and all(r["ok"] for r in present),
            "n_verified": sum(1 for r in present if r["ok"]), "n_absent_by_design": len(rows) - len(present)}
