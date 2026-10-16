"""Read-only inspectors: Git state, processes, machine health, tool availability.

Git is called only with fixed, read-only argument lists. Processes are discovered from /proc and
never signalled, traced or attached to.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

READ_ONLY_GIT = {
    ("rev-parse",), ("status",), ("branch",), ("log",), ("show",), ("ls-remote",), ("worktree", "list"),
    ("cat-file",), ("ls-files",), ("remote",), ("hash-object",),
}


class GitError(RuntimeError):
    pass


def git(repo: Path, *args: str, timeout: int = 20) -> str:
    if not any(tuple(args[: len(k)]) == k for k in READ_ONLY_GIT):
        raise GitError(f"git command not in the read-only allowlist: {args[:2]}")
    if args[0] == "hash-object" and "-w" in args:
        raise GitError("hash-object -w writes to the object store")
    if args[0] == "remote" and len(args) > 1 and args[1] not in ("-v", "get-url"):
        raise GitError("only `git remote -v` / `get-url` are allowed")
    try:
        out = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, timeout=timeout,
                             env={**os.environ, "GIT_TERMINAL_PROMPT": "0", "GIT_OPTIONAL_LOCKS": "0"})
    except (OSError, subprocess.TimeoutExpired) as e:
        raise GitError(str(e)) from e
    if out.returncode != 0:
        raise GitError(out.stderr.strip() or f"git {' '.join(args)} failed")
    return out.stdout


def repo_state(repo: Path, remote: bool = False) -> dict:
    d: dict = {"repo": str(repo)}
    d["head"] = git(repo, "rev-parse", "HEAD").strip()
    d["branch"] = git(repo, "rev-parse", "--abbrev-ref", "HEAD").strip()
    status = git(repo, "status", "--porcelain=v1").splitlines()
    d["dirty_paths"] = len(status)
    d["clean"] = not status
    d["worktrees"] = [ln.split()[0] for ln in git(repo, "worktree", "list").splitlines() if ln.strip()]
    if remote:
        try:
            heads = git(repo, "ls-remote", "--heads", "origin", timeout=30)
            d["remote_heads"] = {ln.split("\trefs/heads/")[1]: ln.split()[0] for ln in heads.splitlines() if "\t" in ln}
        except GitError as e:
            d["remote_heads_error"] = str(e)
    return d


def file_commit(repo: Path, relpath: str) -> str | None:
    """Last commit that touched relpath, or None if untracked."""
    out = git(repo, "log", "-1", "--format=%H", "--", relpath).strip()
    return out or None


def blob_sha256_at(repo: Path, commit: str, relpath: str) -> str:
    import hashlib
    try:
        out = subprocess.run(["git", "-C", str(repo), "show", f"{commit}:{relpath}"], capture_output=True, timeout=20)
    except (OSError, subprocess.TimeoutExpired) as e:
        raise GitError(str(e)) from e
    if out.returncode != 0:
        raise GitError(out.stderr.decode(errors="replace").strip())
    return hashlib.sha256(out.stdout).hexdigest()


def processes(patterns: list[str]) -> list[dict]:
    """Processes whose command line contains any pattern. Reads /proc only."""
    found = []
    mine = _ancestors(os.getpid())  # this CLI and the shells that launched it are not "jobs"
    for pid in os.listdir("/proc") if os.path.isdir("/proc") else []:
        if not pid.isdigit() or int(pid) in mine:
            continue
        try:
            raw = Path(f"/proc/{pid}/cmdline").read_bytes()
        except OSError:
            continue
        cmd = raw.replace(b"\0", b" ").decode(errors="replace").strip()
        if not cmd:
            continue
        hit = [p for p in patterns if p in cmd]
        if hit:
            found.append({"pid": int(pid), "match": hit, "cmd": cmd[:200]})
    return found


def _ancestors(pid: int) -> set[int]:
    out = set()
    while pid > 1 and pid not in out:
        out.add(pid)
        try:
            stat = Path(f"/proc/{pid}/stat").read_text()
            pid = int(stat.rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            break
    return out


def health() -> dict:
    h: dict = {"cpus": os.cpu_count()}
    try:
        h["load1"], h["load5"], h["load15"] = os.getloadavg()
    except OSError:
        h["load1"] = None
    try:
        mem = {}
        for ln in Path("/proc/meminfo").read_text().splitlines():
            k, v = ln.split(":", 1)
            mem[k] = int(v.split()[0])
        h["mem_total_gb"] = round(mem["MemTotal"] / 2**20, 1)
        h["mem_available_gb"] = round(mem["MemAvailable"] / 2**20, 1)
    except (OSError, KeyError, ValueError):
        h["mem_available_gb"] = None
    du = shutil.disk_usage(Path.home())
    h["disk_free_gb"] = round(du.free / 2**30, 1)
    if shutil.which("nvidia-smi"):
        try:
            out = subprocess.run(["nvidia-smi", "--query-gpu=name,utilization.gpu,memory.used,memory.total",
                                  "--format=csv,noheader"], capture_output=True, text=True, timeout=10)
            h["gpus"] = [ln.strip() for ln in out.stdout.splitlines() if ln.strip()]
        except (OSError, subprocess.TimeoutExpired):
            h["gpus"] = "nvidia-smi failed"
    else:
        h["gpus"] = []
    return h


def tools() -> dict:
    t = {}
    for name in ("git", "gh", "claude", "sbatch", "squeue", "nvidia-smi"):
        t[name] = shutil.which(name) is not None
    return t
