"""Path allowlist, untrusted-text screening and report redaction."""
from __future__ import annotations

import os
import re
from pathlib import Path


class PathBlocked(PermissionError):
    pass


def resolve_allowed(path: str | Path, roots: list[Path], base: Path, deny: list[Path] = ()) -> Path:
    """Resolve `path` (relative to `base`) and require that it, after following every symlink,
    lies inside one of `roots` and inside none of `deny`. A symlink inside the repo that points
    outside is blocked."""
    p = Path(os.path.expanduser(str(path)))
    if not p.is_absolute():
        p = base / p
    real = Path(os.path.realpath(p))
    for d in deny:
        dd = Path(os.path.realpath(d))
        if real == dd or dd in real.parents:
            raise PathBlocked(f"path is in a private location: {path}")
    for r in roots:
        rr = Path(os.path.realpath(r))
        if real == rr or rr in real.parents:
            return real
    raise PathBlocked(f"path outside allowlist: {path} -> {real}")


# Command-like or instruction-like content inside documents, logs or task notes. It is never executed;
# it is flagged so the report can say that untrusted text tried to direct actions.
_SUSPICIOUS = [
    (r"\brm\s+-[a-z]*r[a-z]*f?\b", "recursive delete"),
    (r"git\s+push\s+.*(--force|-f\b|\+)", "force push"),
    (r"git\s+(reset\s+--hard|filter-(repo|branch)|branch\s+-D|push\s+--delete)", "destructive git"),
    (r"(curl|wget)[^|\n]*\|\s*(ba)?sh", "pipe-to-shell"),
    (r"\bsudo\b", "privilege escalation"),
    (r"\b(kill|pkill|killall)\b", "process kill"),
    (r"ignore (all |the )?(previous|prior|above) (instructions|rules)", "prompt injection"),
    (r"(you are|act as) (now )?(the )?owner|owner (has )?approved", "approval impersonation"),
    (r"chmod\s+777|>\s*/dev/sd", "dangerous system change"),
]


def scan_untrusted(text: str) -> list[dict]:
    hits = []
    for pat, label in _SUSPICIOUS:
        for m in re.finditer(pat, text or "", flags=re.IGNORECASE):
            hits.append({"label": label, "excerpt": text[max(0, m.start() - 20): m.end() + 20].replace("\n", " ")})
    return hits


_SECRET_ENV = re.compile(r"(TOKEN|SECRET|PASSWORD|PASSWD|API_KEY|APIKEY|ACCESS_KEY|PRIVATE_KEY|CREDENTIAL)", re.I)
_EMAIL = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
_TOKENLIKE = re.compile(r"\b(gh[pousr]_[A-Za-z0-9]{20,}|sk-[A-Za-z0-9_-]{20,}|xox[abpr]-[A-Za-z0-9-]{10,})\b")


def redact(text: str, env: dict | None = None) -> str:
    env = os.environ if env is None else env
    out = text
    for k, v in env.items():
        if v and len(v) >= 6 and _SECRET_ENV.search(k):
            out = out.replace(v, f"<redacted:{k}>")
    home = env.get("HOME")
    if home and len(home) > 1:
        out = out.replace(home, "~")
    out = _TOKENLIKE.sub("<redacted:token>", out)
    out = _EMAIL.sub("<redacted:email>", out)
    return out
