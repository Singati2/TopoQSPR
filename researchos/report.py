"""Markdown research decision for one task: sources, hashes, commit, findings, limitations, approvals."""
from __future__ import annotations

import datetime as dt
import json

from researchos import __version__, inspectors
from researchos import executor as X
from researchos.safety import redact, scan_untrusted

LIMITATIONS = [
    "ResearchOS re-derives and cross-checks archived numbers; it does not re-run models. "
    "`consistency_check` evidence shows archived files agree with each other, not that the analysis is right.",
    "`independent_recomputation` here refers to an earlier separate script over predictions held outside Git; "
    "ResearchOS verified the agreement, not the predictions themselves.",
    "The approval gate stops accidental or model-originated execution inside ResearchOS. It is not a security "
    "boundary against a local user with shell access.",
    "No language model was consulted unless an evidence row says `qualitative_llm_review`.",
]


def _fmt(x, nd=4):
    return f"{x:+.{nd}f}" if isinstance(x, (int, float)) else str(x)


def task_report(store, cfg, task_id: str) -> str:
    t = store.task(task_id)
    spec = t["spec"]
    rs = {r["step_id"]: r for r in store.step_results(task_id)}
    st = inspectors.repo_state(cfg.repo)
    dry = X.dry_run(store, cfg, task_id)
    L = []
    w = L.append
    w(f"# ResearchOS decision record: {spec['title']}")
    w("")
    w(f"- Task: `{task_id}` · state **{t['state']}** · attempt {t['attempt']} · spec hash `{dry['spec_hash'][:16]}`")
    tree = "clean" if st["clean"] else f"DIRTY ({st['dirty_paths']} paths)"
    w(f"- Repository: `{st['repo']}` · branch `{st['branch']}` · HEAD `{st['head']}` · working tree {tree}")
    w(f"- Generated {dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')} by ResearchOS {__version__} "
      "(deterministic; no language model, no training)")
    w(f"- Approved scopes: {', '.join(dry['approved_scopes']) or 'none'}")
    w("")
    w("## Question")
    w("")
    w(spec["question"])
    w("")
    w("## Steps")
    w("")
    w("| step | adapter | scope | status | input hash |")
    w("|---|---|---|---|---|")
    for p in dry["plan"]:
        r = rs.get(p["step"])
        status = r["status"] if r else ("awaiting approval" if any("approval" in b for b in p["blockers"]) else "not run")
        ih = r["input_hash"][:12] if r else "-"
        w(f"| {p['step']} | {p['adapter']} | {p['scope']} | {status} | `{ih}` |")
    w("")
    w("## Evidence inputs (sha256 and commit provenance)")
    w("")
    w("| path | sha256 | last commit | verified |")
    w("|---|---|---|---|")
    for r in rs.values():
        for a in (r["result"].get("data") or {}).get("artifacts", []):
            w(f"| `{a['path']}` | `{a.get('sha256', '-')[:16]}…` | `{(a.get('commit') or '-')[:9]}` | "
              f"{'yes' if a['ok'] else 'NO: ' + str(a.get('error'))} |")
    w("")
    w("## Findings")
    for sid, r in rs.items():
        w("")
        w(f"**{sid}** ({r['status']})")
        for f in r["result"].get("findings", []):
            w(f"- {f}")
        comps = (r["result"].get("data") or {}).get("comparisons")
        if comps:
            d = r["result"]["data"]
            w("")
            w(f"| context | R² {d['set_a']} | R² {d['set_b']} | gap | learner {d['set_a']} | learner {d['set_b']} | fit |")
            w("|---|---|---|---|---|---|---|")
            for c in comps:
                w(f"| {c['context']} | {c['r2_a']:.4f} | {c['r2_b']:.4f} | {_fmt(c['gap_a_minus_b'], 3)} | "
                  f"{json.dumps(c['family_a'])} | {json.dumps(c['family_b'])} | {c['fit']} |")
    for cid in spec.get("claims", []):
        c = store.claim(cid)
        if not c:
            continue
        w("")
        w(f"## Claim `{cid}`: status **{c['status']}** ({c['stage']})")
        w("")
        w(f"> {c['text']}")
        w("")
        w(f"Scope: {c['scope']}")
        if c["limitations"]:
            w("")
            w("Limitations: " + "; ".join(c["limitations"]))
        w("")
        w("| evidence | independence | numerical? | supports | method | summary |")
        w("|---|---|---|---|---|---|")
        from researchos.schema import Independence
        for e in c["evidence"]:
            num = "yes" if Independence(e["independence"]).is_numerical else "no (opinion)"
            w(f"| `{e['evidence_id']}` | {e['independence']} | {num} | {e['supports']} | {e['method']} | {e['summary']} |")
    for r in rs.values():
        nw = (r["result"].get("data") or {}).get("narrowed_wording")
        if nw:
            w("")
            w("## Proposed narrowed wording (owner decides)")
            w("")
            w(f"> {nw}")
    pending = [p for p in dry["plan"] if not p["would_run"]]
    if pending:
        w("")
        w("## Awaiting the owner (nothing below has been executed)")
        for p in pending:
            w("")
            w(f"**{p['step']}** ({p['scope']}): " + "; ".join(p["blockers"]))
            if "compute" in p:
                c = p["compute"]
                w("")
                w("Protocol draft (not approved, not frozen):")
                w("")
                for ln in str(c["protocol_draft"]).splitlines():
                    w(f"> {ln}")
                w("")
                w(f"Resource estimate: {c['resource_estimate']}")
                w("")
                w("Input checks: " + "; ".join(f"`{i['path']}` {'ok' if i['ok'] else 'NOT OK: ' + str(i.get('error'))}"
                                              for i in c["input_checks"]))
    flags = scan_untrusted(t["untrusted_notes"])
    w("")
    w("## Untrusted content")
    w("")
    if flags:
        w(f"The task file's prose contained {len(flags)} instruction-like passage(s). They were stored as text and "
          "had no effect:")
        for f in flags:
            w(f"- {f['label']}: `{f['excerpt'].strip()}`")
    else:
        w("No instruction-like text found in the task file's prose.")
    w("")
    w("## Limitations of this record")
    w("")
    for l_ in LIMITATIONS:
        w(f"- {l_}")
    w("")
    w("## Audit trail (latest first)")
    w("")
    for a in store.audit_log(25):
        if task_id in a["target"] or a["target"] in spec.get("claims", []):
            w(f"- `{dt.datetime.fromtimestamp(a['ts'], dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}` "
              f"{a['actor']} {a['action']} {a['target']} {a['detail']}")
    return redact("\n".join(L) + "\n")
