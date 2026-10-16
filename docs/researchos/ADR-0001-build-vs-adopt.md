# ADR-0001: Build a small evidence-first layer; borrow ideas, not code, from FirstMate and OrchestraOS

Date: 2026-10-16 · Status: accepted for v0.1 (owner may revisit)

## Context

The owner asked for a research orchestrator "inspired by FirstMate and OrchestraOS". Before
designing anything, both names were checked against primary sources (GitHub API and README, read
2026-10-16).

| Name | What it actually is (verified) | Licence | Status |
|---|---|---|---|
| **FirstMate**: [github.com/kunchenguid/firstmate](https://github.com/kunchenguid/firstmate) | An "agent distro": a cloned repo of `AGENTS.md`, skills and helper scripts that a terminal coding agent (Claude Code, Codex, Grok, Pi, OpenCode, …) follows. One "first mate" dispatches and supervises "crewmates" in tmux/Herdr panes, each in its own Git worktree (via `treehouse`), and returns PRs or scout reports (`data/<id>/report.md`). It describes itself as "not a model, not a harness, not a skill, not an MCP server, and not a CLI". | MIT (GitHub API) | Active: created 2026-06-12, pushed 2026-10-16, ~8k stars, no release tags |
| **OrchestraOS**: the closest match is [github.com/Tulum-DAO/orchestraos](https://github.com/Tulum-DAO/orchestraos) ("an open harness for running a fleet of coding agents as a team") | A server-hosted harness: agent "seats" with generations and handoff-based context rotation, a durable inter-agent message store, a facts/memory store, an approvals surface of approve/deny "cards" on a web dashboard and iOS/watch app, `orchestra doctor`, systemd units. Its README says: "This is not a finished product." | Apache-2.0 (GitHub API) | Created 2026-09-17, 0 stars, 2 forks |
| Other "OrchestraOS" hits | `QRcode1337/OrchestraOS` (a dashboard theme, no licence), `gilldilshaan/orchestraos` (MIT, "enterprise AI operations"), several ShawCole repos (no licence), a Devpost hackathon entry (page returned 403). | various / none | Not the intended software, or unverifiable |

So **"OrchestraOS" is ambiguous**. The Tulum-DAO project is the only one whose description matches
"fleet of coding agents with human approvals", but it is three weeks old and self-described as
unfinished. Neither project is a scientific-evidence tool. Both orchestrate *coding* agents.

## Decision

Build `researchos/` as a small, standard-library-only Python package inside TopoQSPR. **No code is
copied from either project**, so neither licence imposes obligations. Neither project is installed.
Their setup scripts were not run: FirstMate's first mate offers to install tools, and OrchestraOS
installs systemd units and needs a VPS and Tailscale.

**Borrowed concepts (re-implemented):**
- FirstMate:
  - each line of work runs in an isolated Git worktree (this branch was built that way);
  - one coordinator escalates only real decisions;
  - a watcher that costs no tokens;
  - restart-proof state on disk;
  - investigation reports saved as Markdown;
  - a self-repository boundary: the orchestrator does not publish changes to its own home without consent.
- OrchestraOS:
  - every human decision is an explicit approve/deny record that resumes the waiting work;
  - a `doctor` command;
  - durable, append-only state;
  - secrets only in the environment, never in config.

**Rejected for v0.1:**
- multi-pane live agent fleets and autonomous merges (FirstMate `+yolo`);
- context rotation, phone/watch apps, Telegram/voice, systemd services and VPS deployment (OrchestraOS);
- any mode in which an agent executes LLM-originated shell strings.

These add risk and infrastructure without serving the scientific need: claims tied to verifiable
evidence, with fail-closed gates.

## Consequences

- v0.1 is deterministic and works with zero API access. Agent roles beyond coordinator,
  verifier/reviewer and steward are declared interfaces only (`researchos agents list` says so).
- A future adapter could run FirstMate as the *coding* crew for implementation tasks, while
  ResearchOS stays the evidence and approval authority. That would need its own ADR.

## Acceptance criteria for v0.1

See `IMPLEMENTATION_REPORT.md` §Acceptance. In short: `doctor` and `status` are read-only on the
real repo; one real E3 artifact is ingested with a verified checksum and commit; the worked example
produces a source-linked report without training; unapproved execution fails closed under test;
crash/resume and the negative paths are tested; the work is on its own branch only.
