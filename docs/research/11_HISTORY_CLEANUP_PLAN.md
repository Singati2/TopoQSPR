# 11 — Git history cleanup plan for the copyrighted PDF (FOR APPROVAL; not executed)

**Status 2026-10-13:** the file
`related_papers/Arockiaraj_2023_SARQSAR_hybrid_GHZ_descriptors_photos.pdf` (blob
`e5382ef118e82b7c7d6fff959d43a3de066063d1`, 2,972,693 bytes, sha256 00e804ae…) is removed
from the **current tree** of all five branches by ordinary commits: `main` 3ef0b50,
`v1.0-dev` 11bc1e3, `audit-tests-20261012` 8aad4d7, `audit-2026-10-12` 5dfdaa4,
`research-2026-10-13` 8b9b2e4. It remains in **history**: reachable from every commit since
`ec19cb1` (2026-10-03), the tag `plan-freeze-v1`, the PR refs `refs/pull/1/head` and
`refs/pull/2/head` (read-only on GitHub), any fork or clone, and GitHub's caches. The owner's
research copy is at `~/topoqspr_private/related_papers/` (outside every repository) and must
not be touched by this plan.

## Option A — do nothing further (accept residual exposure)
Cheapest. The file stays downloadable through old commit URLs. Not recommended for a public repo.

## Option B — make the repository private (reversible, no history change)
Stops all public access, including history and PR refs, immediately. It can be made public
again after Option C. Recommended as the **first** step if a full purge is wanted.

## Option C — full history purge (destructive; needs explicit approval of each step)

1. Freeze: announce to co-authors that no one pushes. Wait for E3 to finish (its run record
   cites `fc69883`, which will change).
2. Back up: `git clone --mirror https://github.com/Singati2/TopoQSPR.git topoqspr-backup.git`
   plus a `git bundle create … --all`. Store both privately (they contain the PDF).
3. Rewrite on a fresh mirror clone:
   `git filter-repo --path related_papers/Arockiaraj_2023_SARQSAR_hybrid_GHZ_descriptors_photos.pdf --invert-paths`.
   Check `git log --all -- <path>` is empty and that every tree is otherwise unchanged
   (`git diff` old ↔ new tip of each branch must be empty).
4. Save `.git/filter-repo/commit-map` as `docs/sha_map_<date>.tsv` (second map after the
   2026-10-12 trailer rewrite). Chain it with `sha_map_2026-10-12.tsv` so every run record
   (`659fac3…`, `fe92837…`, `fc69883…`, …) resolves to the final IDs.
5. Force-push all five branches and the tag: `git push --force --mirror` is **not** used
   (it would delete branches). Push the named refs only. Re-create `plan-freeze-v1` on the
   mapped commit and record the change in `docs/EXPERIMENT_PLAN.md`.
6. Ask GitHub Support to (a) delete cached views of the blob, (b) remove the PR refs or
   rewrite them, (c) run garbage collection. Check whether forks exist (`gh api
   repos/Singati2/TopoQSPR/forks`).
7. All collaborators re-clone. Old clones still contain the file.

**Side effects:** every commit ID after 2026-10-03 changes again. Run records, paper macros
(`\runCommit`) and audit documents that cite IDs need the map. PR #2's history changes. Do it
**once**, before any release tag or Zenodo DOI.

## Recommendation
B now (if acceptable), then C once E3 has finished and co-authors agree, in a single
maintenance window that also covers any other pre-release history issues.
