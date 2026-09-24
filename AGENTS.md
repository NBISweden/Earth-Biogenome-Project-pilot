# Working with AI agents on this repository

Process conventions for AI coding agents (Claude Code or similar) working in
this repository. For DSL2/Groovy code style, see
[Nextflow code conventions](#nextflow-code-conventions).

## Issue first, implementation later

Multiple species-assembly repos across several institutions run against this
pipeline. Don't implement a change straight from an idea, a design doc, or a
"this needs fixing" observation. Instead:

1. Draft the exact issue title and body.
2. Show it to the maintainer for review.
3. File it (`gh issue create`) once approved.
4. Implement only when explicitly asked, as a separate step.

This applies to small, clearly correct changes too (e.g. a one-line gating
refactor). The issue provides the paper trail and a review point before code
lands.

## All changes go through a PR

`main` is branch-protected. Never push to it directly, even with an admin
bypass - including release commits (version bumps, CHANGELOG entries, tags).

- **Squash merge.** `main` has one commit per PR.
- **Amend, don't stack.** While iterating on a PR, amend existing commits
  instead of adding "fix" commits on top.
- **CHANGELOG entries cite the PR number.** Open the PR first, then amend the
  CHANGELOG entry to add the number (e.g. `... (#123)`) as the last change
  before merge.

## Design decisions go in ADRs

Record nontrivial design choices (a rejected alternative, a chosen tradeoff)
as numbered ADRs in `docs/decisions/NNNN-short-title.md` (create the directory
for the first one). Propose an ADR rather than leaving the decision only in a
PR description, issue comment, or standalone planning doc.

## Multi-part work uses an Epic with sub-issues

For work that splits into independent pieces, file one tracking issue (Issue
Type: Epic) and link each piece as a native GitHub sub-issue, not a checklist
of issue numbers. Link already-filed issues instead of duplicating them.

## Anonymize species in evidence from other BGE/ERGA repos

When evidence comes from other BGE species-assembly repos or the
ERGA-consortium/EARs report corpus, report by institution/team and count
("3 reports from institution X"), never by species name. Those assemblies are
other teams' unpublished work.

## Nextflow code conventions

Before writing or reviewing any `.nf` file, load the
`nextflow-implementation-patterns` skill if available. Otherwise, follow
[nf-core module conventions](https://nf-co.re/docs/guidelines/components/modules)
and match the existing code in `modules/` and `subworkflows/`.

## Environment and validation

- The repo is `pixi`-managed. Run tools via `pixi run <task>`; don't assume
  `nextflow` or `nf-test` are on `PATH`.
- Minimum validation for a Nextflow change is a stub run:
  `pixi run nextflow run main.nf -profile test -stub-run`. Full validation is
  `pixi run nftest-docker` (or `nftest-singularity`).
- Do exploratory or bulk data work (cloning other repos, extracting data from
  many files) in a scratch directory outside this checkout, with any one-off
  CLI tools installed via `pixi` there.
