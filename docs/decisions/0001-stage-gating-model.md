# 0001: Stage gating via flat `enable_*` params + a `requires` map

Status: accepted

## Context

`main.nf` selects stages with `params.steps`, a comma-separated string
tokenized into a flat membership set and checked with `'stage' in workflow_steps`.
The check and its surrounding boilerplate repeat at each gated stage.

There is no dependency model. Prerequisites (Hi-C to CRAM conversion, FastK and
Meryl database builds) run on every invocation, and a requested subset of
stages is not checked for the inputs it needs (#329).

## Decision

- User-facing parameters are flat, individually named booleans, following
  nf-core's `skip_<name>`/`run_<name>` convention: `enable_inspect`,
  `enable_preprocess`, `enable_assemble`, `enable_screen`, `enable_purge`,
  `enable_polish`, `enable_scaffold`, `enable_curate`, `enable_alignRNA`, and a
  new `enable_prep_data` for the CRAM conversion and FastK/Meryl database builds.
- The workflow collapses them into one internal `enable` map to thread through
  subworkflow calls. The CLI and `params.yml` shape does not change.
  ```nextflow
  def enable = params
      .findAll { k, _v -> k.startsWith('enable_') }
      .collectEntries { k, v -> [(k - 'enable_'): v] }
  ```
- A `requires` map holds computed flags only. Users never set them.
  ```nextflow
  def requires = [
      prep_data: enable.inspect || enable.assemble || enable.purge || enable.polish || ...
  ]
  ```
  The prep-data gate becomes `if (enable.prep_data || requires.prep_data) { ... }`.
- Each stage keeps its current gate shape, `if (enable.<stage>) { ... } else { ... }`.
  The existing if/else and variable reassignment already route a disabled
  stage's consumer to the nearest enabled producer. A stage registry with
  `consumes`/`produces` edges and a resolver is not needed for this.

### Alternatives rejected

**A nested `params.stages` map validated against a declared dependency graph**
(`dependsOn`/`consumes`/`produces` per stage, resolved by one generic function).
This was the first design. It solves a rerouting problem the if/else shape
already solves, and adds registry and resolver code for no behavioural gain at
this size. Reconsider it if `requires` outgrows a flat block of `&&`/`||`
expressions.

**Target-driven gating** (`target: scaffold`, resolving dependencies backward).
- It cannot express skipping an optional stage mid-chain without a second
  override mechanism.
- It cannot validate always-on side-branch stages (`evaluate`, `report`), whose
  inputs are not on the main lineage.
- It does not help resume across clusters or runs, which is structural under
  any gating model.

## Consequences

- Adding a stage or a derived prerequisite means one `enable_X` param, plus a
  `requires.X` line if it is derived. The repeated four-part boilerplate goes away.
- The CLI and `params.yml` stay flat (`--enable_purge false`).
- There is no registry or DSL to maintain until `requires` outgrows a flat block
  of boolean expressions.
- Typed params (#413) do not force a change here: `params.findAll` works with a
  typed `params {}` block. biovat#132 instead builds `enable` and `requires` as
  records in `main.nf`; either form keeps the user-facing flags unchanged.
