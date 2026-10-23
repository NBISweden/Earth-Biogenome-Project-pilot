# 0002: Typed records and stage pipelines

Status: accepted

## Context

Each stage of the workflow should run on its own and the whole workflow should
be includable in another pipeline. Nextflow 26.10 supports this through pipeline
composition: a pipeline's `params` block acts as `take:` and its `output` block
as `emit:`. Composition needs static typing, `params` read in one place, no
project-level assets (`projectDir`, root `bin/`), and tool arguments passed as
process inputs instead of `conf/modules.config` selectors.

Today processes and subworkflows pass `(meta, file)` tuples and a nested `meta`
Map, `configs/modules.config` holds 72 selectors, and `params` is read in 12
subworkflow and module files. The record design on #338 kept modules on tuples
behind permanent adapters. biovat#132 converts every module to records instead.

## Decision

- **Records everywhere.** Every process, including nf-core modules, takes one
  record and emits one wide record. nf-core modules are converted and tracked
  with `nf-core modules patch`. This replaces the permanent adapters proposed
  on #338; adapters exist only inside a conversion commit.
- **Module spec** (from biovat#132): required `id`, optional `meta` (never read
  in process code), `<tool>_args`, `output_prefix`, `ext [ ... ]`,
  `tag "${id}"`, versions through `topic:`. No `when:` or `ext.when`.
- **Typed params.** Each stage declares its own `params {}` and reads `params`
  only in its entry workflow. Related params are bundled in one record type per
  tool. A record param does not merge field defaults (#415), so fields are
  optional and the entry workflow applies defaults. User-facing names such as
  `fastk.kmer_size` stay.
- **Wide outputs.** Every workflow joins related outputs into one record per
  sample or assembly, and the `output {}` block publishes that record. The
  published directory tree changes in the next major release.
- **Stage input.** Each stage accepts the current nested `--input` YAML and uses
  the parts it needs. Pipelines that include a stage pass typed channel params.
- **Order** (agent-skills `migrate-nextflow-code`): strict syntax, topic
  channels, static typing (processes, workflows, params), workflow outputs,
  pipeline composition. Each step is checked against the previous one.
- **Floor and release.** Nextflow >= 26.10 (composition is experimental).
  Released as v2.0.0 with a migration guide (#426).

### Alternatives rejected

**Permanent `(meta, path)` adapters (#338).** Keeps nf-core modules unpatched,
but every call site carries an adapter and composition still needs typed
channels at the stage boundary.

**Flat CSV per stage (biovat).** Cleaner types, but a new input format for every
stage and a migration for existing YAML users.

**Per-tool outputs with an unchanged tree.** Avoids a layout change, but needs
one `emit` per tool and one `output {}` entry per channel.

**v1.x with a long-lived 1.x branch.** Supports sites without Nextflow 26.10 at
the cost of maintaining two lines.

## Consequences

- 47 nf-core modules are patched; `nf-core modules update` applies each patch on
  top of the new upstream version and reports conflicts.
- The published tree, process names (`STAGE:TOOL`) and several params change.
  Downstream configs with fully qualified selectors stop matching without an
  error.
- New modules follow the spec from the day this ADR lands.
- Included pipelines need Nextflow >= 26.10, and composition may change in later
  Nextflow releases.
