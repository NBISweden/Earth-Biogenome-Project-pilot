# Changelog

All notable changes to this pipeline will be documented in this file, going
forward from this release. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Changed (breaking)

- Requires Nextflow >= 26.10.0. The workflow is statically typed: every
  process, local and nf-core, takes and emits one record, and every
  subworkflow and `main.nf` use records instead of tuples. See the
  [migration guide](docs/migration-typed-records.md) and ADR 0002.
- `when:` and `ext.when` are removed. Stages are switched with `steps`.
- `ext.args` and `ext.prefix` for the workflow's own tool calls are record
  fields set by the workflows, not selectors in `configs/modules.config`.
  A user selector still overrides them. Process names carry the stage alias, so
  selectors with stage names that no longer exist match nothing.
- Results are published with `publish:` and an `output {}` block instead of
  `publishDir`. The published folders are unchanged, except that the HiFi
  coverage tracks, which matched no file before, are no longer published.
- Pipeline params are declared in `main.nf`. `nextflow.config` keeps only
  config params (`outdir`, `publish_mode`, `cache`). The nested `stage { ... }`
  config scope is now the param `stage_dirs`.
- ENA is queried in the workflow, and `TOL_SEARCH` writes `DToL.tsv`
  (`ENA_TAXQUERY` and `REPORT_DTOL` are removed).
- 47 nf-core modules are patched with `nf-core modules patch` and tracked in
  `modules.json`.
- The upstream tests of vendored nf-core modules (`modules/nf-core/**/tests/`)
  are removed from git and gitignored. They are tested upstream.

### Added

- `pipelines/assembly_qc`: evaluate assemblies (MerquryFK, Merqury, BUSCO,
  gfastats, QUAST) as a pipeline of its own.
- `stub:` blocks for 25 processes. `-stub-run` now runs the whole workflow.
- ADR 0001 (stage gating) and ADR 0002 (typed records and stage pipelines).

### Fixed

- `genescopefk` declared `prefix` twice and did not compile under Nextflow
  >= 26.04.
- `SMUDGEPLOT` is joined to `GENESCOPEFK` by sample id, not by channel order.
- hifiasm ultra-long reads are passed comma-separated.

## [1.0.0] - 2026-09-23

First tagged release. This formalizes the pipeline's current state rather
than marking a fresh rewrite: NBIS/SciLifeLab already runs this pipeline in
production for de novo genome assembly and curation handoff, ahead of it
having a version number attached.

### Included

- Staged assembly workflow driven by `params.steps`: read inspection,
  assembly (HiFi, with optional Hi-C-guided hifiasm and organelle assembly),
  contamination screening, duplicate purging, error polishing, scaffolding,
  and rapid curation handoff (Pretext/CurationPretext).
- RNAseq alignment as an optional side branch.
- Assembly evaluation (BUSCO, Merqury/MerquryFK, k-mer spectra) and a
  Quarto-based assembly report, with MultiQC-collected logs.
- Resume-friendly design: any stage's output can be supplied directly as
  input, letting a run start partway through.

### Known limitations, tracked for post-1.0 work

- Stage gating is a flat, unordered `params.steps` membership list with
  repeated per-stage boilerplate, and some prerequisites (Hi-C to CRAM
  conversion, FastK/Meryl database builds) run unconditionally regardless of
  which stages are requested - see the architecture refactor epic (#381) and
  ADR 0001 (#378) for the planned replacement.
- The assembly data model doesn't yet cleanly support multiple
  assemblies/haplotypes per sample (#246, #338) - a record-type redesign is
  proposed on #338.
- No CI is currently configured for this repository.

See the [architecture refactor epic](https://github.com/NBISweden/Earth-Biogenome-Project-pilot/issues/381)
for the tracked plan going forward.
