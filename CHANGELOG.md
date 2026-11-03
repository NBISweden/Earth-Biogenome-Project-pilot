# Changelog

All notable changes to this pipeline will be documented in this file, going
forward from this release. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Changed (breaking)

- Requires Nextflow >= 26.10.3 (#428). The workflow is statically typed: every
  process, local and nf-core, takes and emits one record, and every
  subworkflow and `main.nf` use records instead of tuples. See the
  [migration guide](docs/migration-typed-records.md) and ADR 0002.
- `when:` and `ext.when` are removed. Stages are switched with `steps`.
- Nested params (`fastk`, `meryl`, `genescopefk`, `busco`, `fcs`) are records. A
  value in a params file replaces the whole record, so `fastk: {kmer_size: 21}`
  drops the other fields of `fastk`. The unused `diamond.db` and `blast.db`
  params are removed.
- Tools write temporary files to the task work directory (`env.TMPDIR = '$PWD'`
  in `nextflow.config`), not to the node's `TMPDIR`. Override it with the `env`
  scope of your own config.
- `stage {}` is the `stage_dirs` param; a partial map keeps the default folders of the other stages.
- The stage pipelines take the QC options as one `qc` param. The `alt_gfa` key of
  the input assemblies is gone.
- `--input` takes one YAML file, not a glob. `chunk_size` is a memory unit
  (`1.5 GB`, no quotes). With an Arima `hic_type`, hifiasm reads the trimmed
  Hi-C pairs that the later stages use.
- Hi-C reads are optional. A run without them stops at the first step that
  needs them (`scaffold`, `curate`, `hifiasm_with_hic`) with an error at the start.
- Resources of `HIFIASM`, `BUSCO`, `MERQURY`, `MERQURYFK`, `GFASTATS`, `FASTK_FASTK`, `MERYL_COUNT`
  and `MERYL_UNIONSUM` are measured on a real genome and scale with the genome size. Killed tasks
  are retried twice with more memory and time. The trace records requested and used resources
  (`docs/development.md`, "Resources").
- Resource and option selectors are plain process names. A selector on the original name
  (`MERQURY_MERQURY`) also matches its aliases, as before; one on an alias (`ASSEMBLY_QC_PURGED`)
  matches only that instance.
- `ext.args` and `ext.prefix` for the workflow's own tool calls are record
  fields set by the workflows, not selectors in `configs/modules.config`.
  A user selector still overrides them. Process names carry the stage alias, so
  selectors with stage names that no longer exist match nothing.
- Results are published with `publish:` and an `output {}` block instead of
  `publishDir`. **The published tree changed:** each stage folder has one
  subfolder per sample or assembly build with the tools inside, and the QC of
  an assembly is in `<stage dir>/<build>/qc/`. The moved paths are listed in
  the migration guide. The HiFi coverage tracks, which matched no file before,
  are no longer published.
- Stages are pipelines under `pipelines/` that take one wide record per sample
  or assembly (`types/records.nf`). Assembly QC runs per stage: evaluation
  starts when an assembly exists. The old workflow held QC back until all stages
  had finished, so that processing had priority; QC now overlaps with
  processing (accepted, ADR 0003). `ASSEMBLE` waits for the HiFi k-mer databases. `PURGE_DUPLICATES` no
  longer waits for the k-mer coverage estimate. See ADR 0003.
- HiFi reads (`hifi`) are required in the input.
- Pipeline params are declared in `main.nf`. `nextflow.config` keeps only
  config params (`outdir`, `publish_mode`, `cache`). The nested `stage { ... }`
  config scope is now the param `stage_dirs`.
- `TOL_SEARCH` writes `DToL.tsv` (`REPORT_DTOL` is removed).
- 47 nf-core modules are patched with `nf-core modules patch` and tracked in
  `modules.json`.
- The upstream tests of vendored nf-core modules (`modules/nf-core/**/tests/`)
  are removed from git and gitignored. They are tested upstream.

### Added

- `pipelines/assembly_qc`: evaluate assemblies (MerquryFK, Merqury, BUSCO,
  gfastats) as a pipeline of its own, included by the stage pipelines.
- `pipelines/compare_stages`: QUAST over the assemblies of all stages, one run per
  sample, as before (`10_report/quast`).
- Stub-based nf-tests for the local modules (`pixi run nftest-modules`).
- `stub:` blocks for 23 more processes (24 of the 29 local processes and all 47 nf-core ones have one). `-stub-run` now runs the whole workflow.
- ADR 0001 (stage gating), ADR 0002 (typed records and stage pipelines) and
  ADR 0003 (wide records and per-stage pipelines).

### Fixed

- The typed conversion of `merqury` and `merquryfk` first declared one file for
  outputs that they write once per assembly (`*_only.bed`, `*_only.wig`, `*.qv`
  and the `*.spectra-cn.*` files). The real-tool tests found it. They return
  all of them, as the legacy `path()` did.
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
