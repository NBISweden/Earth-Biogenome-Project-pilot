# 0003: Wide records and per-stage pipelines

Status: accepted. Replaces the "Wide outputs" and "Stage input" decisions of
[0002](0002-typed-records-and-stage-pipelines.md), which deferred them.

## Context

0002 converted every process and subworkflow to records but kept the stage
interfaces as they were: several channels in (reads, k-mer databases,
assemblies, coverage) and one channel per tool out. Each stage re-joins its
inputs by `id`, and `main.nf` calls the evaluation subworkflows directly.
Pipeline composition passes one record of params to an included pipeline, so a
stage with five inputs needs five params and the caller has to do the joins.

## Decision

- **One wide record per subject, one input per stage.** A stage takes a single
  `Channel` of wide records and returns the same records with its outputs added
  (`record + record(...)`). The subject is a sample for the read stages and an
  assembly build for the assembly stages. Joins happen once, where the
  subjects are created, not inside each stage.
- **Stage pipelines.** Each stage is a pipeline under `pipelines/` with a typed
  `params` block (its input channel plus its tool options) and an `output`
  block. `main.nf` includes them as `workflow as <ALIAS>` and chains them.
  Stages stay runnable on their own from a list of records in YAML.
- **Shared record types** live in `types/records.nf` and are imported with
  `include { Assembly } from ...`. A type declared in each script makes
  `Channel<Assembly>` from two scripts incompatible.
- **Assembly QC is a pipeline, not part of `main.nf`.** `pipelines/assembly_qc`
  is included by each stage pipeline that makes assemblies (assemble,
  decontaminate, purge, polish, scaffold) under one alias per stage, and by `curate`
  for the curated assemblies that were supplied (`ASSEMBLY_QC_CURATED`), for
  example `ASSEMBLY_QC_PURGED`. Assemblies supplied as input are evaluated by the
  stage that matches their `stage`. The k-mer databases travel in the wide
  record: task hashes include the alias, so a database built inside each
  included copy is rebuilt per stage.
- **Only `main.nf` publishes** when stages are included. Included pipelines
  return channels and the entry workflow decides what to publish.
- **Published tree, one folder per subject.** Read stages publish
  `<stage dir>/<sample>/<tool>/`. Assembly stages publish
  `<stage dir>/<build>/<tool>/`, with the QC of that assembly in
  `<stage dir>/<build>/qc/<tool>/`. This is a breaking change to every
  published path and is listed in the migration guide.

### Alternatives rejected

**Keep tool-first folders and move only the QC outputs.** Smaller break for
downstream tooling, but it keeps one `output` entry per tool and no output
grouped by subject, which is what wide records provide.

**QC inside `main.nf` with a stage loop.** Keeps one copy of the QC wiring but
leaves evaluation in the root workflow, so a stage cannot be run with its QC
on its own.

**Building the k-mer databases inside `assembly_qc`.** Standalone use is
simpler, but each alias rebuilds the databases (hours on a 17 Gb HiFi set). The
databases are part of the wide record, and `assembly_qc` builds them only for
assemblies that arrive without them.

## Consequences

- QUAST compares the assemblies of all stages of a sample in one run, to show how the contiguity
  changes from stage to stage. It is not part of the per-stage QC. `pipelines/compare_stages` runs
  it from `main.nf` with the assemblies of every stage, so it is the one evaluation that waits for the
  last stage. It is published in `10_report/quast/<sample>`, as before.
- Evaluation starts as soon as an assembly exists. **This drops a deliberate
  behaviour**: the old `main.nf` delayed all assembly QC until every stage had
  finished (`toList().flatMap()`), so that processing had priority on the
  cluster. The maintainers accept that QC overlaps with processing. A gate would
  have to come from a later stage than the pipeline that includes the QC, and
  Nextflow 26.10 has no documented way to express that.
- Stage pipelines depend on the experimental composition feature
  (Nextflow >= 26.10), as in 0002.
- Selectors for process config carry the alias, e.g. `PURGE_DUPLICATES:ASSEMBLY_QC_PURGED:EVALUATE_ASSEMBLY:BUSCO`.
