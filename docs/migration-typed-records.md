# Migration to typed records and stage pipelines

This page describes the change from tuple-based processes and subworkflows to statically typed Nextflow with records, and from one workflow to pipelines that can run on their own. It is written for people who run the workflow and for contributors who know the code from before this change. The decision is recorded in [ADR 0002](decisions/0002-typed-records-and-stage-pipelines.md).

- [Why](#why)
- [What changed at a glance](#what-changed-at-a-glance)
- [If you run the workflow](#if-you-run-the-workflow)
- [Processes and modules](#processes-and-modules)
- [Data between stages](#data-between-stages)
- [The main workflow](#the-main-workflow)
- [Pipelines that run on their own](#pipelines-that-run-on-their-own)
- [Publishing](#publishing)
- [Configuration](#configuration)
- [Patched nf-core modules](#patched-nf-core-modules)
- [Behaviour changes](#behaviour-changes)
- [Tests and verification](#tests-and-verification)
- [Known limitations and follow-ups](#known-limitations-and-follow-ups)
- [How to add a module or a stage](#how-to-add-a-module-or-a-stage)

## Why

Nextflow 26.10 can include a whole pipeline (its `params` block, entry `workflow` and `output` block) from another pipeline. The `params` block acts as `take:` and the `output` block as `emit:`. That lets a stage run on its own and lets the workflow be one step in a longer analysis.

Composition needs static typing, `params` read in one place, results published through the `output` block, and tool arguments passed as process inputs instead of through `configs/modules.config` selectors. This change does that. See the Nextflow docs on [pipeline composition](https://docs.seqera.io/nextflow/workflow-typed#pipeline-composition) and [typed parameters](https://docs.seqera.io/nextflow/typed-parameters). Composition is experimental in Nextflow 26.10 and may change.

## What changed at a glance

| Before | After |
| --- | --- |
| Nextflow >= 25.10.2 | Nextflow >= 26.10.0 |
| Tuples such as `[ meta, files ]`; `meta` is a Map | Records such as `record(id: ..., meta: ..., reads: ...)`; `meta` is a Record |
| `when:` and `ext.when` on 75 processes | No `when:`. The workflow decides |
| `ext.args` and `ext.prefix` selectors in `configs/modules.config` | `<tool>_args` and `output_prefix` fields set by the workflows |
| `publishDir` selectors in `configs/modules.config` | A `publish:` section and an `output {}` block in `main.nf` |
| Params defaults in `nextflow.config`, read anywhere through `params` | Params declared and typed in `main.nf`; subworkflows get values as inputs |
| `.out` and legacy operators (`branch`, `multiMap`, `groupTuple`, `combine(by:)`) | Variables, `filter`/`map`/`join(by:)`/`groupBy` |
| Helper functions `combineByMetaKeys`, `joinByMetaKeys`, `deepMergeMaps`, ... | Removed. Joins use named record fields |
| `ENA_TAXQUERY` and `REPORT_DTOL` processes | ENA is queried in `PREPARE_INPUT`; `TOL_SEARCH` also writes `DToL.tsv` |
| `collectFile` for the versions file and report params | `WRITE_FILE` process |
| One workflow | The workflow, and `pipelines/assembly_qc` as a pipeline of its own |

## If you run the workflow

- **Install Nextflow >= 26.10.0.** The manifest stops older versions.
- **The input and the params are the same.** `--input` is still the nested YAML (`sample`, `hifi`, `hic`, `assembly`, ...). `-params-file` keys keep their names, including nested ones such as `fastk.kmer_size` and `busco.lineages`.
- **Nested params are records without field defaults.** If you give `fastk: { kmer_size: 21 }` you override the whole `fastk` record. Fields you leave out are `null`, and the workflow applies its default (31 for the k-mer sizes). The `busco` and `fcs` records work the same way.
- **The published tree is the same** (`01_read_inspection`, `03_assembly/<build>`, ...). One file set is no longer published: the HiFi coverage tracks matched a pattern (`*.hitile`) that no module produces.
- **Your `nextflow.config` selectors:**
  - `ext.args` and `ext.prefix` selectors still work and **replace** the value the workflow passes. If you set `ext.args` for a tool, include the flags the workflow relies on. A selector for `FASTK_FASTK` that sets `ext.args` replaces `-t1 -k31 -v`.
  - Process names carry the alias of the stage workflow, for example `ASSEMBLE:ASSEMBLE_HIFI:HIFIASM`. A selector with the plain name (`HIFIASM`) still matches. A selector with a stage name that no longer exists (for example `EVALUATE_PURGED_ASSEMBLY:MERQURYFK_MERQURYFK`) matches nothing and gives **no error**. Check the names in `.nextflow.log` with `nextflow log -f process`.
  - `ext.when = false` no longer switches a process off. Use `steps`, or `process.when` to drop samples (the closure sees the input record fields, for example `id`).
  - `stage { ... }` in `configs/modules.config` is now the param `stage_dirs`. `stage` is a built-in config scope.
- **`publish_mode`, `outdir` and `cache`** are config params. Set them in `-params-file` or `--outdir` as before.
- **A pipeline that includes this workflow** has to enable static typing and pass the params as a record. See [Pipelines that run on their own](#pipelines-that-run-on-their-own).

## Processes and modules

Every process, local and nf-core, follows the same rules.

1. One record as input and one record as output.
2. Versions go to the `versions` topic as `tuple(task.process, 'tool', version) >> 'versions'`.
3. The input record has a required `id` String, an optional `meta` Record, an optional `<tool>_args` String for each tool call, an optional `output_prefix` String, and then the input files and values. Optional files use `?`. Values that control the process (for example `single_end`) are fields of their own. Process code never reads `meta`.
4. The output record has `id`, `meta` and every output of the process. Optional outputs use `file(..., optional: true)`; they are `null` when missing.
5. `ext` values are independent closure directives in an `ext [ ... ]` list, so the script reads `task.ext.args` and `task.ext.prefix`. The default prefix is built from `id`.
6. No `when:` section and no `ext.when`.
7. File inputs are `Path`, `List<Path>` when order matters and `Bag<Path>` when it does not. Outputs that are `files(...)` globs are `Set<Path>`; the workflow converts them to a sorted `List` when order matters (for example fastp reads).

The module `modules/nf-core/seqkit/stats/main.nf` shows the spec. `modules/local/gfatools/gfa2fa.nf` is a local one.

`id` identifies one task. For a sample-level process it is the sample id. For an assembly-level process it is the assembly build (`hifiasm-raw-default`) or a haplotype key (`<sample>:<build>:hap1`). The output prefix is set explicitly where a file name used to depend on the old `ext.prefix` selector.

Typed Nextflow has no `.out`, `set`, `tap`, `|`, implicit closure parameters, `Record.keySet()`, regex `Matcher`s or `List.minus`. The converted code uses explicit parameters, `findAll`, `collectMany` and so on. A `List<Path>` interpolated in a script is rendered as `[a, b]`; join it, as in `${reads.join(' ')}`.

## Data between stages

Records are defined next to the workflow that uses them. They are structural: a workflow accepts any record with the declared fields.

| Name | Fields |
| --- | --- |
| `Reads` | `id`, `meta`, `reads: List<Path>` |
| `Assembly` | `id` (the build), `meta` (includes `meta.assembly`), `pri_fasta`, `alt_fasta?`, `pri_gfa?`, `alt_gfa?` |
| `FastkDatabase` | `id`, `meta`, `hist`, `ktab: List<Path>`, `prof: List<Path>` |
| `MerylDatabase` | `id`, `meta`, `meryl_db` |
| `KmerCoverage` | `id`, `kmer_cov?` |

`meta` is a record with `id` (the sample id with spaces replaced by `_`), `sample` (`name`, `genome_size`, `ploidy`, `haploid_number`, `tax_id`, ...) and `settings` (`busco.lineages`). Reads add `single_end`, and Hi-C pairs add `pair_id`. Assemblies add `assembly` (`assembler`, `stage`, `id`, `build`) and, while a haplotype is processed, `haplotype`. `meta` behaves like a Map at runtime.

Records are compared and hashed by value, and `a + record(x: 1)` adds a field. Join on a named field: `a.join(b, by: 'id')`. Project the right-hand side to the join key and the new fields first, or fields with the same name collide. Convert a Map with `record(map)` (shallow).

An assembly stage changes the build, because the build is used in file names. `setAssemblyStage` in `main.nf` rewrites `meta.assembly.stage` and `meta.assembly.build`.

## The main workflow

`main.nf` is typed. It declares the params (with defaults) and calls the stage workflows in order:

```
PREPARE_INPUT -> CONVERT_FASTQ_CRAM, BUILD_FASTK_*, BUILD_MERYL_*
              -> INSPECT_DATA -> ASSEMBLE -> DECONTAMINATE -> PURGE_DUPLICATES
              -> DVPOLISH -> SCAFFOLD -> SCAFFOLD_CURATION
              -> (ALIGN_RNASEQ) -> COMPARE_ASSEMBLIES, EVALUATE_ASSEMBLY
              -> REPORT_VERSIONS -> ASSEMBLY_REPORT
```

Stages are gated by `params.steps` as before. A disabled stage gets an empty input, so it runs nothing, and the next stage receives the assemblies of the stage before it. `main.nf` builds the `enable` record from `steps`. The planned replacement by flat `enable_*` params is ADR 0001 (#379).

`PREPARE_INPUT` reads the YAML with SnakeYAML, completes the sample with the ENA taxonomy (in the workflow body) and with GOAT (the `GOAT_TAXONSEARCH` process), and emits typed channels.

The k-mer databases, the CRAM conversion and the report still run for every run, as before (#329).

## Pipelines that run on their own

`pipelines/assembly_qc` evaluates assemblies on its own. It reads typed params from a `-params-file` (sample, assembly list, HiFi reads for the k-mer databases, BUSCO and k-mer options), builds on the subworkflows the main workflow uses, and publishes with an `output {}` block. See `pipelines/assembly_qc/README.md`.

```bash
nextflow run pipelines/assembly_qc -profile singularity -params-file qc_parameters.yml --outdir results
```

To include a pipeline in another one, both enable static typing, and the caller passes the params as one record:

```groovy
nextflow.enable.types = true

include { workflow as ASSEMBLY_QC } from './pipelines/assembly_qc/main.nf'

workflow {
    main:
    qc = ASSEMBLY_QC(record(sample: ..., assembly: ..., hifi: ..., busco: ...))
    // qc.merquryfk, qc.busco, qc.gfastats, qc.quast
}
```

The caller sets `validate_params: false` if it uses nf-schema, and loads the container and resource configuration itself: the `nextflow.config` of an included pipeline is not loaded. Including a pipeline in another has not been tested yet; `assembly_qc` has only been run on its own. Other stages are not separate pipelines yet (#424).

## Publishing

`publishDir` is gone. `main.nf` assigns each result channel in `publish:` and places it in the `output {}` block with `path`. The paths reproduce the old tree: static paths such as `${params.stage_dirs.inspect}/fastk`, closures for per-build folders (`03_assembly/<build>`), and `file >> ("dir/")` where only some files were published (`*.stat` of pairtools, `*.mcool`, `*.beddb`, `*.bw`, `*_wTracks.pretext`, BUSCO summaries). The result folder for the evaluation depends on `meta.assembly.stage`.

`outdir` and `publish_mode` are config params: `outputDir = params.outdir` and `workflow.output.mode = params.publish_mode` are set in `nextflow.config`. `storeDir` for the FCS-GX database is unchanged (#365).

## Configuration

`nextflow.config` keeps the config params (`outdir`, `publish_mode`, `cache`, `trace_timestamp`), the container registries and the profiles. `configs/modules.config` only has process settings that are not inputs: `scratch`, `executor` and `storeDir`. `configs/tool_resources.config` is unchanged and its closures read `meta.sample.genome_size`.

A user selector that sets `ext.args` replaces the argument string the workflow passes. Do not add a selector that sets `ext.prefix` or `ext.args` to a config the workflow loads without a reason: it overrides the input field for every user and hides the value the workflow passes.

## Patched nf-core modules

nf-core modules have no typed version upstream, so the converted modules are tracked with `nf-core modules patch`. Each has a `<tool>-<name>.diff` next to `main.nf` and a `patch` entry in `modules.json`. 47 modules are patched.

- Edit `main.nf`, delete the existing `.diff`, then run `pixi run nf-core modules patch <name>`.
- `nf-core modules update <name>` applies the patch on top of the new upstream version and reports a conflict if they overlap.
- `nf-core modules patch` ends with a traceback about `conf/containers_*.config` in a pipeline without `conf/`. The `.diff` and `modules.json` are written correctly before that.

## Behaviour changes

- Requires Nextflow >= 26.10.0. Composition is experimental.
- `ext.when` is not supported and `when:` is removed from every process.
- Process names have a stage alias prefix. Selectors with old stage names match nothing, without an error.
- `SMUDGEPLOT` is joined to `GENESCOPEFK` by sample id instead of by channel order.
- hifiasm ultra-long reads are passed comma-separated (`--ul a,b`); the old code printed the list as `[a, b]`.
- `OATK_SELECTHMM` output is joined to its input by id.
- The FCS-GX database reaches each screen as a dataflow value.
- The merqury, merquryfk, hifiasm, genescopefk and samtools import stubs produce what the module declares, and 25 processes have a `stub:` block. `-stub-run` now runs the whole workflow (106 tasks), and it needs a container engine because `eval()` outputs run on the host in stub mode.
- `SAMTOOLS_MERGE` inputs are sorted by name.
- The HiFi coverage tracks are not published (see above).
- A missing optional output is `null` in a typed workflow, where legacy processes emitted nothing. A process whose optional output is missing fails at the next step that needs it instead of being dropped by `join`.

## Tests and verification

- **Lint:** `pixi run lint` (`nextflow lint`) has no errors and one warning (`projectDir` in `merquryfk`, #422).
- **Stub run:** `pixi run nextflow run main.nf -profile docker,test -stub-run` completes 106 tasks. `pipelines/assembly_qc` has its own stub run on the tiny data (13 tasks), documented in its README.
- **Real data:** a wasp data set (*Ancistrocerus nigricornis*, HiFi and Hi-C) ran on nac under SLURM and Singularity through `inspect` with the typed modules. A later run adds `assemble`.
- **Default nf-test:** `tests/main.nf.test` and its snapshot have not been re-checked after this change. It runs about an hour per test on a laptop and needs an amd64 host.
- **Module tests:** the tests under `modules/*/tests` are still the upstream tuple form. `nf-test.config` only runs `tests/`, so they are not run. Converting them needs regenerated snapshots.
- nf-test needs `nft-utils@0.0.5`. If the download fails with HTTP 406, download the jar with `curl -L` into `.nf-test/plugins/nft-utils/0.0.5/`.

## Known limitations and follow-ups

- The stage gating is still `params.steps` (ADR 0001, #379). The k-mer databases, the CRAM conversion and the report run for every run (#329).
- Only `assembly_qc` is a pipeline of its own; the other stages are subworkflows (#424). The root workflow does not include `assembly_qc` yet.
- `PREPARE_INPUT` reads the YAML with SnakeYAML and parses the GOAT table in the workflow. `Channel<Record>` params only accept a list of records, not the nested YAML.
- The published outputs are still one channel per tool, to keep the tree. Joining the outputs of a stage into one record per sample or assembly is a follow-up.
- `merquryfk` uses `$projectDir/bin/false_duplications.awk`, and the Hi-C curation modules call scripts from the root `bin/` (#422).
- `EXTRACT_CHROMSIZES`, `TOL_SEARCH`, `REPORT_GENOMETRAITS` and `WRITE_FILE` are `exec:` processes.
- The report notebook reads `params.yml` written by `WRITE_FILE` and the `log_files` staged by name.

## How to add a module or a stage

**A module.** Install it with `nf-core modules install`, convert it to the rules above, and run `nf-core modules patch`. For a local module, copy `modules/local/gfatools/gfa2fa.nf`.

**A stage workflow.**

1. Create `subworkflows/local/<name>/main.nf` with `nextflow.enable.types = true`, record types for its inputs and a `workflow` with typed `take:` and `emit:`.
2. Build the process records in the workflow, with `<tool>_args` and `output_prefix` as fields.
3. Emit the result channels that should be published.
4. Call it in `main.nf` behind its `enable` flag, and add the `publish:` assignment and the `output {}` entry.
5. Check it with `pixi run lint` and `-stub-run`.
