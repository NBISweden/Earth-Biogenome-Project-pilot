# Migration to typed records and stage pipelines

This page describes the change from tuple-based processes and subworkflows to statically typed Nextflow with records, and from one workflow to pipelines that can run on their own. It is written for people who run the workflow and for contributors who know the code from before this change. The decisions are recorded in [ADR 0002](decisions/0002-typed-records-and-stage-pipelines.md) and [ADR 0003](decisions/0003-wide-records-and-stage-pipelines.md).

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
| Nextflow >= 25.10.2 | Nextflow >= 26.10.3 |
| Tuples such as `[ meta, files ]`; `meta` is a Map | Records such as `record(id: ..., meta: ..., reads: ...)`; `meta` is a Record |
| `when:` and `ext.when` on 75 processes | No `when:`. The workflow decides |
| `ext.args` and `ext.prefix` selectors in `configs/modules.config` | `<tool>_args` and `output_prefix` fields set by the workflows |
| `publishDir` selectors in `configs/modules.config` | A `publish:` section and an `output {}` block in `main.nf` |
| Params defaults in `nextflow.config`, read anywhere through `params` | Params declared and typed in `main.nf`; subworkflows get values as inputs |
| `.out` and legacy operators (`branch`, `multiMap`, `groupTuple`, `combine(by:)`) | Variables, `filter`/`map`/`join(by:)`/`groupBy` |
| Helper functions `combineByMetaKeys`, `joinByMetaKeys`, `deepMergeMaps`, ... | Removed. Joins use named record fields |
| `REPORT_DTOL` process | `TOL_SEARCH` also writes `DToL.tsv`. `ENA_TAXQUERY` stays a process (typed, `exec:`), so the lookup is cached with `-resume` |
| `collectFile` for the versions file and report params | `WRITE_FILE` process |
| One workflow; stages are subworkflows with several inputs | The workflow, and one pipeline per stage in `pipelines/`, each taking one wide record per sample or assembly |
| Assembly QC evaluated all stages at the end, inside `main.nf` | `pipelines/assembly_qc`, included by each stage pipeline under an alias, so QC runs per stage |
| Result folders by tool (`03_assembly/busco/<build>`) | Result folders by subject (`03_assembly/<build>/qc/busco`) |

## If you run the workflow

- **Install Nextflow >= 26.10.3.** The manifest stops older versions.
- **The input and the params are the same.** `--input` is still the nested YAML (`sample`, `hifi`, `hic`, `assembly`, ...). `-params-file` keys keep their names, including nested ones such as `fastk.kmer_size` and `busco.lineages`.
- **Nested params are records without field defaults.** If you give `fastk: { kmer_size: 21 }` you override the whole `fastk` record. Fields you leave out are `null`, and the workflow applies its default (31 for the k-mer sizes). The `busco` and `fcs` records work the same way.
- **The published tree changed.** Every stage folder has one subfolder per sample or assembly build, and the tools are inside it. See [Publishing](#publishing) for the moved paths. The HiFi coverage tracks were never published, because they matched a pattern (`*.hitile`) that no module produces.
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

Record types are shared. They live in `types/records.nf` and scripts import them (`include { Assembly } from '../../types/records'`). A type declared in each script makes `Channel<Assembly>` from two scripts incompatible.

| Name | Fields |
| --- | --- |
| `Reads` | `id`, `meta`, `reads: List<Path>`. One file or pair, used inside the read subworkflows |
| `Sample` | `id`, `meta`, `hifi: List<Path>`, `hifi_merged`, `hic: List<List<Path>>`, `rnaseq`, `isoseq`, the HiFi FastK and Meryl databases, the Hi-C FastK database. One per sample |
| `Assembly` | `id` (the build), `meta` (includes `meta.assembly`), `pri_fasta`, `alt_fasta?`, `pri_gfa?` (the contig graph of the raw assembly), plus the sample data the stages need: `hifi`, `hifi_merged`, `hic`, `rnaseq`, the HiFi FastK and Meryl databases. One per build |
| `FastkDatabase`, `MerylDatabase` | `id`, `meta` and the database files. Used by the database subworkflows |

`Sample` and `Assembly` are **wide**: a stage takes one record per subject with everything it needs and does not join its inputs. The joins happen once in `main.nf`, where the sample data is added to the assemblies. Fields are flat (paths, lists of paths, strings), so a stage can run on its own from a YAML samplesheet.

`meta` is a record with `id` (the sample id with spaces replaced by `_`), `sample` (`name`, `genome_size`, `ploidy`, `haploid_number`, `tax_id`, ...) and `settings` (`busco.lineages`). Reads add `single_end`, and Hi-C pairs add `pair_id`. Assemblies add `assembly` (`assembler`, `stage`, `id`, `build`). `meta` behaves like a Map at runtime.

Rules for code that handles these records:

- Join on a named field: `a.join(b, by: 'id')`. Project the right-hand side to the join key and the new fields first, or fields with the same name collide.
- `join(by:)` on records is many-to-many: every record with the key on one side pairs with every record with the key on the other. (It pairs once per key for tuples, which is not what typed code uses.) `combine` has no `by:` for records.
- A function with a declared record return type (`-> Assembly`) returns an instance that `join`, `combine` and `+` reject. The helpers in `types/records.nf` (`withFastas`) and `main.nf` do not declare a return type. They take an `Assembly` and return `a + record(...)`, which the type checker types as the full `Assembly`.
- A function that takes an `Assembly` accepts a record with more fields, and `a + record(...)` replaces fields and keeps the others. `record + record` is how a stage changes an assembly. A bare `null` in a record has no type for the checker: write `null as Path`.
- A stage returns its assemblies as the wide records of its input, changed with `+` (`withFastas`), and the tool results in separate output channels. Tool results are not added to the wide record.
- `Channel<Assembly>` params of a pipeline run on their own accept a samplesheet file (YAML, JSON or CSV) with flat rows. Nested named records do not load.
- An assembly stage changes the build, because the build is used in file names. `setAssemblyStage` in `main.nf` rewrites `meta.assembly.stage` and `meta.assembly.build`.

## The main workflow

`main.nf` is typed. It declares the params (with defaults), prepares one `Sample` per sample, and chains the stage pipelines:

```
PREPARE_INPUT -> CONVERT_FASTQ_CRAM, BUILD_FASTK_*, BUILD_MERYL_* -> INSPECT_DATA
              -> ASSEMBLE -> DECONTAMINATE -> PURGE_DUPLICATES -> DVPOLISH -> SCAFFOLD
              -> SCAFFOLD_CURATION
              -> ALIGN_RNASEQ
              -> REPORT_VERSIONS -> ASSEMBLY_REPORT
```

Each stage is a pipeline in `pipelines/` (`assemble`, `decontaminate`, `purge_duplicates`, `polish`, `scaffold`, `curate`, `align_rnaseq`) included as `workflow as <NAME>`. Every stage that produces assemblies includes `pipelines/assembly_qc` under an alias (`ASSEMBLY_QC_RAW`, `ASSEMBLY_QC_DECONTAMINATED`, `ASSEMBLY_QC_PURGED`, `ASSEMBLY_QC_POLISHED`, `ASSEMBLY_QC_SCAFFOLDED`) and returns the QC results next to the assemblies. `main.nf` does not call the evaluation itself.

Stages are gated by `params.steps` as before. A disabled stage gets an empty input, so it runs nothing, and the next stage receives the assemblies of the stage before it. Assemblies supplied in the input YAML enter the stage that matches their `stage`. They are evaluated by that stage's QC (`ASSEMBLY_QC_CURATED` in `curate`, which makes tracks and no assemblies) and passed on unchanged. `main.nf` builds the `enable` record from `steps`. The planned replacement by flat `enable_*` params is ADR 0001 (#379).

`PREPARE_INPUT` reads the YAML with SnakeYAML, completes the sample with the ENA taxonomy (the `ENA_TAXQUERY` process, only when a taxonomy key is missing) and with GOAT (the `GOAT_TAXONSEARCH` process), and emits one record per sample.

The k-mer databases, the CRAM conversion and the report still run for every run, as before (#329).

## Pipelines that run on their own

Every stage in `pipelines/` is an ordinary pipeline. It declares typed `params` (its input channel, an optional `provided` channel for assemblies supplied at that stage, the tool options and one `qc` record with the options of the assembly QC: `busco`, `fastk`, `meryl` and `use_phased`) and an `output {}` block. A caller passes the params as one record, and the outputs come back as channels:

```groovy
nextflow.enable.types = true

include { workflow as PURGE_DUPLICATES } from './pipelines/purge_duplicates'

workflow {
    main:
    purge = PURGE_DUPLICATES(record(assemblies: assemblies, qc: qc))
    // purge.assemblies, purge.purge_dups, purge.busco, purge.logs
}
```

Only the calling pipeline publishes. A stage run on its own publishes with its own `output {}` block, to `<tool>/<build>`.

A stage runs on its own from a samplesheet of flat `Assembly` rows. `pipelines/<stage>/tests/` has an example for the tiny data set:

```bash
nextflow run pipelines/purge_duplicates -profile docker,test -stub-run \
    -params-file pipelines/purge_duplicates/tests/params.yml -c configs/test_resources.config
```

`pipelines/assembly_qc` also takes a friendlier input for a run on its own: `sample`, an `assembly` list and `hifi` reads in a `-params-file`. See its README. It builds the HiFi k-mer databases only for assemblies that arrive without them, because the task hash includes the include alias and a database built inside each alias would be rebuilt for every stage.

A pipeline that includes these has to enable static typing. The `nextflow.config` of an included pipeline is not loaded, so the caller loads the container and resource configuration itself.

## Publishing

`publishDir` is gone. `main.nf` assigns each result channel in `publish:` and places it in the `output {}` block with `path`. Each stage folder has one subfolder per subject: the sample for the read stages, the assembly build for the assembly stages. The tools are inside it, and the QC of an assembly is in `qc/`:

```
01_read_inspection/<sample>/{goat,seqkit_hifi_stats,seqkit_hic_stats,fastqc_hic,fastk,genescopefk,smudgeplot,katgc,kat_comp,dtol_search}
02_read_preprocessing/<sample>/{hic_cram,converted_bam_to_fasta}
03_assembly/<build>/{hifiasm,gfa2fa,organelle/*,qc/{merquryfk,merqury,busco,gfastats}}
04_contamination_screen/<build>/{fcsgx,qc/...}
05_duplicate_purging/<build>/{purge_dups,qc/...}
06_error_polishing/<build>/{dvpolish,qc/...}
07_scaffolding/<build>/{pairtools,yahs,merged_haps,qc/...}
08_curation/<build>/{higlass,pretext}
09_RNA_alignment/<build>/star
10_report/quast/<sample>            QUAST of the assemblies of all stages together
10_report/{assembly_report.html,assembly_report.md,multiqc-summary.html,versions.yml}
```

Moved paths, with `<stage>` the folder of the stage that made the assembly:

| Before | After |
| --- | --- |
| `01_read_inspection/<tool>/` | `01_read_inspection/<sample>/<tool>/` |
| `02_read_preprocessing/hi-c_cram/` | `02_read_preprocessing/<sample>/hic_cram/` |
| `03_assembly/<build>/` (hifiasm files and `gfa2fa`) | `03_assembly/<build>/hifiasm/` and `03_assembly/<build>/gfa2fa/` |
| `03_assembly/organelle/{mitohifi,oatk,dotplots}` | `03_assembly/<build or sample>/organelle/{mitohifi,oatk,dotplots}` |
| `04_contamination_screen/fcsgx/` | `04_contamination_screen/<build>/fcsgx/` |
| `05_duplicate_purging/purge_dups/` | `05_duplicate_purging/<build>/purge_dups/` |
| `06_error_polishing/dvpolish/<build>/` | `06_error_polishing/<build>/dvpolish/` |
| `06_error_polishing/merqury/<build>/<sample>_consensus*`, `…_<hap>.purged*` (the Merqury runs inside DVPOLISH) | `06_error_polishing/<build>/dvpolish/{merqury_input,merqury_polished}/` |
| `07_scaffolding/{pairtools,yahs,merged_haps}/<build>/` | `07_scaffolding/<build>/{pairtools,yahs,merged_haps}/` |
| `08_curation/{higlass,pretext}/<build>/` | `08_curation/<build>/{higlass,pretext}/` |
| `09_RNA_alignment/star/<build>/` | `09_RNA_alignment/<build>/star/` |
| `<stage>/{merquryfk,merqury,busco,gfastats}/<build>/` | `<stage>/<build>/qc/{merquryfk,merqury,busco,gfastats}/` |
| `03_assembly/mitohifi/references` | `03_assembly/<sample>/organelle/mitohifi/references` |
| `10_report/quast/` (all stages together) | `10_report/quast/<sample>/` (unchanged: one run per sample over all stages, from `pipelines/compare_stages`) |

The paths are built from closures such as `${subjectDir(params.stage_dirs.purge, r.meta)}/purge_dups`, and `file >> ("dir/")` where only some files are published (`*.stat` of pairtools, `*.mcool`, `*.beddb`, `*.bw`, `*_wTracks.pretext`, BUSCO summaries).

`outdir` and `publish_mode` are config params: `outputDir = params.outdir` and `workflow.output.mode = params.publish_mode` are set in `nextflow.config`. `storeDir` for the FCS-GX database is unchanged (#365).

## Configuration

`nextflow.config` keeps the config params (`outdir`, `publish_mode`, `cache`, `trace_timestamp`), the container registries and the profiles. `configs/modules.config` only has process settings that are not inputs: `scratch`, `executor` and `storeDir`. `configs/tool_resources.config` keeps its selectors. Its closures read `meta.sample.genome_size`, and the memory and time of eight measured tools now scale with it (see `docs/development.md`, "Resources"). Killed tasks are retried twice.

A user selector that sets `ext.args` replaces the argument string the workflow passes. Do not add a selector that sets `ext.prefix` or `ext.args` to a config the workflow loads without a reason: it overrides the input field for every user and hides the value the workflow passes.

## Patched nf-core modules

nf-core modules have no typed version upstream, so the converted modules are tracked with `nf-core modules patch`. Each has a `<tool>-<name>.diff` next to `main.nf` and a `patch` entry in `modules.json`. 47 modules are patched.

- Edit `main.nf`, delete the existing `.diff`, then run `pixi run nf-core modules patch <name>`.
- `nf-core modules update <name>` applies the patch on top of the new upstream version and reports a conflict if they overlap.
- `nf-core modules patch` ends with a traceback about `conf/containers_*.config` in a pipeline without `conf/`. The `.diff` and `modules.json` are written correctly before that.

## Behaviour changes

- Requires Nextflow >= 26.10.3. Composition is experimental.
- `ext.when` is not supported and `when:` is removed from every process.
- Process names have a stage alias prefix, and the QC processes sit under the stage that includes them (`PURGE_DUPLICATES:ASSEMBLY_QC_PURGED:EVALUATE_ASSEMBLY:BUSCO`). Selectors with old stage names match nothing, without an error.
- The published tree changed (see [Publishing](#publishing)).
- HiFi reads are required. The input YAML needs `hifi`, and `assembly_qc` fails for an assembly that has neither HiFi k-mer databases nor HiFi reads. Before, an assembly without reads was evaluated without the k-mer tools.
- **QC is no longer held back.** The old workflow delayed all assembly QC (MerquryFK, Merqury, BUSCO, GFAStats) until every assembly stage had finished, on purpose, so that the processing steps had priority on the cluster. This conversion does not do that: each stage pipeline evaluates its assemblies as soon as they exist, so QC can compete with processing. Only the QUAST comparison across the stages waits for the last stage. The overlap is accepted (ADR 0003).
- `ASSEMBLE` waits for the HiFi k-mer databases and the Hi-C conversion, because the assemblies carry them and the QC reads them. The databases used to be built alongside the assembly.
- `PURGE_DUPLICATES` does not wait for the k-mer coverage estimate of `INSPECT_DATA` any more. The value was never used. 
- hifiasm (`hifiasm_with_hic`) reads the converted Hi-C pairs, like the later stages. With an Arima `hic_type` these are trimmed (`--trim_front1/2=5`). This is intended. `INSPECT_DATA` reads the raw pairs.
- Hi-C reads are optional. A run goes as far as it can without them. If `--steps` includes `scaffold`, `curate` or `assemble` with `hifiasm_with_hic`, the run fails at the start with a message. The standalone `scaffold` and `curate` pipelines fail for an assembly without `hic`.
- FCS-GX runs from a ramdisk by default (`/dev/shm/gxdb`). Set `fcs.ramdisk_path: ''` to turn it off.
- Resource selectors in `configs/tool_resources.config` are plain process names, so they also apply when a stage pipeline runs on its own. A selector on the original process name matches every alias of it (tested), so `BWAMEM2_MEM` still covers `BWAMEM2_MEM_SCAFFOLD` and `BWAMEM2_MEM_CURATION`.
- `--input` takes one YAML file. It used to accept a glob that matched several files, which was not documented.
- `chunk_size` (polish) is a memory unit param. Give it as `1.5 GB` or `100.MB` without quotes. `MemoryUnit.of` is not used.
- `SMUDGEPLOT` is joined to `GENESCOPEFK` by sample id instead of by channel order.
- hifiasm ultra-long reads are passed comma-separated (`--ul a,b`); the old code printed the list as `[a, b]`.
- `OATK_SELECTHMM` output is joined to its input by id.
- The FCS-GX database reaches each screen as a dataflow value.
- The merqury, merquryfk, hifiasm, genescopefk and samtools import stubs produce what the module declares, and 24 of the 29 local processes have a `stub:` block (23 more than on `main`; the 47 nf-core processes all have one). `-stub-run` runs the whole workflow, and it needs a container engine, because the `eval()` outputs of the version topics run in the container.
- The stage pipelines take the QC options as one `qc` param (`qc.busco`, `qc.fastk`, `qc.meryl`, `qc.use_phased`) instead of four. `main.nf` keeps the flat `busco`, `fastk`, `meryl` and `use_phased` params. A params file for a stage on its own, as in `pipelines/<stage>/tests/params.yml`, puts them under `qc:`.
- The `alt_gfa` field and the `alt_gfa` key of the input assemblies are gone: nothing read them. `pri_gfa` is kept for the raw assembly, where the organelle step reads it.
- The stage pipelines share one test resource config, `configs/test_resources.config` (run a stage with `-c configs/test_resources.config`).
- `OATK` runs once per hifiasm build. The old workflow zipped its inputs by position, so with several builds and no plastid HMM it ran once.
- STAR writes one set of files per RNA-seq library (`<sample>_<library>`); a run with several libraries would have written the same file names.
- A `stage_dirs` param that names only some stages keeps the default folders of the others.
- The inputs of `SAMTOOLS_MERGE` are sorted by file name. The lists of FastK and Meryl tables, pairtools and curation BAMs, and polish VCFs are sorted with `toSorted()`, which orders Paths by full path. That is stable between identical runs, but is not file-name order.
- The HiFi coverage tracks are not published (see above).
- A missing optional output is `null` in a typed workflow, where legacy processes emitted nothing. A process whose optional output is missing fails at the next step that needs it instead of being dropped by `join`.

## Tests and verification

- **Lint:** `pixi run lint` (`nextflow lint`) has no errors and one warning (`projectDir` in `merquryfk`, #422). The type checker at launch (`Type checking found N error(s)`) can report errors that lint does not. Check `.nextflow.log` for `Type error`.
- **Stub tests:** `pixi run nftest-stub` runs six stub tests of `main.nf` on the tiny data: every assembly and polishing region reaches every stage, one QUAST per sample, QC of a supplied curated assembly, Hi-C handling and an unknown stage. They check wiring, not tool output. Each stage pipeline also has a stub run on the tiny data in `pipelines/<stage>/tests/`.
- **Local module tests:** `pixi run nftest-modules` runs 13 nf-tests of 10 local modules in `tests/modules/local/`. Eight modules are tested with `-stub-run` and assert on output file names (a snapshot would hold the absolute work directory). `extract/chromsizes` is an `exec:` process and runs for real. `ena/taxquery` has no stub: its two tests query the ENA service (tag `network`). `pixi run nftest-functions` tests the offline input helpers.
- **Real data:** a wasp data set (*Ancistrocerus nigricornis*, HiFi, Hi-C and RNA-seq) runs on nac under SLURM and Singularity through every stage except `screen` (FCS-GX needs more memory than the nac nodes have). The first round ran `inspect`, `preprocess` and `assemble`. The result of the full run is in the PR description.
- **Default nf-test:** `tests/main.nf.test` passes with the `apptainer` profile on an amd64 cluster node (pelle, 12 CPUs, about 25 minutes for the four tests) and its snapshots are regenerated for the new tree. The tests run 107, 76, 16 and 55 tasks. Compared with `main`, the published file names differ only by the paths, the number of hidden FastK table parts (one per thread, and FastK now runs with 8 CPUs), more `unionsum.hist.ploidy` and polish-stage Merqury files (the Merqury runs inside DVPOLISH are published again in `dvpolish/merqury_input` and `merqury_polished`), and the sample folder in `10_report/quast/<sample>`. QUAST is one run per sample over all stages, as on `main`. It runs about an hour per test on a laptop. Real tools exposed bugs that stubs hide: `merqury` and `merquryfk` declared one file for outputs they write once per assembly.
- **Vendored nf-core tests:** the upstream tests of the vendored nf-core modules (`modules/nf-core/**/tests/`, and `subworkflows/nf-core/**/tests/` if any are installed) are no longer tracked and are gitignored. They are tested upstream, are in the tuple form, and `nf-test.config` never ran them. `nf-core modules install` and `update` still write them locally, and `nf-core modules patch` records them as "removed" in the `.diff` trailer, which has no effect on patching.
- nf-test needs `nft-utils@0.0.5`. If the download fails with HTTP 406, download the jar with `curl -L` into `.nf-test/plugins/nft-utils/0.0.5/`.

## Known limitations and follow-ups

- Launching `main.nf` prints `Type checking found 3 error(s)`. They come from the `chunk_size: MemoryUnit` param (the type checker of Nextflow 26.10.3 and 26.10.4 treats the declared `MemoryUnit` and the `100.MB` literal as different types and does not know `toBytes`). `nextflow lint` reports none and the values are correct at runtime (`1.5 GB` is 1,610,612,736 bytes).
- The stage gating is still `params.steps` (ADR 0001, #379). The k-mer databases, the CRAM conversion and the report run for every run (#329).
- `PREPARE_INPUT` reads the YAML with SnakeYAML and parses the GOAT table in the workflow. `Channel<Record>` params only accept a samplesheet file with flat rows, not the nested YAML.
- The sample-level subworkflows (`CONVERT_FASTQ_CRAM`, the database builds, `INSPECT_DATA`) are not pipelines of their own, and the report is not either.
- A stage run on its own needs the sample data (reads, databases) in each assembly row, as `pipelines/<stage>/tests/` shows.
- The stage tools publish through one output channel per tool. Joining them into one record per assembly is a follow-up.
- `merquryfk` uses `$projectDir/bin/false_duplications.awk`, and the Hi-C curation modules call scripts from the root `bin/` (#422).
- `EXTRACT_CHROMSIZES`, `TOL_SEARCH`, `REPORT_GENOMETRAITS` and `WRITE_FILE` are `exec:` processes.
- The report notebook reads `params.yml` written by `WRITE_FILE` and the `log_files` staged by name.
- `configs/tool_resources.config` reads `meta.sample.genome_size` in resource closures, so a process that uses those labels needs a `meta` with `sample.genome_size` even though its code never reads `meta`.

## How to add a module or a stage

**A module.** Install it with `nf-core modules install`, convert it to the rules above, and run `nf-core modules patch`. For a local module, copy `modules/local/gfatools/gfa2fa.nf`.

**A stage pipeline.**

1. Copy `pipelines/purge_duplicates`: `main.nf` with `nextflow.enable.types = true`, a `params` block with `assemblies: Channel<Assembly>`, `provided: Channel<Assembly>?` and the tool and QC options, an entry workflow and an `output {}` block, and `tests/{params.yml,assemblies.yml,test.config}`.
2. Build the process records from the fields of the wide record, with `<tool>_args` and `output_prefix` as fields. Do not join on the sample inside the stage.
3. Return the assemblies as `a + record(...)` (see `withFastas`) and the tool results as separate output channels. Include `pipelines/assembly_qc` under an alias if the stage makes assemblies.
4. In `main.nf`, include it as `workflow as <NAME>`, call it with a record behind its `enable` flag, and add the `publish:` assignment and the `output {}` entry with the path of its subject folder.
5. Check it with `pixi run lint` and a `-stub-run` of the pipeline and of `main.nf`.
