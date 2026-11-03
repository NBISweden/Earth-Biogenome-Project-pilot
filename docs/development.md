# Developer information

Improvements in the workflow here are intended to be contributed back to
[nf-core/genomeassembler](https://github.com/nf-core/genomeassembler).
As such, it follows the nf-core way of using modules, and parameter passing.

This means, where possible, tools should be included as modules in the
[nf-core/modules](https://github.com/nf-core/modules) repository.

Channels carry typed records. A record has named fields, so the channel
content is described by its fields and not by its position in a tuple. See
[the migration guide](migration-typed-records.md) for the rules that every
process and workflow follows, and for how to add a module or a stage.

Most records have an `id`, a `meta` record and the files of a step. The record
types shared by the workflow and the pipelines are in `types/records.nf`. A
stage takes one wide record per subject (`Sample` for the read stages,
`Assembly` for the assembly stages) and does not join its inputs:

```
record Assembly {
    id             : String    // meta.assembly.build
    meta           : Record    // sample meta data, including meta.assembly
    pri_fasta      : Path
    alt_fasta      : Path?
    pri_gfa        : Path?
    hifi           : List<Path>?       // sample data the stages need
    hifi_merged    : Path?
    hic            : List<List<Path>>?
    rnaseq         : List<List<Path>>?
    fastk_hifi_hist: Path?  // and fastk_hifi_ktab, fastk_hifi_prof, meryl_hifi_db
}
```

The `meta` record of a sample has this structure. Nothing in a process reads
it; it is passed through, so workflows and configs can use it.

```
id = {sample.name/ /_}            # Required by nf-core
sample = [
    name: <species name>,
    genome_size: <GOAT>,
    ploidy: <GOAT>,
    haploid_number: <GOAT>
]
settings = [
    busco: [ lineages: <GOAT or params.busco.lineages> ]
]
single_end = {*_reads.single_end} # Added to the meta data of reads
pair_id = {index}                 # Added to the meta data of Hi-C pairs
assembly = [                      # Added to the meta data of assemblies
    assembler: {hifiasm},
    stage: {raw,decontaminated,purged,polished,scaffolded,curated},
    id: {user defined, e.g. opts01}
    build: assembly.assembler-assembly.stage-assembly.id
]
haplotype = {hap1,hap2}           # Added while one haplotype is processed
```

Join channels on a named field (`a.join(b, by: 'id')`). `join` on records is
many-to-many, so a key with several records on one side pairs with all of them.
The shared record types are in `types/records.nf`.

## Tests

Run everything through `pixi run <task>`. Docker or Singularity/Apptainer is
needed, including for stub runs, because `eval()` outputs run in the container.

| Task | What it runs | Real tools |
|---|---|---|
| `lint` | `nextflow lint` on the workflow, pipelines and subworkflows | n/a |
| `nftest-modules` | One test per local module (`tests/modules/local`) | No, stubs. Checks the wiring and output file names. `ena/taxquery` has no stub and queries the ENA service (tag `network`) |
| `nftest-functions` | The offline helpers of `PREPARE_INPUT` (ENA completion of the sample, GOAT decision), `tests/functions` | n/a |
| `nftest-stub` | Fan-out checks of `main.nf` (`tests/main_stub.nf.test`): every assembly and region reaches every stage | No, stubs |
| `nftest-docker`, `nftest-singularity` | `tests/main.nf.test`: four runs on the tiny data with a file-list snapshot | Yes |

The real run takes about 25 minutes in total on a cluster node with 12 CPUs and
about an hour per test on a laptop. Regenerate the snapshot with
`--update-snapshot` and read the diff: a shorter file list is lost output, not
a path change. Run it on a cluster with Singularity or Apptainer (an amd64 host).
The upstream tests of the vendored nf-core modules are not kept (`.gitignore`).

## Resources

`configs/tool_resources.config` sets the resources of the tools. Memory and time of the measured
tools are proportional to the genome size of the sample, at least the measured value with a margin
and at most the value they had before for genomes over 3 Gb. The measurements are from one run on a
0.22 Gb genome (a wasp, HiFi and Hi-C, nac, 16 CPUs and 120 GB per task at most):

| Tool | Peak memory | Run time | Cores used | Note |
|---|---|---|---|---|
| `HIFIASM` | 21.6 GB | 3.8 h | 14 of 16 | |
| `BUSCO` | 5.1 GB | 5.5 min | 7.5 | |
| `MERQURY` | 2.6 GB | 19 min | 1.2 | |
| `MERQURYFK_MERQURYFK` | 4.1 GB | 4.4 min | 2.2 | |
| `GFASTATS` | 0.3 GB | 8 s | 1 | |
| `FASTK_FASTK` | 14 GB (HiFi), 34 GB (Hi-C) | 29 min, 3.5 h | under 1 | the data volume sets the memory |
| `MERYL_COUNT` | 80 GB | 44 min, 4.2 h | 4 | uses the memory it is given |
| `MERYL_UNIONSUM` | 3.1 GB | 19 min | 1 | |

Tools that did not run at this scale (pairtools, YAHS, purge_dups, bwa-mem2, PretextMap, FCS-GX,
DeepVariant, the polish and curation tools) keep their earlier values. One genome is one data
point, so a task killed for memory or time is retried twice with more of both. The trace file
(`pipeline_info/execution_trace_*.txt`) has the requested next to the used resources. Collect
traces of more genomes there, with the genome size, clade and data volume, to fit the allocation
(a model of the resources from the input is future work).
