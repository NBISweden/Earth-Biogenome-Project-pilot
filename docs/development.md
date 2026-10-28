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
    alt_gfa        : Path?
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

Join channels on a named field (`a.join(b, by: 'id')`). The record types are
defined next to the workflow that uses them.
