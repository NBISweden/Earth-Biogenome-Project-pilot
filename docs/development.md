# Developer information

Improvements in the workflow here are intended to be contributed back to
[nf-core/genomeassembler](https://github.com/nf-core/genomeassembler).
As such, it follows the nf-core way of using modules, and parameter passing.

This means, where possible, tools should be included as modules in the
[nf-core/modules](https://github.com/nf-core/modules) repository.

Most channels will pass meta data through channels along with associated files
so an input `tuple` will look like `[ meta, file(s) ... ]`, where `meta` is
a `Map` ( or dictionary if you prefer ). The meta data map used in this workflow
has the following nested structure:

```
id = {sample.name/ /_}            # Required by nf-core
single_end = {*_reads.single_end} # Required by nf-core
sample = [
    name: <species name>,
    genome_size: <GOAT>,
    ploidy: <GOAT>,
    haploid_number: <GOAT>
]
hifi_reads = [
    single_end: true,
    kmer_cov: <GENESCOPEFK>
]
hic_reads = [
    single_end: false,
    kmer_cov: <GENESCOPEFK>
]
rnaseq_reads = [
    single_end: false
]
isoseq_reads = [
    single_end: false
]
assembly = [
    assembler: {hifiasm},
    stage: {raw,decontaminated,deduplicated,polished,scaffolded,curated},
    id: UUID
    build: assembly.assembler-assembly.stage-assembly.id
]
```

Channels containing assemblies also use a special `Map` to pass around files with
meta data.

```
assembly = [
    assembler: {hifiasm},
    stage: {raw,decontaminated,deduplicated,polished,scaffolded,curated},
    id: UUID
    build: assembly.assembler-assembly.stage-assembly.id
    pri_fasta: <ASSEMBLER>,
    alt_fasta: null/<ASSEMBLER>,
    pri_gfa: <ASSEMBLER>,
    alt_gfa: null/<ASSEMBLER>
]
```