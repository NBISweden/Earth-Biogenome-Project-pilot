nextflow.enable.types = true

include { ENA_TAXQUERY           } from "../../../modules/local/ena/taxquery/main"
include { GOAT_TAXONSEARCH        } from "../../../modules/nf-core/goat/taxonsearch/main"
include { SAMTOOLS_FASTA          } from "../../../modules/local/samtools/fasta/main"
include { CAT_CAT as MERGE_PACBIO } from "../../../modules/nf-core/cat/cat/main"

/* The input file is a YAML file (--input), for example:
```yaml
sample:
  name: Species name
  # Optional. Looked up in ENA and GOAT when missing:
  # tax_id, genetic_code, mito_code, domain, lineage, genome_size, haploid_number, ploidy
assembly:
  - assembler: hifiasm
    stage: raw
    id: uuid
    pri_fasta: /path/to/primary_asm.fasta
    alt_fasta: /path/to/alternate_asm.fasta
    pri_gfa: /path/to/primary_asm.gfa
hic:
  - read1: /path/to/read1.fastq.gz
    read2: /path/to/read2.fastq.gz
hifi:
  - reads: /path/to/reads
rnaseq:
  - reads: /path/to/reads        # or read1 and read2
isoseq:
  - reads: /path/to/reads
```
*/

/*
The meta record of a sample:
```
id: "Species_name",
sample: [ name: "Species name", genome_size: 1234, ploidy: 2, haploid_number: 13, tax_id: ..., ... ],
settings: [ busco: [ lineages: "metazoa_odb10,..." ] ]
```
Meta data of the reads adds single_end (and pair_id for Hi-C). Meta data of the assemblies adds
assembly: [ assembler: "hifiasm", stage: "raw", id: uuid, build: "hifiasm-raw-uuid" ].
*/

include { Reads } from "../../../types/records"

def readYaml(yaml_file: Path) -> Map {
    return new org.yaml.snakeyaml.Yaml().load(yaml_file.text) as Map
}

// The ENA taxonomy is looked up unless the user gave all of its keys
def needsEna(sample: Map) -> Boolean {
    def ena_keys = [ 'tax_id', 'genetic_code', 'mito_code', 'domain', 'lineage' ]
    return !ena_keys.every { key -> sample[key] != null }
}

// The sample with the ENA taxonomy for the keys the user did not give
def sampleWithEna(sample: Map, ena: Map) -> Map {
    def lineage = ena['lineage'] as String
    return sample + [
        tax_id      : sample.tax_id ?: ena['taxId'],
        genetic_code: sample.genetic_code ?: ena['geneticCode'],
        mito_code   : sample.mito_code ?: ena['mitochondrialGeneticCode'],
        domain      : sample.domain ?: lineage.tokenize(';')[0],
        lineage     : sample.lineage ?: lineage
    ]
}

// Look up the genome properties and BUSCO lineages in GOAT for eukaryotes, unless all is user defined
def useGoat(sample: Map, busco_lineages: String?) -> Boolean {
    def goat_keys = [ 'genome_size', 'haploid_number', 'ploidy' ]
    return sample.domain == "Eukaryota" && ( !goat_keys.every { key -> sample[key] != null } || !busco_lineages )
}

def sampleId(sample: Map) -> String {
    return (sample.name as String).replace(" ", "_")
}

// The sample meta record, with the genome properties and BUSCO lineages from GOAT when it was queried
def sampleMeta(sample: Map, goat_tsv: Path?, busco_lineages: String?) -> Record {
    def id = (sample.name as String).replace(" ", "_")
    if ( goat_tsv == null ) {
        return record(id: id, sample: sample, settings: [ busco: [ lineages: busco_lineages ] ])
    }
    def rows = goat_tsv.splitCsv( sep: "\t", header: true, quote: '"' )
    def species = rows.find { row -> row.scientific_name == sample.name } as Map
    assert species != null : "Failed to parse species '${sample.name}' in GOAT_TAXONSEARCH output .tsv"
    def goat_lineages = rows.findAll { row -> row.odb10_lineage }.collect { row -> row.odb10_lineage }.join(',')
    return record(
        id      : id,
        sample  : sample + [
            genome_size   : sample.genome_size ?: species['genome_size'],
            haploid_number: sample.haploid_number ?: species['haploid_number'],
            ploidy        : sample.ploidy ?: species['ploidy']
        ],
        settings: [ busco: [ lineages: busco_lineages ?: goat_lineages ] ]
    )
}

workflow PREPARE_INPUT {

    take:
    input         : Path    // YAML file with the sample and the data
    busco_lineages: String? // Comma separated BUSCO lineages, or null to use GOAT
    steps         : Map<String,Boolean> // The requested steps: assemble, assemble_hic (hifiasm with Hi-C), scaffold, curate

    main:
    def data = readYaml(input)
    def given_sample = data.sample as Map

    // Taxonomy from ENA for the keys the user did not give. A process, so the result is cached
    def needs_ena = needsEna(given_sample)
    ena = ENA_TAXQUERY(
        channel.of( record(id: sampleId(given_sample), species_name: given_sample.name as String) )
            .filter { _x -> needs_ena }
    )
    samples_with_taxonomy = ena
        .map { e -> sampleWithEna(given_sample, e.json as Map) }
        .mix( channel.of(given_sample).filter { _x -> !needs_ena } )

    goat = GOAT_TAXONSEARCH(
        samples_with_taxonomy
            .filter { s -> useGoat(s, busco_lineages) }
            .map { s -> record(id: sampleId(s), taxon: s.name as String, taxa_file: null, goat_args: "--lineage --variables odb10_lineage,haploid_number,genome_size,ploidy") }
    )
    metadata = goat
        .map { g -> record(id: g.id, taxonsearch: g.taxonsearch) }
        .join(
            samples_with_taxonomy
                .filter { s -> useGoat(s, busco_lineages) }
                .map { s -> record(id: sampleId(s), sample: s) },
            by: 'id'
        )
        .map { g -> sampleMeta(g.sample as Map, g.taxonsearch, busco_lineages) }
        .mix(
            samples_with_taxonomy
                .filter { s -> !useGoat(s, busco_lineages) }
                .map { s -> sampleMeta(s, null, busco_lineages) }
        )

    // Assemblies supplied as input
    def assembly_inputs = ( data.assembly ?: [] ) as List
    assemblies = metadata.flatMap { m ->
        assembly_inputs.collect { a ->
            def asm = a as Map
            def build = "${asm.assembler}-${asm.stage}-${asm.id}"
            record(
                id       : build,
                sample_id: m.id,
                meta     : m + record(assembly: [ assembler: asm.assembler, stage: asm.stage, id: asm.id, build: build ]),
                pri_fasta: file(asm.pri_fasta as String, checkIfExists: true),
                alt_fasta: asm.alt_fasta ? file(asm.alt_fasta as String, checkIfExists: true) : null,
                pri_gfa  : asm.pri_gfa ? file(asm.pri_gfa as String, checkIfExists: true) : null
            )
        }.toList()
    }

    // Hi-C read pairs
    def hic_pairs = ( data.hic ?: [] ).collect { pair ->
        [ file((pair as Map).read1 as String, checkIfExists: true), file((pair as Map).read2 as String, checkIfExists: true) ]
    }.toList()

    // Supplied assemblies must name a known stage: the stage decides where they enter the workflow
    def known_stages = [ 'raw', 'decontaminated', 'purged', 'polished', 'scaffolded', 'curated' ]
    assembly_inputs.each { a ->
        def stage = (a as Map).stage
        if ( !(stage in known_stages) ) {
            error "The assembly ${(a as Map).id} in ${input} has the stage '${stage}'. The stage is one of: ${known_stages.join(', ')}."
        }
    }

    // Hi-C reads are optional, but a step that will run and needs them cannot go on without.
    // Scaffolding works on assemblies of the stages before it, and curation on those before it. The
    // assemble step makes assemblies that reach every later stage.
    def supplied_stages = assembly_inputs.collect { a -> (a as Map).stage }
    def scaffold_needs_hic = steps.scaffold && ( steps.assemble || supplied_stages.any { s -> s in [ 'raw', 'decontaminated', 'purged', 'polished' ] } )
    def curate_needs_hic = steps.curate && ( steps.assemble || supplied_stages.any { s -> s in [ 'raw', 'decontaminated', 'purged', 'polished', 'scaffolded' ] } )
    def hic_steps = [
        steps.assemble_hic ? 'assemble (hifiasm_with_hic)' : null,
        scaffold_needs_hic ? 'scaffold' : null,
        curate_needs_hic ? 'curate' : null
    ].findAll { step -> step != null }
    if ( !data.hic && hic_steps ) {
        error "The input ${input} has no Hi-C reads. The steps ${hic_steps.join(', ')} need them. Add the reads under `hic`, or leave those steps out of --steps."
    }

    // PacBio HiFi reads are required: the k-mer databases and the assembly QC are built from them
    if ( !data.hifi ) {
        error "The input ${input} has no HiFi reads. Add them under `hifi`."
    }

    // PacBio HiFi reads. Convert BAMs to FastQ
    def hifi_files = ( data.hifi ?: [] ).collect { h -> file((h as Map).reads as String, checkIfExists: true) }.toList()
    hifi_bam = metadata.flatMap { m ->
        hifi_files
            .findAll { f -> f.name.endsWith(".bam") }
            .collect { bam -> record(id: "${m.id}:${bam.baseName}", meta: m + record(single_end: true), single_end: true, bam: bam, output_prefix: bam.baseName) }
            .toList()
    }
    converted_bams = SAMTOOLS_FASTA( hifi_bam )
    hifi_converted = converted_bams
        .flatMap { r -> r.fasta.collect { fasta -> record(id: r.meta.id, meta: r.meta, reads: [ fasta ]) }.toList() }
    hifi = metadata
        .flatMap { m ->
            hifi_files
                .findAll { f -> !f.name.endsWith(".bam") } // assume everything else is fastx
                .collect { f -> record(id: m.id, meta: m + record(single_end: true), reads: [ f ]) }
                .toList()
        }
        .mix(hifi_converted)

    // Merge the HiFi reads of a sample
    hifi_by_sample = hifi
        .map { r -> tuple(r.id, r) }
        .groupBy()
        .map { sample_id, results ->
            def sorted = results.toList().toSorted { a, b -> a.reads[0].name <=> b.reads[0].name }
            record(id: sample_id, meta: sorted[0].meta, reads: sorted.collect { r -> r.reads[0] }.toList())
        }
    merged = MERGE_PACBIO(
        hifi_by_sample
            .filter { s -> s.reads.size() > 1 }
            .map { s -> record(id: s.id, meta: s.meta, files_in: s.reads, output_prefix: "${s.id}_mergedpb.fasta.gz") }
    )
    hifi_merged = hifi_by_sample
        .filter { s -> s.reads.size() == 1 }
        .map { s -> record(id: s.id, meta: s.meta, reads: s.reads) }
        .mix(merged.map { m -> record(id: m.id, meta: m.meta, reads: [ m.file_out ]) })

    // RNAseq reads: a pair, or a single file
    def rnaseq_pairs = ( data.rnaseq ?: [] ).collect { r ->
        def rec = r as Map
        rec.reads
            ? [ file(rec.reads as String, checkIfExists: true) ]
            : [ file(rec.read1 as String, checkIfExists: true), file(rec.read2 as String, checkIfExists: true) ]
    }.toList()

    // Isoseq reads
    def isoseq_files = ( data.isoseq ?: [] ).collect { r -> file((r as Map).reads as String, checkIfExists: true) }.toList()

    // One record per sample with its reads
    samples = metadata
        .map { m -> record(id: m.id, meta: m, hic: hic_pairs, rnaseq: rnaseq_pairs, isoseq: isoseq_files) }
        .join(hifi_by_sample.map { h -> record(id: h.id, hifi: h.reads) }, by: 'id', remainder: true)
        .join(hifi_merged.map { h -> record(id: h.id, hifi_merged: h.reads[0]) }, by: 'id', remainder: true)

    emit:
    goat       : Channel<Record> = goat
    converted  : Channel<Record> = converted_bams
    metadata   : Channel<Record> = metadata
    samples    = samples
    assemblies = assemblies
}
