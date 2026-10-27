nextflow.enable.types = true

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
    alt_gfa: /path/to/alternate_asm.gfa
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

record Reads {
    id   : String
    meta : Record
    reads: List<Path>
}

record Assembly {
    id       : String  // meta.assembly.build
    meta     : Record  // sample meta data, including meta.assembly
    pri_fasta: Path
    alt_fasta: Path?
    pri_gfa  : Path?
    alt_gfa  : Path?
}

def readYaml(yaml_file: Path) -> Map {
    return new org.yaml.snakeyaml.Yaml().load(yaml_file.text) as Map
}

// Species taxonomy from ENA
def fetchEnaTaxonomy(species_name: String) -> Map {
    def encoded_name = java.net.URLEncoder.encode(species_name, 'UTF-8')
    def response = file("https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/${encoded_name}").text
    def entries = new groovy.json.JsonSlurper().parseText(response) as List
    return entries[0] as Map
}

// Complete the sample with the ENA taxonomy for the keys the user did not give
def sampleFromEna(sample: Map) -> Map {
    def ena_keys = [ 'tax_id', 'genetic_code', 'mito_code', 'domain', 'lineage' ]
    if ( ena_keys.every { key -> sample[key] != null } ) {
        return sample
    }
    def ena = fetchEnaTaxonomy(sample.name as String)
    def lineage = ena['lineage'] as String
    return sample + [
        tax_id      : sample.tax_id ?: ena['taxId'],
        genetic_code: sample.genetic_code ?: ena['geneticCode'],
        mito_code   : sample.mito_code ?: ena['mitochondrialGeneticCode'],
        domain      : sample.domain ?: lineage.tokenize(';')[0],
        lineage     : sample.lineage ?: lineage
    ]
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

    main:
    def data = readYaml(input)
    def sample = sampleFromEna(data.sample as Map)

    // Look up the genome properties and BUSCO lineages in GOAT for eukaryotes, unless all is user defined
    def goat_keys = [ 'genome_size', 'haploid_number', 'ploidy' ]
    def use_goat = sample.domain == "Eukaryota" && ( !goat_keys.every { key -> sample[key] != null } || !busco_lineages )
    goat = GOAT_TAXONSEARCH(
        channel.of( record(id: (sample.name as String).replace(" ", "_"), taxon: sample.name as String, taxa_file: null, goat_args: "--lineage --variables odb10_lineage,haploid_number,genome_size,ploidy") )
            .filter { _x -> use_goat }
    )
    metadata = goat
        .map { g -> sampleMeta(sample, g.taxonsearch, busco_lineages) }
        .mix(
            channel.of(sample.name)
                .filter { _x -> !use_goat }
                .map { _x -> sampleMeta(sample, null, busco_lineages) }
        )

    // Assemblies supplied as input
    def assembly_inputs = ( data.assembly ?: [] ) as List
    assemblies = metadata.flatMap { m ->
        assembly_inputs.collect { a ->
            def asm = a as Map
            def build = "${asm.assembler}-${asm.stage}-${asm.id}"
            record(
                id       : build,
                meta     : m + record(assembly: [ assembler: asm.assembler, stage: asm.stage, id: asm.id, build: build ]),
                pri_fasta: file(asm.pri_fasta as String, checkIfExists: true),
                alt_fasta: asm.alt_fasta ? file(asm.alt_fasta as String, checkIfExists: true) : null,
                pri_gfa  : asm.pri_gfa ? file(asm.pri_gfa as String, checkIfExists: true) : null,
                alt_gfa  : asm.alt_gfa ? file(asm.alt_gfa as String, checkIfExists: true) : null
            )
        }.toList()
    }

    // Hi-C read pairs, numbered
    def hic_inputs = ( data.hic ?: [] ) as List
    hic = metadata.flatMap { m ->
        (0..<hic_inputs.size()).collect { index ->
            def pair = hic_inputs[index] as Map
            record(
                id   : m.id,
                meta : m + record(single_end: false, pair_id: index),
                reads: [ file(pair.read1 as String, checkIfExists: true), file(pair.read2 as String, checkIfExists: true) ]
            )
        }.toList()
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

    // RNAseq reads: paired or single end
    def rnaseq_inputs = ( data.rnaseq ?: [] ) as List
    rnaseq = metadata.flatMap { m ->
        rnaseq_inputs.collect { r ->
            def rec = r as Map
            def files = rec.reads
                ? [ file(rec.reads as String, checkIfExists: true) ]
                : [ file(rec.read1 as String, checkIfExists: true), file(rec.read2 as String, checkIfExists: true) ]
            record(id: m.id, meta: m + record(single_end: files.size() == 1), reads: files)
        }.toList()
    }

    // Isoseq reads
    def isoseq_inputs = ( data.isoseq ?: [] ) as List
    isoseq = metadata.flatMap { m ->
        isoseq_inputs.collect { r -> record(id: m.id, meta: m + record(single_end: true), reads: [ file((r as Map).reads as String, checkIfExists: true) ]) }.toList()
    }

    emit:
    goat       : Channel<Record>   = goat
    converted  : Channel<Record>   = converted_bams
    metadata   : Channel<Record>   = metadata
    assemblies : Channel<Assembly> = assemblies
    hic        : Channel<Reads>    = hic
    hifi       : Channel<Reads>    = hifi
    hifi_merged: Channel<Reads>    = hifi_merged
    rnaseq     : Channel<Reads>    = rnaseq
    isoseq     : Channel<Reads>    = isoseq
}
