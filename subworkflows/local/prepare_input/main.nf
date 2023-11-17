#! /usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

include { SAMTOOLS_FASTA } from "$projectDir/modules/local/samtools/fasta/main"
include { GOAT_TAXONSEARCH } from "$projectDir/modules/nf-core/goat/taxonsearch/main"

/* params.input example sample sheet (samplesheet.yml)
```yaml
sample:
  name: Species name
assembly:
  - assembler: hifiasm
    stage: raw
    id: uuid
    pri_fasta: /path/to/primary_asm.fasta
    alt_fasta: /path/to/alternate_asm.fasta
    pri_gfa: /path/to/primary_asm.gfa
    alt_gfa: /path/to/alternate_asm.gfa
  - assembler: ipa
    stage: raw
    id: uuid
    pri_fasta: /path/to/primary_asm.fasta
    alt_fasta: /path/to/alternate_asm.fasta
hic:
  - read1: /path/to/read1.fastq.gz
    read2: /path/to/read2.fastq.gz
hifi:
  - reads: /path/to/reads
rnaseq:
  - reads: /path/to/reads
isoseq:
  - reads: /path/to/reads
```
leads to the following YAML data structure
```
[
    sample:[
        name: Species name
    ],
    assembly:[
        [
            assembler: hifiasm,
            stage: raw,
            id: uuid,
            pri_fasta: /path/to/primary_asm.fasta,
            alt_fasta: /path/to/alternate_asm.fasta,
            pri_gfa  : /path/to/primary_asm.gfa,
            alt_gfa  : /path/to/alternate_asm.gfa
        ],
        [
            assembler: ipa,
            stage: raw,
            id: uuid,
            pri_fasta: /path/to/primary_asm.fasta,
            alt_fasta: /path/to/alternate_asm.fasta
        ]
    ],
    hic:[
        read1: /path/to/read1.fastq.gz,
        read2: /path/to/read2.fastq.gz
    ],
    hifi:[
        reads: /path/to/reads
    ],
    rnaseq:[
        reads: /path/to/reads
    ],
    isoseq:[
        reads: /path/to/reads
    ]
]
```
*/

/*
Output meta map structure:
```
[
    id: "Species_name", // Needed for nf-core modules functionality
    single_end: true,   // Needed for nf-core modules functionality
    sample: [
        name: "Species name",
        genome_size: 1234,     // GOAT
        ploidy: 2,             // GOAT
        chr_count: 13          // GOAT
    ],
    reads_hifi: [
        single_end: true,
        kmer_cov: 25
    ],
    reads_hic: [
        single_end: false,
        kmer_cov: 50
    ]
    settings: [
        busco: [
            lineages: auto     // GOAT
        ]
    ]
    assembly: [
        assembler: "hifiasm",
        stage: "raw",
        id: uuid
    ]
]
```
*/

workflow PREPARE_INPUT {

    take:
    infile

    main:
    // Read in YAML
    ch_input = Channel.fromPath( infile )
        .map { file -> readYAML( file ) }
    
    // Update meta with GOAT before propagating
    // TODO: Bypass if certain fields are present
    GOAT_TAXONSEARCH( ch_input.map { data -> [ data, data.sample.name, [] ] } ).taxonsearch
        .map { meta, tsv ->
            def busco_lineages = tsv.splitCsv( sep:"\t", header: true ).findAll { it.odb10_lineage }.collect { it.odb10_lineage }.join(',') 
            def species = tsv.splitCsv( sep:"\t", header: true ).find { it.scientific_name == meta.sample.name }
            // Update meta in place since there should be no concurrent access here.
            meta.id = meta.sample.name.replace(" ","_")
            meta.sample = meta.sample + [ genome_size: species.genome_size, chromosome_number: species.chromosome_number, ploidy: species.ploidy ]
            if( ! meta.settings ) {
                meta = meta + [ settings: [ busco: [ lineages: busco_lineages ] ] ]
            } else if ( ! meta.settings.busco ) {
                meta.settings = meta.settings + [ busco: [ lineages: busco_lineages ] ]
            } else if ( ! meta.settings.busco.lineages ) {
                meta.settings.busco = meta.settings.busco + [ lineages: busco_lineages ]
            } else {
                meta // Leave settings unchanged
            }
            meta
        }
        .multiMap { data ->
            assembly_ch : ( data.assembly ? [ data.subMap('id','sample','settings') , data.assembly ] : [] )
            hic_ch      : ( data.hic      ? [ data.subMap('id','sample','settings') + [ single_end: false ], data.hic.collect { [ file( it.read1, checkIfExists: true ), file( it.read2, checkIfExists: true ) ] } ] : [] )
            hifi_ch     : ( data.hifi     ? [ data.subMap('id','sample','settings') + [ single_end: true ], data.hifi.collect { file( it.reads, checkIfExists: true ) } ] : [] )
            rnaseq_ch   : ( data.rnaseq   ? [ data.subMap('id','sample','settings'), data.rnaseq.collect { it.reads ? file( it.reads, checkIfExists: true ) : [ file( it.read1, checkIfExists: true ), file( it.read2, checkIfExists: true ) ] } ] : [] )
            isoseq_ch   : ( data.isoseq   ? [ data.subMap('id','sample','settings') + [ single_end: true ], data.isoseq.collect { file( it.reads, checkIfExists: true ) } ] : [] )
        }
        .set{ input }

    // Convert assembly filename to files for correct staging
    input.assembly_ch
        .filter { !it.isEmpty() }
        .transpose()     // Data is now [ sample, [ id:'assemblerX_build1', pri_fasta: '/path/to/primary_asm.fasta', alt_fasta: '/path/to/alternate_asm.fasta', pri_gfa: '/path/to/primary_asm.gfa', alt_gfa: '/path/to/alternate_asm.gfa' ]]
        .map { meta, assembly ->
            def asm = [
                assembler: assembly.assembler,
                stage: assembly.stage,
                id: assembly.id,
                build: "$assembly.assembler-$assembly.stage-$assembly.id",
                pri_fasta: file( assembly.pri_fasta, checkIfExists: true ),
                alt_fasta: ( assembly.alt_fasta ? file( assembly.alt_fasta, checkIfExists: true ) : null ),
                pri_gfa  : ( assembly.pri_gfa ? file( assembly.pri_gfa, checkIfExists: true ) : null ),
                alt_gfa  : ( assembly.alt_gfa ? file( assembly.alt_gfa, checkIfExists: true ) : null )
            ]
            [ meta + [ assembly: asm.subMap(['assembler','stage','id','build']) ], asm ]
        }
        .set { assembly_ch }

    // Combine Hi-C channels
    input.hic_ch.filter { !it.isEmpty() }
        .transpose()
        .set { hic_fastx_ch }

    // Prepare PacBio HiFi channel
    // Convert HiFi BAMS to FastQ
    input.hifi_ch
        .filter { !it.isEmpty() }
        .transpose()   // Transform to [ [ id: 'sample_name'], file('/path/to/read')  ]
        .branch { meta, filename ->
            bam_ch: filename.toString().endsWith(".bam")
            fastx_ch: true // assume everything else is fastx
        }.set { hifi }
    SAMTOOLS_FASTA ( hifi.bam_ch )
    hifi.fastx_ch.mix( SAMTOOLS_FASTA.out.fasta )
        .set { hifi_fastx_ch }

    // Prepare RNAseq channel
    input.rnaseq_ch.filter { !it.isEmpty() }
        .transpose()
        .map { meta, reads -> [ meta + [ single_end: reads instanceof Path ], reads ] }
        .set { rnaseq_fastx_ch }

    // Prepare Isoseq channel
    input.isoseq_ch.filter { !it.isEmpty() }
        .transpose()
        .set { isoseq_fastx_ch }

    emit:
    assemblies = assembly_ch.dump( tag: 'Input: Assemblies' )
    hic        = hic_fastx_ch.dump( tag: 'Input: Hi-C' )
    hifi       = hifi_fastx_ch.dump( tag: 'Input: PacBio HiFi' )
    rnaseq     = rnaseq_fastx_ch.dump( tag: 'Input: Illumina RnaSeq' )
    isoseq     = isoseq_fastx_ch.dump( tag: 'Input: PacBio IsoSeq' )
}

def readYAML( yamlfile ) {
    // TODO: Validate sample file
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}
