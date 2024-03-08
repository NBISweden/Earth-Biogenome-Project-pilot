#! /usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

include { UNTAR as UNTAR_TAXONOMY } from "$projectDir/modules/nf-core/untar/main"
include { TAXONKIT_NAME2LINEAGE   } from "$projectDir/modules/local/taxonkit/name2lineage"
include { GOAT_TAXONSEARCH        } from "$projectDir/modules/nf-core/goat/taxonsearch/main"
include { SAMTOOLS_FASTA          } from "$projectDir/modules/local/samtools/fasta/main"
include { CAT_CAT as MERGE_PACBIO } from "$projectDir/modules/nf-core/cat/cat/main"

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
    taxdb

    main:
    // Read in YAML
    ch_input = Channel.fromPath( infile, checkIfExists: true )
        .map { file -> readYAML( file ) }

    UNTAR_TAXONOMY( Channel.fromPath( taxdb, checkIfExists: true ).map{ tar -> [ [ id: 'taxdb' ], tar ] } )
    TAXONKIT_NAME2LINEAGE( ch_input, UNTAR_TAXONOMY.out.untar.map{ meta, archive -> archive }.collect() ).tsv
        .branch { meta, tsv_f -> def sv = tsv_f.splitCsv( sep:"\t" )
        def new_meta = meta.deepMerge( [ id: sv[0][0].replace(" ","_"), sample: [ taxid: sv[0][1], kingdom: sv[0][2] ] ] )
            eukaryota: sv[0][2] == 'Eukaryota'
                return new_meta
            other: true
                return new_meta
        }.set { ch_input_wTaxID }
    // Update meta with GOAT before propagating (eukaryotes only)
    GOAT_TAXONSEARCH( ch_input_wTaxID.eukaryota.map { data -> [ data, data.sample.name, [] ] } ).taxonsearch
        .map { meta, tsv ->
            def busco_lineages = tsv.splitCsv( sep:"\t", header: true ).findAll { it.odb10_lineage }.collect { it.odb10_lineage }.join(',')
            def species = tsv.splitCsv( sep:"\t", header: true ).find { it.scientific_name == meta.sample.name }
            meta.deepMerge([
                sample: [
                    genome_size: meta.sample.genome_size ?: species.genome_size,
                    haploid_number: meta.sample.haploid_number ?: species.haploid_number,
                    ploidy: meta.sample.ploidy ?: species.ploidy
                ],
                settings: [ busco: [ lineages: params.busco.lineages?: busco_lineages ] ]
            ])
        }
        .mix( ch_input_wTaxID.other )
        .tap { sample_meta_ch }
        .dump( tag: 'Input: Meta', pretty: true )
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
                build: "${assembly.assembler}-${assembly.stage}-${assembly.id}",
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
    sample_fastx = hifi_fastx_ch.groupTuple()
        .branch { meta, fastx ->
            single: fastx.size() == 1
                return [ meta, *fastx ]
            multi: true
                return [ meta, fastx ]
        }
    MERGE_PACBIO( sample_fastx.multi )

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
    sample_meta = sample_meta_ch.map { data -> data.subMap('sample') }
    assemblies  = assembly_ch.dump( tag: 'Input: Assemblies' )
    hic         = hic_fastx_ch.dump( tag: 'Input: Hi-C' )
    hifi        = hifi_fastx_ch.dump( tag: 'Input: PacBio HiFi' )
    hifi_merged = sample_fastx.single.mix( MERGE_PACBIO.out.file_out )
    rnaseq      = rnaseq_fastx_ch.dump( tag: 'Input: Illumina RnaSeq' )
    isoseq      = isoseq_fastx_ch.dump( tag: 'Input: PacBio IsoSeq' )
}

def readYAML( yamlfile ) {
    // TODO: Validate sample file
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}
