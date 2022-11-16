#! /usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

include { SAMTOOLS_FASTA } from "$projectDir/modules/local/samtools/fasta/main"

/* params.input example sample sheet (samplesheet.yml)
```yaml
sample:
  id: Awesome_Species
  kmer_size: 31
  ploidy: 2
assembly:
  - id: assemblerX_build1
    pri_fasta: /path/to/primary_asm.fasta
    alt_fasta: /path/to/alternate_asm.fasta
    pri_gfa: /path/to/primary_asm.gfa
    alt_gfa: /path/to/alternate_asm.gfa
  - id: assemblerY_build1
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
[
    sample:[
        id: Awesome_Species
        kmer_size: 31
        ploidy: 2
    ],
    assembly:[
        [
            id:assemblerX_build1,
            pri_fasta: /path/to/primary_asm.fasta,
            alt_fasta: /path/to/alternate_asm.fasta,
            pri_gfa  : /path/to/primary_asm.gfa,
            alt_gfa  : /path/to/alternate_asm.gfa
        ],
        [
            id:assemblerY_build1,
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
*/

workflow PREPARE_INPUT {

    take:
    infile

    main:
    // Read in YAML
    Channel.fromPath( infile )
        .map { file -> readYAML( file ) }
        .multiMap { data ->
            assembly_ch : ( data.assembly ? [ data.sample, data.assembly ] : [] )
            hic_ch      : ( data.hic      ? [ data.sample + [ single_end: false ], data.hic.collect { [ file( it.read1, checkIfExists: true ), file( it.read2, checkIfExists: true ) ] } ] : [] )
            hifi_ch     : ( data.hifi     ? [ data.sample + [ single_end: true ], data.hifi.collect { file( it.reads, checkIfExists: true ) } ] : [] )
            rnaseq_ch   : ( data.rnaseq   ? [ data.sample, data.rnaseq.collect { it.reads ? file( it.reads, checkIfExists: true ) : [ file( it.read1, checkIfExists: true ), file( it.read2, checkIfExists: true ) ] } ] : [] )
            isoseq_ch   : ( data.isoseq   ? [ data.sample + [ single_end: true ], data.isoseq.collect { file( it.reads, checkIfExists: true ) } ] : [] )
        }
        .set{ input }

    // Convert assembly filename to files for correct staging
    input.assembly_ch
        .filter { !it.isEmpty() }
        .transpose()     // Data is [ sample, [ id:'assemblerX_build1', pri_fasta: '/path/to/primary_asm.fasta', alt_fasta: '/path/to/alternate_asm.fasta', pri_gfa: '/path/to/primary_asm.gfa', alt_gfa: '/path/to/alternate_asm.gfa' ]]
        .map { sample, assembly -> 
            [
                sample,
                [
                    id: assembly.id,
                    pri_fasta: file( assembly.pri_fasta, checkIfExists: true ),
                    alt_fasta: ( assembly.alt_fasta ? file( assembly.alt_fasta, checkIfExists: true ) : null ),
                    pri_gfa  : ( assembly.pri_gfa ? file( assembly.pri_gfa, checkIfExists: true ) : null ),
                    alt_gfa  : ( assembly.alt_gfa ? file( assembly.alt_gfa, checkIfExists: true ) : null )
                ]
            ]
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
