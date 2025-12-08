#! /usr/bin/env nextflow

// include { UNTAR as UNTAR_TAXONOMY } from "../../../modules/nf-core/untar/main"
// include { TAXONKIT_NAME2LINEAGE   } from "../../../modules/local/taxonkit/name2lineage"
include { FETCH_SAMPLE_METADATA   } from "./fetch_sample_metadata"
// include { GOAT_TAXONSEARCH        } from "../../../modules/nf-core/goat/taxonsearch/main"
include { SAMTOOLS_FASTA          } from "../../../modules/local/samtools/fasta/main"
include { CAT_CAT as MERGE_PACBIO } from "../../../modules/nf-core/cat/cat/main"

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
mito_hmm:
  - fam: /path/to/mito_hmm.fam
    h3f: /path/to/mito_hmm.fam.h3f
    h3i: /path/to/mito_hmm.fam.h3i
    h3m: /path/to/mito_hmm.fam.h3m
    h3p: /path/to/mito_hmm.fam.h3p
plastid_hmm:
  - fam: /path/to/plastid_hmm.fam
    h3f: /path/to/plastid_hmm.fam.h3f
    h3i: /path/to/plastid_hmm.fam.h3i
    h3m: /path/to/plastid_hmm.fam.h3m
    h3p: /path/to/plastid_hmm.fam.h3p
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
    // taxdb

    main:
    // Read in YAML
    ch_input = channel.fromPath( infile, checkIfExists: true )
        .map { file -> readYAML( file ) }

    // Populate meta data
    FETCH_SAMPLE_METADATA( ch_input ).metadata
        .dump( tag: 'Input: Meta', pretty: true )
        .multiMap { data ->
            assembly_ch    : ( data.assembly ? [ data.subMap('id','sample','settings') , data.assembly ] : [] )
            hic_ch         : ( data.hic      ? [ data.subMap('id','sample','settings') + [ single_end: false ], data.hic.collect { [ file( it.read1, checkIfExists: true ), file( it.read2, checkIfExists: true ) ] } ] : [] )
            hifi_ch        : ( data.hifi     ? [ data.subMap('id','sample','settings') + [ single_end: true ], data.hifi.collect { file( it.reads, checkIfExists: true ) } ] : [] )
            rnaseq_ch      : ( data.rnaseq   ? [ data.subMap('id','sample','settings'), data.rnaseq.collect { it.reads ? file( it.reads, checkIfExists: true ) : [ file( it.read1, checkIfExists: true ), file( it.read2, checkIfExists: true ) ] } ] : [] )
            isoseq_ch      : ( data.isoseq   ? [ data.subMap('id','sample','settings') + [ single_end: true ], data.isoseq.collect { file( it.reads, checkIfExists: true ) } ] : [] )
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
        .flatMap { meta, hic_pairs -> hic_pairs.withIndex().collect{ pair, index -> [ meta + [pair_id: index], pair ] } }
        .set { hic_fastx_ch }

    // Prepare PacBio HiFi channel
    // Convert HiFi BAMS to FastQ
    input.hifi_ch
        .filter { !it.isEmpty() }
        .transpose()   // Transform to [ [ id: 'sample_name'], file('/path/to/read')  ]
        .branch { _meta, filename ->
            bam_ch: filename.toString().endsWith(".bam")
            fastx_ch: true // assume everything else is fastx
        }.set { hifi }
    SAMTOOLS_FASTA ( hifi.bam_ch )
    hifi.fastx_ch.mix( SAMTOOLS_FASTA.out.fasta )
        .set { hifi_fastx_ch }
    sample_fastx = hifi_fastx_ch.groupTuple()
        .branch { meta, fastx ->
            single: fastx.size() == 1
                return [ meta ] + fastx
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

    // versions
    FETCH_SAMPLE_METADATA.out.versions.mix(
        SAMTOOLS_FASTA.out.versions.first(),
        MERGE_PACBIO.out.versions.first(),
    ).set { versions }

    emit:
    sample_meta = FETCH_SAMPLE_METADATA.out.metadata.map { data -> data.subMap('sample') }
    assemblies  = assembly_ch.dump( tag: 'Input: Assemblies', pretty: true )
    hic         = hic_fastx_ch.dump( tag: 'Input: Hi-C', pretty: true )
    hifi        = hifi_fastx_ch.dump( tag: 'Input: PacBio HiFi', pretty: true )
    hifi_merged = sample_fastx.single.mix( MERGE_PACBIO.out.file_out )
    rnaseq      = rnaseq_fastx_ch.dump( tag: 'Input: Illumina RnaSeq', pretty: true )
    isoseq      = isoseq_fastx_ch.dump( tag: 'Input: PacBio IsoSeq', pretty: true )
    versions
}

def readYAML( yamlfile ) {
    // TODO: Validate sample file - nf-schema
    return new org.yaml.snakeyaml.Yaml().load( new FileReader( yamlfile.toString() ) )
}
