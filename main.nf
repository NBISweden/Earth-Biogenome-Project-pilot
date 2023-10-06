#! /usr/bin/env nextflow

include { PREPARE_INPUT } from "$projectDir/subworkflows/prepare_input/main"

include { BUILD_DATABASES as BUILD_HIFI_DATABASES } from "$projectDir/subworkflows/build_databases/main"
include { BUILD_DATABASES as BUILD_HIC_DATABASES  } from "$projectDir/subworkflows/build_databases/main"

include { GENOME_PROPERTIES } from "$projectDir/subworkflows/genome_properties/main"
include { COMPARE_LIBRARIES } from "$projectDir/subworkflows/compare_libraries/main"
include { SCREEN_READS      } from "$projectDir/subworkflows/screen_read_contamination/main"

include { PURGE_DUPLICATES } from "$projectDir/subworkflows/purge_dups/main"

include { COMPARE_ASSEMBLIES } from "$projectDir/subworkflows/compare_assemblies/main"
include { EVALUATE_ASSEMBLY  } from "$projectDir/subworkflows/evaluate_assembly/main"
include { ALIGN_RNASEQ       } from "$projectDir/subworkflows/align_rnaseq/main"

workflow {

    // Define constants
    def workflow_permitted_stages = ['data_qc','preprocess','assemble','validate','curate']

    // Check input
    def workflow_steps = params.steps.tokenize(",")
    if ( ! workflow_steps.every { it in workflow_permitted_stages } ) {
        error "Unrecognised workflow step in $params.steps ( $workflow_permitted_stages )"
    }

    // The primary workflow for the Earth Biogenome Project Pilot
    log.info("""
    Running NBIS Earth Biogenome Project Assembly workflow.
    """)

    // Read in data
    PREPARE_INPUT ( params.input )

    // Build necessary databases
    if ( ['data_qc','validate'].any{ it in workflow_steps}) {
        BUILD_HIFI_DATABASES ( PREPARE_INPUT.out.hifi )
        BUILD_HIC_DATABASES ( PREPARE_INPUT.out.hic )
    }
    
    // Data inspection
    if ( 'data_qc' in workflow_steps ) {
        // QC Steps
        GENOME_PROPERTIES ( 
            BUILD_HIFI_DATABASES.out.fastk_histogram.join( BUILD_HIFI_DATABASES.out.fastk_ktab )
            // BUILD_HIFI_DATABASES.out.meryl_histogram
        )
        GENOME_PROPERTIES.out.model
            .map { meta, file ->
                // Parse kmer coverage from GenomeScope model
                def covline = file.readLines().collect { it.startsWith('kmercov') ? it : '' }
                [ meta , [ kmercov: new BigDecimal( covline.join('').tokenize(' ')[1] ).round(2) ] ]
            }.set { ch_hifi_kmercov }
        COMPARE_LIBRARIES (
            BUILD_HIFI_DATABASES.out.fastk_histogram.join( BUILD_HIFI_DATABASES.out.fastk_ktab ).join(
            BUILD_HIC_DATABASES.out.fastk_histogram.join( BUILD_HIC_DATABASES.out.fastk_ktab ) )
        )
        SCREEN_READS ( 
            PREPARE_INPUT.out.hifi,
            // TODO:: Allow custom database ala nf-core/genomeassembler.
            file( params.mash_screen_db, checkIfExists: true )
        )
    }


    // Preprocess data
    if ( 'preprocess' in workflow_steps ) {
        // Adapter filtering etc
    }
    
    // Assemble
    if ( 'assemble' in workflow_steps ) {
        // Run assemblers
    }

    // Curate assemblies 
    if ( 'curate' in workflow_steps ) {
        ch_topurge = PREPARE_INPUT.out.hifi
                .map { meta, reads -> [ meta.findAll { ! (it.key in [ 'single_end' ]) }, reads ] }
                .combine( PREPARE_INPUT.out.assemblies, by:0 )
        if ( 'data_qc' in workflow_steps ) {
            // Add kmer coverage from GenomeScope model
            ch_topurge.combine( ch_hifi_kmercov, by: 0 )
                .map { meta, reads, assemblies, kmer_cov -> [ meta + kmer_cov, reads, assemblies ] }
                .set { ch_topurge }
        }
        PURGE_DUPLICATES ( ch_topurge.dump( tag: 'Purge duplicates: input' ) )
        // Break and reassemble misassemblies, separate organelles, etc
            // MitoHiFi
            // PurgeDups
            // Kraken2
            // Blobtoolkit
            // FCS-Genome
    }

    // Assess assemblies
    if ( 'validate' in workflow_steps ) {
        COMPARE_ASSEMBLIES (
            PREPARE_INPUT.out.assemblies,
            params.reference ? file( params.reference, checkIfExists: true ) : []
        )

        EVALUATE_ASSEMBLY (
            PREPARE_INPUT.out.assemblies,
            PREPARE_INPUT.out.hifi,
            BUILD_HIFI_DATABASES.out.fastk_histogram.join( BUILD_HIFI_DATABASES.out.fastk_ktab ),
            params.reference ? file( params.reference, checkIfExists: true ) : [],
            Channel.of( params.busco_lineages.tokenize(',') ),
            params.busco_lineage_path ? file( params.busco_lineage_path, checkIfExists: true ) : []
        )

        // TODO: Run only if RNAseq data
        ALIGN_RNASEQ ( 
            PREPARE_INPUT.out.rnaseq,
            PREPARE_INPUT.out.assemblies
                .map { meta, assembly -> [ meta + [ build: assembly.id ], assembly.pri_fasta ] }
        )
    }

}

workflow.onComplete {
    if( workflow.success ){
        log.info("""
        Thank you for using the NBIS Earth Biogenome Project Assembly workflow.

        Results are located in the folder: $params.outdir
        """)
    } else {
        log.info("""
        The workflow completed unsuccessfully.

        Please read over the error message. If you are unable to solve it, please
        post an issue at https://github.com/NBISweden/Earth-Biogenome-Project-pilot/issues
        where we will do our best to help.
        """)
    }
}
