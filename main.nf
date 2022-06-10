#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPARE_INPUT } from "$projectDir/subworkflows/prepare_input/main"

include { BUILD_DATABASES as BUILD_HIFI_DATABASES } from "$projectDir/subworkflows/build_databases/main"
include { BUILD_DATABASES as BUILD_HIC_DATABASES  } from "$projectDir/subworkflows/build_databases/main"

include { GENOME_PROPERTIES } from "$projectDir/subworkflows/genome_properties/genome_properties"

include { ASSEMBLY_VALIDATION } from "$projectDir/subworkflows/assembly_validation/assembly_validation"

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
    BUILD_HIFI_DATABASES ( PREPARE_INPUT.out.hifi )
    BUILD_HIC_DATABASES ( PREPARE_INPUT.out.hic )
    
    // Data inspection
    if ( 'data_qc' in workflow_steps ) {
        // QC Steps
        GENOME_PROPERTIES( 
            BUILD_HIFI_DATABASES.out.fastk_histex,
            BUILD_HIFI_DATABASES.out.meryl_histogram
        )
    }

    // Preprocess data
    if ( 'preprocess' in workflow_steps ) {
        // Adapter filtering etc
    }
    
    // Assemble
    if( 'assemble' in workflow_steps ) {
        // Run assemblers
    }

    // Curate assemblies 
    if ( 'curate' in workflow_steps ) {
        // Break and reassemble misassemblies, separate organelles, etc
    }

    // Assess assemblies
    if ( 'validate' in workflow_steps ) {
        ASSEMBLY_VALIDATION(
            PREPARE_INPUT.out.assemblies,
            PREPARE_INPUT.out.hifi,
            BUILD_HIFI_DATABASES.out.fastk_histogram.join( BUILD_HIFI_DATABASES.out.fastk_ktab ),
            BUILD_HIFI_DATABASES.out.meryl_uniondb,
            params.reference ? file( params.reference, checkIfExists: true ) : [],
            params.busco_lineages.tokenize(','),
            params.busco_lineage_path ? file( params.busco_lineage_path, checkIfExists: true ) : [],
            Channel.fromPath( params.diamond_db, checkIfExists: true ),
            Channel.fromPath( params.blast_db, checkIfExists: true ),
            file( params.ncbi_taxonomy, checkIfExists: true )
        )
    }

}

// Deprecated workflow - TODO: remove workflow
// workflow VALIDATE_ASSEMBLIES {

//     log.info("""
//     Running NBIS Earth Biogenome Project Assembly validation workflow.
//     """)

//     PREPARE_INPUT ( params.input )
//     ASSEMBLY_VALIDATION(
//         PREPARE_INPUT.out.assemblies,
//         PREPARE_INPUT.out.hifi,
//         params.reference ? file( params.reference, checkIfExists: true ) : [],
//         params.busco_lineages.tokenize(','),
//         params.busco_lineage_path ? file( params.busco_lineage_path, checkIfExists: true ) : [],
//         file( params.uniprot_db, checkIfExists: true ),
//         file( params.ncbi_nt_db, checkIfExists: true ),
//         file( params.ncbi_taxonomy, checkIfExists: true )
//     )

// }

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
