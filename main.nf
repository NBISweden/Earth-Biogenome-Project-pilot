#! /usr/bin/env nextflow

// Include Map.deepMerge() function
evaluate(new File("$projectDir/lib/MapExtended.groovy"))

include { combineByMetaKeys                        } from "$projectDir/modules/local/functions"
include { assembliesFromStage as preassembledInput } from "$projectDir/modules/local/functions"
include { setAssemblyStage                         } from "$projectDir/modules/local/functions"

include { PREPARE_INPUT } from "$projectDir/subworkflows/local/prepare_input/main"

include { BUILD_DATABASES as BUILD_HIFI_DATABASES } from "$projectDir/subworkflows/local/build_databases/main"
include { BUILD_DATABASES as BUILD_HIC_DATABASES  } from "$projectDir/subworkflows/local/build_databases/main"

include { INSPECT_DATA } from "$projectDir/subworkflows/local/inspect_data/main"

include { ASSEMBLE                                   } from "$projectDir/subworkflows/local/assemble/main"
include { ASSEMBLE_ORGANELLES                        } from "$projectDir/subworkflows/local/assemble_organelles/main"
include { COMPARE_ASSEMBLIES                         } from "$projectDir/subworkflows/local/compare_assemblies/main"
include { EVALUATE_ASSEMBLY as EVALUATE_RAW_ASSEMBLY } from "$projectDir/subworkflows/local/evaluate_assembly/main"

include { DECONTAMINATE } from "$projectDir/modules/local/decontaminate/main"

include { PURGE_DUPLICATES } from "$projectDir/subworkflows/local/purge_dups/main"

include { EVALUATE_ASSEMBLY as EVALUATE_PURGED_ASSEMBLY } from "$projectDir/subworkflows/local/evaluate_assembly/main"

include { ALIGN_RNASEQ       } from "$projectDir/subworkflows/local/align_rnaseq/main"

/*
 * Development: See docs/development to understand the workflow programming model and
 * how channel contents are structured.
 */

workflow {

    // Define constants
    def workflow_permitted_stages = [
        'inspect',      // 01 - Read inspection
        'preprocess',   // 02 - Read preprocessing
        'assemble',     // 03 - Assembly
        'screen',       // 04 - Contamination screening
        'purge',        // 05 - Duplicate purging
        'polish',       // 06 - Error polishing
        'scaffold',     // 07 - Scaffolding
        'curate',       // 08 - Rapid curation
        'alignRNA'      // 09 - Align RNAseq data
    ]

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
    PREPARE_INPUT (
        params.input,
        params.ncbi.taxdb
    )

    // Build necessary databases
    if ( ['inspect','preprocess','assemble','purge','polish','screen','scaffold','curate'].any{ it in workflow_steps}) {
        BUILD_HIFI_DATABASES ( PREPARE_INPUT.out.hifi )
        BUILD_HIC_DATABASES ( PREPARE_INPUT.out.hic )
    }

    // Data inspection
    if ( 'inspect' in workflow_steps ) {
        // QC Steps
        INSPECT_DATA(
            PREPARE_INPUT.out.hifi,
            BUILD_HIFI_DATABASES.out.fastk_histogram.join( BUILD_HIFI_DATABASES.out.fastk_ktab ),
            BUILD_HIC_DATABASES.out.fastk_histogram.join( BUILD_HIC_DATABASES.out.fastk_ktab )
        )
    }

    // Preprocess data
    if ( 'preprocess' in workflow_steps ) {
        // Adapter filtering etc
        // Subsampling
        // Host contamination filter
    }

    // Assemble
    ch_raw_assemblies = preassembledInput( PREPARE_INPUT.out.assemblies, 'raw' )
    if ( 'assemble' in workflow_steps ) {
        // Run assemblers
        ASSEMBLE ( PREPARE_INPUT.out.hifi_merged )
        ch_raw_assemblies = ch_raw_assemblies.mix( ASSEMBLE.out.raw_assemblies )
    } else {
        // Nothing more than evaluate
    }
    // TODO: Add organelle assembly from reads
    ASSEMBLE_ORGANELLES ( raw_assemblies )
    // TODO: filter organelles from assemblies

    // Assess assemblies
    COMPARE_ASSEMBLIES ( ch_raw_assemblies )
    EVALUATE_RAW_ASSEMBLY (
        ch_raw_assemblies,
        PREPARE_INPUT.out.hifi,
        BUILD_HIFI_DATABASES.out.fastk_histogram.join( BUILD_HIFI_DATABASES.out.fastk_ktab ),
    )

    // Contamination screen
    ch_to_screen = setAssemblyStage (
        ch_raw_assemblies,
        'decontaminated' // Set assembly stage now for filenaming
    )
    if ( 'screen' in workflow_steps ) {
        DECONTAMINATE( ch_to_screen )
        ch_cleaned_assemblies = DECONTAMINATE.out.assemblies
    } else {
        ch_cleaned_assemblies = ch_to_screen
    }
    ch_cleaned_assemblies = ch_cleaned_assemblies.mix (
        preassembledInput( PREPARE_INPUT.out.assemblies, 'decontaminated' )
    )

    // Purge duplicates
    ch_to_purge = setAssemblyStage (
        ch_cleaned_assemblies,
        'purged' // Set assembly stage now for filenaming
    )
    if ( 'purge' in workflow_steps ) {
        PURGE_DUPLICATES (
            ch_to_purge,
            'inspect' in workflow_steps ? INSPECT_DATA.out.hifi : PREPARE_INPUT.out.hifi
        )
        ch_purged_assemblies = PURGE_DUPLICATES.out.assemblies
    } else {
        ch_purged_assemblies = ch_to_purge
    }
    ch_purged_assemblies = ch_purged_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'purged' )
    )
    EVALUATE_PURGED_ASSEMBLY (
        ch_purged_assemblies,
        PREPARE_INPUT.out.hifi,
        BUILD_HIFI_DATABASES.out.fastk_histogram.join( BUILD_HIFI_DATABASES.out.fastk_ktab )
    )

    // Polish
    ch_to_polish = setAssemblyStage (
        ch_purged_assemblies,
        'polished' // Set assembly stage now for filenaming
    )
    if ( 'polish' in workflow_steps ) {
        // Run polishers
        ch_polished_assemblies = ch_to_polish
    } else {
        ch_polished_assemblies = ch_to_polish
    }
    ch_polished_assemblies = ch_polished_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'polished' )
    )

    // Scaffold
    ch_to_scaffold = setAssemblyStage (
        ch_polished_assemblies,
        'scaffolded' // Set assembly stage now for filenaming
    )
    if ( 'scaffold' in workflow_steps ) {
        // Run scaffolder
        ch_scaffolded_assemblies = ch_to_scaffold
    } else {
        ch_scaffolded_assemblies = ch_to_scaffold
    }
    ch_scaffolded_assemblies = ch_scaffolded_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'scaffolded' )
    )

    // Curate
    ch_to_curate = setAssemblyStage (
        ch_scaffolded_assemblies,
        'curated' // Set assembly stage now for filenaming
    )
    if ( 'curate' in workflow_steps ) {
        // Run assemblers
        ch_curated_assemblies = ch_to_curate
    } else {
        ch_curated_assemblies = ch_to_curate
    }
    ch_curated_assemblies = ch_curated_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'curated' )
    )

    // Align RNAseq
    if( 'alignRNA' in workflow_steps ) {
        ALIGN_RNASEQ (
            PREPARE_INPUT.out.rnaseq,
            PREPARE_INPUT.out.assemblies // TODO: Select assembly stage
                .map { meta, assembly -> [ meta, assembly.pri_fasta ] }
        )
    }
}

workflow.onComplete {
    if( workflow.success ){
        log.info("""
        Thank you for using the NBIS Earth Biogenome Project Assembly workflow.
        The workflow completed successfully.

        Results are located in the folder: $params.outdir
        """)
    } else {
        log.info("""
        Thank you for using the NBIS Earth Biogenome Project Assembly workflow.
        The workflow completed unsuccessfully.

        Please read over the error message. If you are unable to solve it, please
        post an issue at https://github.com/NBISweden/Earth-Biogenome-Project-pilot/issues
        where we will do our best to help.
        """)
    }
}
