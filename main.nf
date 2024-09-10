#! /usr/bin/env nextflow

// Include Map.deepMerge() function
evaluate(new File("$projectDir/lib/MapExtended.groovy"))

include { combineByMetaKeys                        } from "$projectDir/modules/local/functions"
include { assembliesFromStage as preassembledInput } from "$projectDir/modules/local/functions"
include { setAssemblyStage                         } from "$projectDir/modules/local/functions"

include { PREPARE_INPUT } from "$projectDir/subworkflows/local/prepare_input/main"

include { BUILD_FASTK_DATABASE as BUILD_FASTK_HIFI_DATABASE } from "$projectDir/subworkflows/local/build_fastk_database/main"
include { BUILD_FASTK_DATABASE as BUILD_FASTK_HIC_DATABASE  } from "$projectDir/subworkflows/local/build_fastk_database/main"
include { BUILD_MERYL_DATABASE as BUILD_MERYL_HIFI_DATABASE } from "$projectDir/subworkflows/local/build_meryl_database/main"
include { BUILD_MERYL_DATABASE as BUILD_MERYL_HIC_DATABASE  } from "$projectDir/subworkflows/local/build_meryl_database/main"

include { INSPECT_DATA } from "$projectDir/subworkflows/local/inspect_data/main"

include { ASSEMBLE                                   } from "$projectDir/subworkflows/local/assemble/main"
include { ASSEMBLE_ORGANELLES                        } from "$projectDir/subworkflows/local/assemble_organelles/main"
include { COMPARE_ASSEMBLIES                         } from "$projectDir/subworkflows/local/compare_assemblies/main"
include { EVALUATE_ASSEMBLY as EVALUATE_RAW_ASSEMBLY } from "$projectDir/subworkflows/local/evaluate_assembly/main"

include { DECONTAMINATE } from "$projectDir/subworkflows/local/decontaminate/main"

include { PURGE_DUPLICATES } from "$projectDir/subworkflows/local/purge_dups/main"

include { EVALUATE_ASSEMBLY as EVALUATE_PURGED_ASSEMBLY } from "$projectDir/subworkflows/local/evaluate_assembly/main"

include { SCAFFOLD } from "$projectDir/subworkflows/local/scaffold/main.nf"
include { EVALUATE_ASSEMBLY as EVALUATE_SCAFFOLDED_ASSEMBLY } from "$projectDir/subworkflows/local/evaluate_assembly/main"

include { ALIGN_RNASEQ       } from "$projectDir/subworkflows/local/align_rnaseq/main"

include { ASSEMBLY_REPORT } from "$projectDir/subworkflows/local/assembly_report/main"

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

    // Setup sink channels
    ch_multiqc_files = Channel.value( file(params.multiqc_assembly_report_config, checkIfExists: true) )
    ch_quarto_files  = Channel.empty()
    ch_versions      = Channel.empty()

    // Read in data
    PREPARE_INPUT (
        params.input,
        params.ncbi.taxdb
    )

    // Build necessary databases
    if ( ['inspect','preprocess','assemble','purge','polish','screen','scaffold','curate'].any{ it in workflow_steps}) {
        // TODO: Migrate back to Meryl. Genome inspection missing KATGC and PLOIDYPLOT for meryldb
        BUILD_FASTK_HIFI_DATABASE ( PREPARE_INPUT.out.hifi )
        BUILD_FASTK_HIC_DATABASE ( PREPARE_INPUT.out.hic )
        BUILD_MERYL_HIFI_DATABASE ( PREPARE_INPUT.out.hifi )
        BUILD_MERYL_HIC_DATABASE ( PREPARE_INPUT.out.hic )
    }

    // Data inspection
    ch_hifi = PREPARE_INPUT.out.hifi
    if ( 'inspect' in workflow_steps ) {
        // QC Steps
        INSPECT_DATA(
            PREPARE_INPUT.out.hifi,
            BUILD_FASTK_HIFI_DATABASE.out.fastk_hist_ktab,
            BUILD_FASTK_HIC_DATABASE.out.fastk_hist_ktab
        )
        ch_hifi = INSPECT_DATA.out.hifi // with added kmer coverage
        ch_multiqc_files = ch_multiqc_files.mix( INSPECT_DATA.out.logs )
        ch_quarto_files = ch_quarto_files.mix( INSPECT_DATA.out.quarto_files )
        ch_versions = ch_versions.mix( INSPECT_DATA.out.versions )
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
    ch_raw_assemblies.dump(tag: 'Assemblies: Raw')

    // Organelle assembly
    if ( params.organelle_assembly_mode == 'reads' ) {
        // TODO: Add organelle assembly from reads
    } else if ( params.organelle_assembly_mode == 'contigs' ){
        ASSEMBLE_ORGANELLES ( ch_raw_assemblies )
        // TODO: filter organelles from assemblies
    } // else params.organelle_assembly_mode == 'none'

    // Assess assemblies
    COMPARE_ASSEMBLIES ( ch_raw_assemblies )
    EVALUATE_RAW_ASSEMBLY (
        ch_raw_assemblies,
        BUILD_FASTK_HIFI_DATABASE.out.fastk_hist_ktab,
        BUILD_MERYL_HIFI_DATABASE.out.uniondb
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        EVALUATE_RAW_ASSEMBLY.out.logs,
        COMPARE_ASSEMBLIES.out.logs
    )
    ch_versions = ch_versions.mix( EVALUATE_RAW_ASSEMBLY.out.versions )

    // Contamination screen
    ch_to_screen = setAssemblyStage (
        ch_raw_assemblies,
        'decontaminated' // Set assembly stage now for filenaming
    ).dump(tag: 'Assemblies: to screen')
    if ( 'screen' in workflow_steps ) {
        DECONTAMINATE( ch_to_screen )
        ch_cleaned_assemblies = DECONTAMINATE.out.assemblies
    } else {
        ch_cleaned_assemblies = ch_to_screen
    }
    ch_cleaned_assemblies = ch_cleaned_assemblies.mix (
        preassembledInput( PREPARE_INPUT.out.assemblies, 'decontaminated' )
    ).dump(tag: 'Assemblies: Cleaned')

    // Purge duplicates
    ch_to_purge = setAssemblyStage (
        ch_cleaned_assemblies,
        'purged' // Set assembly stage now for filenaming
    ).dump(tag: 'Assemblies: to purge')
    if ( 'purge' in workflow_steps ) {
        PURGE_DUPLICATES (
            ch_to_purge,
            ch_hifi
        )
        ch_purged_assemblies = PURGE_DUPLICATES.out.assemblies
    } else {
        ch_purged_assemblies = ch_to_purge
    }
    ch_purged_assemblies = ch_purged_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'purged' )
    ).dump(tag: 'Assemblies: Purged')
    EVALUATE_PURGED_ASSEMBLY (
        ch_purged_assemblies,
        BUILD_FASTK_HIFI_DATABASE.out.fastk_hist_ktab,
        BUILD_MERYL_HIFI_DATABASE.out.uniondb
    )
    ch_multiqc_files = ch_multiqc_files.mix( EVALUATE_PURGED_ASSEMBLY.out.logs )
    ch_versions = ch_versions.mix( EVALUATE_PURGED_ASSEMBLY.out.versions )

    // Polish
    ch_to_polish = setAssemblyStage (
        ch_purged_assemblies,
        'polished' // Set assembly stage now for filenaming
    ).dump(tag: 'Assemblies: to polish')
    if ( 'polish' in workflow_steps ) {
        // Run polishers
        ch_polished_assemblies = ch_to_polish
    } else {
        ch_polished_assemblies = ch_to_polish
    }
    ch_polished_assemblies = ch_polished_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'polished' )
    ).dump(tag: 'Assemblies: Polished')

    // Scaffold
    ch_to_scaffold = setAssemblyStage (
        ch_polished_assemblies,
        'scaffolded' // Set assembly stage now for filenaming
    ).dump(tag: 'Assemblies: to scaffold')
    if ( 'scaffold' in workflow_steps ) {
        SCAFFOLD (
            ch_to_scaffold,
            PREPARE_INPUT.out.hic
        )
        ch_scaffolded_assemblies = SCAFFOLD.out.assemblies
    } else {
        ch_scaffolded_assemblies = ch_to_scaffold
    }
    ch_scaffolded_assemblies = ch_scaffolded_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'scaffolded' )
    ).dump(tag: 'Assemblies: Scaffolded')
    EVALUATE_SCAFFOLDED_ASSEMBLY (
        ch_scaffolded_assemblies,
        BUILD_FASTK_HIFI_DATABASE.out.fastk_hist_ktab,
        BUILD_MERYL_HIFI_DATABASE.out.uniondb
    )
    ch_multiqc_files = ch_multiqc_files.mix( EVALUATE_SCAFFOLDED_ASSEMBLY.out.logs )
    ch_versions = ch_versions.mix( EVALUATE_SCAFFOLDED_ASSEMBLY.out.versions )

    // Curate
    ch_to_curate = setAssemblyStage (
        ch_scaffolded_assemblies,
        'curated' // Set assembly stage now for filenaming
    ).dump(tag: 'Assemblies: to curate')
    if ( 'curate' in workflow_steps ) {
        // Run assemblers
        ch_curated_assemblies = ch_to_curate
    } else {
        ch_curated_assemblies = ch_to_curate
    }
    ch_curated_assemblies = ch_curated_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'curated' )
    ).dump(tag: 'Assemblies: Curated')

    // Align RNAseq
    if( 'alignRNA' in workflow_steps ) {
        ALIGN_RNASEQ (
            PREPARE_INPUT.out.rnaseq,
            PREPARE_INPUT.out.assemblies // TODO: Select assembly stage
                .map { meta, assembly -> [ meta, assembly.pri_fasta ] }
        )
    }

    ASSEMBLY_REPORT(
        PREPARE_INPUT.out.sample_meta.map{ meta -> [ meta, file(params.quarto_assembly_report, checkIfExists: true) ] },
        ch_multiqc_files,
        ch_versions
    )
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
