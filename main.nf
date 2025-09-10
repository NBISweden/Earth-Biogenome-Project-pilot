#! /usr/bin/env nextflow

// Functions
include { combineByMetaKeys                                 } from "./modules/local/functions"
include { assembliesFromStage as preassembledInput          } from "./modules/local/functions"
include { setAssemblyStage                                  } from "./modules/local/functions"
// Data import
include { PREPARE_INPUT                                     } from "./subworkflows/local/prepare_input/main"
// Database preparation
include { BUILD_FASTK_DATABASE as BUILD_FASTK_HIFI_DATABASE } from "./subworkflows/local/build_fastk_database/main"
include { BUILD_FASTK_DATABASE as BUILD_FASTK_HIC_DATABASE  } from "./subworkflows/local/build_fastk_database/main"
include { BUILD_MERYL_DATABASE as BUILD_MERYL_HIFI_DATABASE } from "./subworkflows/local/build_meryl_database/main"
include { BUILD_MERYL_DATABASE as BUILD_MERYL_HIC_DATABASE  } from "./subworkflows/local/build_meryl_database/main"
// Data conversion
include { CONVERT_FASTQ_CRAM                                } from "./subworkflows/local/convert_fastq_cram/main"
// Data inspection
include { INSPECT_DATA                                      } from "./subworkflows/local/inspect_data/main"
// Assembly
include { ASSEMBLE                                          } from "./subworkflows/local/assemble/main"
include { ASSEMBLE_ORGANELLES                               } from "./subworkflows/local/assemble_organelles/main"
// Decontamination
include { DECONTAMINATE                                     } from "./subworkflows/local/decontaminate/main"
// Purge duplicates
include { PURGE_DUPLICATES                                  } from "./subworkflows/local/purge_dups/main"
// Scaffold
include { SCAFFOLD                                          } from "./subworkflows/local/scaffold/main.nf"
// Curation
include { SCAFFOLD_CURATION                                 } from "./subworkflows/local/scaffold_curation/main.nf"
// Evaluate assemblies
include { COMPARE_ASSEMBLIES                                } from "./subworkflows/local/compare_assemblies/main"
include { EVALUATE_ASSEMBLY                                 } from "./subworkflows/local/evaluate_assembly/main"
// Align RNAseq
include { ALIGN_RNASEQ                                      } from "./subworkflows/local/align_rnaseq/main"
// Report
include { ASSEMBLY_REPORT                                   } from "./subworkflows/local/assembly_report/main"

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
    // ch_quarto_files  = Channel.empty()
    ch_versions      = Channel.empty()

    // Read in data
    PREPARE_INPUT (
        params.input,
        // params.ncbi.taxdb
    )
    ch_evaluate_assemblies = PREPARE_INPUT.out.assemblies

    // Create Cram and trim
    CONVERT_FASTQ_CRAM ( PREPARE_INPUT.out.hic, params.hic_type.startsWith('arima') )
    // Build necessary databases
    // TODO: Migrate back to Meryl. Genome inspection missing KATGC and PLOIDYPLOT for meryldb
    BUILD_FASTK_HIFI_DATABASE ( PREPARE_INPUT.out.hifi )
    BUILD_FASTK_HIC_DATABASE ( CONVERT_FASTQ_CRAM.out.fastq )
    BUILD_MERYL_HIFI_DATABASE ( PREPARE_INPUT.out.hifi )
    BUILD_MERYL_HIC_DATABASE ( CONVERT_FASTQ_CRAM.out.fastq )
    ch_versions = ch_versions.mix(
        CONVERT_FASTQ_CRAM.out.versions,
        BUILD_FASTK_HIFI_DATABASE.out.versions,
        BUILD_FASTK_HIC_DATABASE.out.versions,
        BUILD_MERYL_HIFI_DATABASE.out.versions,
        BUILD_MERYL_HIC_DATABASE.out.versions,
    )

    // Data inspection
    ch_hifi = PREPARE_INPUT.out.hifi
    if ( 'inspect' in workflow_steps ) {
        // QC Steps
        INSPECT_DATA(
            PREPARE_INPUT.out.hifi,
            PREPARE_INPUT.out.hic,
            BUILD_FASTK_HIFI_DATABASE.out.fastk_hist_ktab,
            BUILD_FASTK_HIC_DATABASE.out.fastk_hist_ktab
        )
        ch_hifi = INSPECT_DATA.out.hifi // with added kmer coverage
        ch_multiqc_files = ch_multiqc_files.mix( INSPECT_DATA.out.logs )
        ch_versions = ch_versions.mix( INSPECT_DATA.out.versions )
    }

    // Preprocess data
    if ( 'preprocess' in workflow_steps ) {
        // Adapter filtering etc
        // Subsampling
        // Host contamination filter
    }

    // Assemble
    if ( 'assemble' in workflow_steps ) {
        // Run assemblers
        ASSEMBLE ( PREPARE_INPUT.out.hifi_merged )
        ch_evaluate_assemblies = ch_evaluate_assemblies.mix( ASSEMBLE.out.raw_assemblies )
        ch_raw_assemblies = ASSEMBLE.out.raw_assemblies
        ch_multiqc_files = ch_multiqc_files.mix( ASSEMBLE.out.logs )
        ch_versions = ch_versions.mix( ASSEMBLE.out.versions )
    } else {
        ch_raw_assemblies = Channel.empty() // No assemblies from a previous stage
    }
    ch_raw_assemblies = ch_raw_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'raw' )
    ).dump(tag: 'Assemblies: Raw', pretty: true)

    // Organelle assembly
    if ( params.organelle_assembly_mode == 'reads' ) {
        // TODO: Add organelle assembly from reads
    } else if ( params.organelle_assembly_mode == 'contigs' ){
        ASSEMBLE_ORGANELLES ( ch_raw_assemblies )
        ch_versions = ch_versions.mix(ASSEMBLE_ORGANELLES.out.versions)
        // TODO: filter organelles from assemblies
    } // else params.organelle_assembly_mode == 'none'

    // Contamination screen
    ch_to_screen = setAssemblyStage (
        ch_raw_assemblies,
        'decontaminated' // Set assembly stage now for filenaming
    ).dump(tag: 'Assemblies: to screen', pretty: true)
    if ( 'screen' in workflow_steps ) {
        DECONTAMINATE( ch_to_screen )
        ch_evaluate_assemblies = ch_evaluate_assemblies.mix( DECONTAMINATE.out.assemblies )
        ch_cleaned_assemblies = DECONTAMINATE.out.assemblies
        ch_multiqc_files = ch_multiqc_files.mix( DECONTAMINATE.out.logs )
        ch_versions = ch_versions.mix(DECONTAMINATE.out.versions)
    } else {
        ch_cleaned_assemblies = ch_to_screen
    }
    ch_cleaned_assemblies = ch_cleaned_assemblies.mix (
        preassembledInput( PREPARE_INPUT.out.assemblies, 'decontaminated' )
    ).dump(tag: 'Assemblies: Cleaned', pretty: true)

    // Purge duplicates
    ch_to_purge = setAssemblyStage (
        ch_cleaned_assemblies,
        'purged' // Set assembly stage now for filenaming
    ).dump(tag: 'Assemblies: to purge', pretty: true)
    if ( 'purge' in workflow_steps ) {
        PURGE_DUPLICATES (
            ch_to_purge,
            ch_hifi
        )
        ch_evaluate_assemblies = ch_evaluate_assemblies.mix( PURGE_DUPLICATES.out.assemblies )
        ch_purged_assemblies = PURGE_DUPLICATES.out.assemblies
        ch_multiqc_files = ch_multiqc_files.mix( PURGE_DUPLICATES.out.logs )
        ch_versions = ch_versions.mix( PURGE_DUPLICATES.out.versions )
    } else {
        ch_purged_assemblies = ch_to_purge
    }
    ch_purged_assemblies = ch_purged_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'purged' )
    ).dump(tag: 'Assemblies: Purged', pretty: true)

    // Polish
    ch_to_polish = setAssemblyStage (
        ch_purged_assemblies,
        'polished' // Set assembly stage now for filenaming
    ).dump(tag: 'Assemblies: to polish', pretty: true)
    if ( 'polish' in workflow_steps ) {
        // Run polishers
        ch_polished_assemblies = ch_to_polish
    } else {
        ch_polished_assemblies = ch_to_polish
    }
    ch_polished_assemblies = ch_polished_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'polished' )
    ).dump(tag: 'Assemblies: Polished', pretty: true)

    // Scaffold
    ch_to_scaffold = setAssemblyStage (
        ch_polished_assemblies,
        'scaffolded' // Set assembly stage now for filenaming
    ).dump(tag: 'Assemblies: to scaffold', pretty: true)
    if ( 'scaffold' in workflow_steps ) {
        SCAFFOLD (
            ch_to_scaffold,
            CONVERT_FASTQ_CRAM.out.fastq
        )
        ch_evaluate_assemblies = ch_evaluate_assemblies.mix( SCAFFOLD.out.assemblies )
        ch_scaffolded_assemblies = SCAFFOLD.out.assemblies
        ch_multiqc_files = ch_multiqc_files.mix( SCAFFOLD.out.logs )
        ch_versions = ch_versions.mix( SCAFFOLD.out.versions )
    } else {
        ch_scaffolded_assemblies = ch_to_scaffold
    }
    ch_scaffolded_assemblies = ch_scaffolded_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'scaffolded' )
    ).dump(tag: 'Assemblies: Scaffolded', pretty: true)

    // Curate
    ch_to_curate = setAssemblyStage (
        ch_scaffolded_assemblies,
        'curated' // Set assembly stage now for filenaming
    ).dump(tag: 'Assemblies: to curate', pretty: true)
    if ( 'curate' in workflow_steps ) {
        SCAFFOLD_CURATION (
            ch_to_curate,
            CONVERT_FASTQ_CRAM.out.fastq,
            PREPARE_INPUT.out.hifi
        )
        ch_versions = ch_versions.mix( SCAFFOLD_CURATION.out.versions )
    }
    preassembledInput( PREPARE_INPUT.out.assemblies, 'curated' ).dump(tag: 'Assemblies: Curated')

    // Align RNAseq
    if( 'alignRNA' in workflow_steps ) {
        ALIGN_RNASEQ (
            PREPARE_INPUT.out.rnaseq,
            PREPARE_INPUT.out.assemblies // TODO: Select assembly stage
                .map { meta, assembly -> [ meta, assembly.pri_fasta ] }
        )
    }

    // Evaluate assemblies
    COMPARE_ASSEMBLIES ( ch_evaluate_assemblies )
    EVALUATE_ASSEMBLY (
        ch_evaluate_assemblies
            .toList().flatMap(), // Introduce bottleneck to delay evaluation until previous steps are performed
        BUILD_FASTK_HIFI_DATABASE.out.fastk_hist_ktab,
        BUILD_MERYL_HIFI_DATABASE.out.uniondb
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        COMPARE_ASSEMBLIES.out.logs,
        EVALUATE_ASSEMBLY.out.logs
    )
    ch_versions = ch_versions.mix(
        COMPARE_ASSEMBLIES.out.versions,
        EVALUATE_ASSEMBLY.out.versions
    )

    ASSEMBLY_REPORT(
        PREPARE_INPUT.out.sample_meta.map{ meta ->
            [
                meta,
                file(params.quarto_assembly_report, checkIfExists: true),
                files(params.quarto_assembly_report_aux_files, checkIfExists: true)
            ]
        },
        ch_multiqc_files,
        ch_versions,
        [ diagnostics: "debug" in workflow.profile.tokenize(",") ] +
            workflow_permitted_stages.collectEntries{ step -> [(step): step in params.steps.tokenize(",")] }
    )

    workflow.onComplete = {
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
}
