#! /usr/bin/env nextflow

// Functions
include { assembliesFromStage as preassembledInput          } from "./modules/local/functions"
include { setAssemblyStage                                  } from "./modules/local/functions"
include { assemblyRecordsToTuples                           } from "./modules/local/functions"
include { assemblyTuplesToRecords                           } from "./modules/local/functions"
// Data import
include { PREPARE_INPUT                                     } from "./subworkflows/local/prepare_input/main"
// Database preparation
include { BUILD_FASTK_DATABASE as BUILD_FASTK_HIFI_DATABASE } from "./subworkflows/local/build_fastk_database/main"
include { BUILD_FASTK_DATABASE as BUILD_FASTK_HIC_DATABASE  } from "./subworkflows/local/build_fastk_database/main"
include { BUILD_MERYL_DATABASE as BUILD_MERYL_HIFI_DATABASE } from "./subworkflows/local/build_meryl_database/main"
include { BUILD_MERYL_DATABASE as BUILD_MERYL_HIC_DATABASE  } from "./subworkflows/local/build_meryl_database/main"
// Data conversion
include { CONVERT_FASTQ_CRAM                                } from "./subworkflows/local/02_convert_fastq_cram/main"
// Data inspection
include { INSPECT_DATA                                      } from "./subworkflows/local/01_inspect_data/main"
// Assembly
include { ASSEMBLE                                          } from "./subworkflows/local/03_assemble/main"
// Decontamination
include { DECONTAMINATE                                     } from "./subworkflows/local/04_decontaminate/main"
// Purge duplicates
include { PURGE_DUPLICATES                                  } from "./subworkflows/local/05_purge_dups/main"
// Polish
include { DVPOLISH                                          } from "./subworkflows/local/06_polish/main"
// Scaffold
include { SCAFFOLD                                          } from "./subworkflows/local/07_scaffold/main"
// Curation
include { SCAFFOLD_CURATION                                 } from "./subworkflows/local/08_scaffold_curation/main"
// Evaluate assemblies
include { COMPARE_ASSEMBLIES                                } from "./subworkflows/local/compare_assemblies/main"
include { EVALUATE_ASSEMBLY                                 } from "./subworkflows/local/evaluate_assembly/main"
// Align RNAseq
include { ALIGN_RNASEQ                                      } from "./subworkflows/local/09_align_rnaseq/main"
// Tool versions report
include { REPORT_VERSIONS                                   } from "./subworkflows/local/version_report/main"
// Assembly report
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
    if ( ! workflow_steps.every { step -> step in workflow_permitted_stages } ) {
        error "Unrecognised workflow step in $params.steps ( $workflow_permitted_stages )"
    }

    // The primary workflow for the Earth Biogenome Project Pilot
    log.info("""
    Running NBIS Earth Biogenome Project Assembly workflow.
    """)

    // Setup sink channels
    ch_multiqc_files = channel.value( file(params.multiqc_assembly_report_config, checkIfExists: true) )

    // Read in data
    PREPARE_INPUT (
        params.input
    )
    ch_evaluate_assemblies = PREPARE_INPUT.out.assemblies

    // Typed inputs: [ meta, reads ] tuples to records
    ch_hifi_reads = PREPARE_INPUT.out.hifi.map { meta, reads -> record(id: meta.id, meta: record(meta), reads: [ reads ].flatten()) }
    ch_hic_reads  = PREPARE_INPUT.out.hic.map { meta, reads -> record(id: meta.id, meta: record(meta), reads: [ reads ].flatten()) }
    def fastk_kmer_size = params.fastk?.kmer_size ?: 31
    def meryl_kmer_size = params.meryl.kmer_size
    // Create Cram and trim
    CONVERT_FASTQ_CRAM ( ch_hic_reads, params.hic_type.startsWith('arima'), params.hic_type )
    ch_hic_fastq = CONVERT_FASTQ_CRAM.out.fastq
    // Build necessary databases
    // TODO: Migrate back to Meryl. Genome inspection missing KATGC and SMUDGEPLOT for meryldb
    ch_fastk_hifi = BUILD_FASTK_HIFI_DATABASE ( ch_hifi_reads, 'hifi', fastk_kmer_size )
    ch_fastk_hic  = BUILD_FASTK_HIC_DATABASE ( ch_hic_fastq, 'hic', fastk_kmer_size )
    ch_meryl_hifi = BUILD_MERYL_HIFI_DATABASE ( ch_hifi_reads, 'hifi', meryl_kmer_size )
    BUILD_MERYL_HIC_DATABASE ( ch_hic_fastq, 'hic', meryl_kmer_size )

    // Data inspection
    ch_hifi_kmer_cov = PREPARE_INPUT.out.hifi_merged.map { meta, _reads -> tuple(meta, []) }
    if ( 'inspect' in workflow_steps ) {
        // QC Steps
        ch_inspect = INSPECT_DATA(
            ch_hifi_reads,
            ch_hic_reads,
            ch_fastk_hifi,
            ch_fastk_hic,
            params.genescopefk?.kmer_size ?: 31
        )
        ch_hifi_kmer_cov = ch_inspect.hifi_kmer_cov.map { r -> tuple([:] + r.meta, r.kmer_cov) }
        ch_multiqc_files = ch_multiqc_files.mix( ch_inspect.logs )
    }

    // Preprocess data
    if ( 'preprocess' in workflow_steps ) {
        // Adapter filtering etc
        // Subsampling
        // Host contamination filter
    }

    // Assemble
    ch_hifi_merged_reads = PREPARE_INPUT.out.hifi_merged.map { meta, reads -> record(id: meta.id, meta: record(meta), reads: [ reads ].flatten()) }
    if ( 'assemble' in workflow_steps ) {
        // Run assemblers
        ch_assembled = ASSEMBLE (
                ch_hifi_merged_reads,
                ch_hic_reads,
                params.oatkdb ? file(params.oatkdb, checkIfExists: true, type: 'dir') : null,
                params.nuclear_assembly_mode,
                params.organelle_assembly_mode,
                params.hifiasm ?: [:],
                params.hifiasm_with_hic,
                params.use_phased
        )
        ch_raw_assemblies = assemblyRecordsToTuples( ch_assembled.raw_assemblies )
        ch_evaluate_assemblies = ch_evaluate_assemblies.mix( ch_raw_assemblies )
        ch_multiqc_files = ch_multiqc_files.mix( ch_assembled.logs )
    } else {
        ch_raw_assemblies = channel.empty() // No assemblies from a previous stage
    }
    ch_raw_assemblies = ch_raw_assemblies.mix(
        preassembledInput( PREPARE_INPUT.out.assemblies, 'raw' )
    ).dump(tag: 'Assemblies: Raw', pretty: true)

    // Contamination screen
    ch_to_screen = setAssemblyStage (
        ch_raw_assemblies,
        'decontaminated' // Set assembly stage now for filenaming
    ).dump(tag: 'Assemblies: to screen', pretty: true)
    if ( 'screen' in workflow_steps ) {
        ch_screened = DECONTAMINATE(
            assemblyTuplesToRecords( ch_to_screen ),
            params.fcs.database ? file(params.fcs.database, checkIfExists: true, type: 'dir') : null,
            params.fcs.manifest,
            params.fcs.ramdisk_path
        )
        ch_cleaned_assemblies = assemblyRecordsToTuples( ch_screened.assemblies )
        ch_evaluate_assemblies = ch_evaluate_assemblies.mix( ch_cleaned_assemblies )
        ch_multiqc_files = ch_multiqc_files.mix( ch_screened.logs )
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
        ch_purged = PURGE_DUPLICATES (
            assemblyTuplesToRecords( ch_to_purge ),
            ch_hifi_merged_reads,
            ch_hifi_kmer_cov.map { meta, kmercov -> record(id: meta.id, kmer_cov: kmercov ? "${kmercov}" : null) }
        )
        ch_purged_assemblies = assemblyRecordsToTuples( ch_purged.assemblies )
        ch_evaluate_assemblies = ch_evaluate_assemblies.mix( ch_purged_assemblies )
        ch_multiqc_files = ch_multiqc_files.mix( ch_purged.logs )
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
        ch_polished = DVPOLISH (
            assemblyTuplesToRecords( ch_to_polish ),
            ch_hifi_merged_reads,
            ch_meryl_hifi.uniondb,
            MemoryUnit.of(params.chunk_size).toBytes(),
            params.use_phased
        )
        ch_polished_assemblies = assemblyRecordsToTuples( ch_polished.assemblies )
        ch_evaluate_assemblies = ch_evaluate_assemblies.mix( ch_polished_assemblies )
        ch_multiqc_files = ch_multiqc_files.mix( ch_polished.logs )
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
        ch_scaffolded = SCAFFOLD (
            assemblyTuplesToRecords( ch_to_scaffold ),
            ch_hic_fastq,
            params.use_phased,
            params.hic_type,
            params.hic_bam_sorting
        )
        ch_scaffolded_assemblies = assemblyRecordsToTuples( ch_scaffolded.assemblies )
        ch_evaluate_assemblies = ch_evaluate_assemblies.mix( ch_scaffolded_assemblies )
        ch_multiqc_files = ch_multiqc_files.mix( ch_scaffolded.logs )
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
            assemblyTuplesToRecords( ch_to_curate ),
            ch_hic_fastq,
            ch_hifi_reads,
            params.hic_map_sort_by,
            params.hic_map_qv,
            params.cooler_bin_size,
            params.hifi_coverage_cap,
            params.telomer_motif
        )
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
    ch_to_evaluate = assemblyTuplesToRecords( ch_evaluate_assemblies )
    ch_compared = COMPARE_ASSEMBLIES (
        ch_to_evaluate,
        params.reference ? file( params.reference, checkIfExists: true ) : null
    )
    ch_evaluated = EVALUATE_ASSEMBLY (
        assemblyTuplesToRecords(
            ch_evaluate_assemblies
                .toList().flatMap() // Introduce bottleneck to delay evaluation until previous steps are performed
        ),
        ch_fastk_hifi,
        ch_meryl_hifi.uniondb,
        params.use_phased,
        params.busco.max_lineages,
        params.busco.lineages_db_path ? file( params.busco.lineages_db_path, checkIfExists: true ) : null
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_compared,
        ch_evaluated
    )

    // Version reporting
    REPORT_VERSIONS(
        channel.topic("versions")
    )

    // Assembly report
    ASSEMBLY_REPORT(
        PREPARE_INPUT.out.sample_meta.map{ meta ->
            [
                meta,
                file(params.quarto_assembly_report, checkIfExists: true),
                files(params.quarto_assembly_report_aux_files, checkIfExists: true)
            ]
        },
        ch_multiqc_files,
        REPORT_VERSIONS.out.multiqc_versions_report_ch,
        [
            diagnostics: "debug" in workflow.profile.tokenize(","),
            nuclear_assembly_mode: params.nuclear_assembly_mode,
            organelle_assembly_mode: params.organelle_assembly_mode in ['contigs', 'reads']
        ] + workflow_permitted_stages.collectEntries{ step ->
            [(step): step in params.steps.tokenize(",")]
        }
    )

    // Completion message
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
