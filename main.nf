#! /usr/bin/env nextflow

include { PREPARE_INPUT } from "$projectDir/subworkflows/local/prepare_input/main"

include { BUILD_DATABASES as BUILD_HIFI_DATABASES } from "$projectDir/subworkflows/local/build_databases/main"
include { BUILD_DATABASES as BUILD_HIC_DATABASES  } from "$projectDir/subworkflows/local/build_databases/main"

include { GENOME_PROPERTIES } from "$projectDir/subworkflows/local/genome_properties/main"
include { COMPARE_LIBRARIES } from "$projectDir/subworkflows/local/compare_libraries/main"
include { SCREEN_READS      } from "$projectDir/subworkflows/local/screen_read_contamination/main"

include { ASSEMBLE_HIFI } from "$projectDir/subworkflows/local/assemble_hifi/main"

include { FCSGX_FETCHDB } from "$projectDir/modules/local/fcs/fetchdb"
include { FCS_FCSGX     } from "$projectDir/modules/nf-core/fcs/fcsgx/main"

include { PURGE_DUPLICATES } from "$projectDir/subworkflows/local/purge_dups/main"

include { MITOHIFI_FINDMITOREFERENCE } from "$projectDir/modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI          } from "$projectDir/modules/nf-core/mitohifi/mitohifi/main"

include { COMPARE_ASSEMBLIES } from "$projectDir/subworkflows/local/compare_assemblies/main"
include { EVALUATE_ASSEMBLY  } from "$projectDir/subworkflows/local/evaluate_assembly/main"
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
    PREPARE_INPUT ( params.input )

    // Build necessary databases
    // if ( ['inspect','preprocess','assemble','purge','polish','screen','scaffold','curate'].any{ it in workflow_steps}) {
    BUILD_HIFI_DATABASES ( PREPARE_INPUT.out.hifi )
    BUILD_HIC_DATABASES ( PREPARE_INPUT.out.hic )
    // }

    // Data inspection
    if ( 'inspect' in workflow_steps ) {
        // QC Steps
        GENOME_PROPERTIES (
            BUILD_HIFI_DATABASES.out.fastk_histogram.join( BUILD_HIFI_DATABASES.out.fastk_ktab )
        )
        COMPARE_LIBRARIES (
            BUILD_HIFI_DATABASES.out.fastk_histogram.join( BUILD_HIFI_DATABASES.out.fastk_ktab ).join(
            BUILD_HIC_DATABASES.out.fastk_histogram.join( BUILD_HIC_DATABASES.out.fastk_ktab ) )
        )
        SCREEN_READS (
            PREPARE_INPUT.out.hifi,
            // TODO:: Allow custom database ala nf-core/genomeassembler.
            file( params.mash.screen_db, checkIfExists: true )
        )
    }

    // Preprocess data
    if ( 'preprocess' in workflow_steps ) {
        // Adapter filtering etc
    }

    ch_assemblies = PREPARE_INPUT.out.assemblies
    // Assemble
    if ( 'assemble' in workflow_steps ) {
        // Run assemblers
        ASSEMBLE_HIFI( PREPARE_INPUT.out.hifi )
        ch_assemblies = ch_assemblies.filter { meta, assembly -> meta.assembly.stage in ['raw'] }
            .mix( ASSEMBLE_HIFI.out.assemblies )

        // Find mitochondria
        // Need to check options to mitohifi modules.
        MITOHIFI_FINDMITOREFERENCE( ch_assemblies.map { meta, assemblies -> [ meta, meta.sample.name ] }.unique() )
        mitohifi_ch = ch_assemblies
            .combine( MITOHIFI_FINDMITOREFERENCE.out.fasta, by: 0 )
            .combine( MITOHIFI_FINDMITOREFERENCE.out.gb, by: 0 )
            .multiMap { meta, assembly, mitofa, mitogb ->
                input: [ meta, assembly.pri_fasta ]
                reference: mitofa
                genbank: mitogb
            }
        MITOHIFI_MITOHIFI(
            mitohifi_ch.input,
            mitohifi_ch.reference,
            mitohifi_ch.genbank,
            'c',
            params.mitohifi.code
        )

        // Assess assemblies
        COMPARE_ASSEMBLIES (
            ch_assemblies,
            params.reference ? file( params.reference, checkIfExists: true ) : []
        )

        EVALUATE_ASSEMBLY (
            ch_assemblies,
            PREPARE_INPUT.out.hifi,
            BUILD_HIFI_DATABASES.out.fastk_histogram.join( BUILD_HIFI_DATABASES.out.fastk_ktab ),
            params.reference ? file( params.reference, checkIfExists: true ) : [],
            params.busco.lineages_db_path ? file( params.busco.lineages_db_path, checkIfExists: true ) : []
        )
    }

    // Contamination screen
    if ( 'screen' in workflow_steps ) {
        FCSGX_FETCHDB ( params.fcs.database ? Channel.empty() : Channel.fromPath( params.fcs.manifest, checkIfExists: true ) )
        ch_fcs_database = params.fcs.database ? Channel.fromPath( params.fcs.database, checkIfExists: true, type: 'dir' ) : FCSGX_FETCHDB.out.database
        // Do we need a separate stage here?
        ch_to_screen = ch_assemblies.filter { meta, assembly -> meta.assembly.stage in ['raw'] }
            .flatMap { meta, assembly ->
                assembly.alt_fasta ? [ [ meta, assembly.pri_fasta ], [ meta, assembly.alt_fasta ] ] : [ [ meta, assembly.pri_fasta ] ]
            }
        FCS_FCSGX( ch_to_screen, ch_fcs_database.collect() )
    }

    // Purge duplicates
    if ( 'purge' in workflow_steps ) {
        ch_topurge = PREPARE_INPUT.out.hifi.combine( ch_assemblies.filter { meta, assembly -> meta.assembly.stage in ['raw','decontaminated'] } , by:0 )
        if ( 'inspect' in workflow_steps ) {
            // Add kmer coverage from GenomeScope model
            ch_topurge.map { meta, reads, assemblies -> [ meta.subMap(['id','sample']), meta, reads, assemblies ] }
                .combine( GENOME_PROPERTIES.out.kmer_cov.map{ meta, cov -> [ meta.subMap(['id','sample']), cov ] } , by: 0 )
                .map { key, meta, reads, assemblies, kmer_cov -> [ meta + [ kmercov: kmer_cov ], reads, assemblies ] }
                .set { ch_topurge }
        }
        PURGE_DUPLICATES ( ch_topurge.dump( tag: 'Purge duplicates: input' ) )
    }

    // Polish
    if ( 'polish' in workflow_steps ) {
        // Run polishers
    }

    // Scaffold
    if ( 'scaffold' in workflow_steps ) {
        // Run scaffolder
    }

    // Curate
    if ( 'curate' in workflow_steps ) {
        // Run assemblers
    }

    // Align RNAseq
    if( 'alignRNA' in workflow_steps ) {
        ALIGN_RNASEQ (
            PREPARE_INPUT.out.rnaseq,
            PREPARE_INPUT.out.assemblies
                .map { meta, assembly -> [ meta, assembly.pri_fasta ] }
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

// Custom method for Map class called merge
// returns a new Map with entries merged.
Map.metaClass.merge = { Map rhs ->
    def lhs = delegate
    lhs.collectEntries { k, v -> rhs[k] instanceof Map ? [ (k): v.merge( rhs[k] ) ] : [ (k): v ] } + rhs.subMap(rhs.keySet()-lhs.keySet())
}
