#! /usr/bin/env nextflow

include { PREPARE_INPUT } from "$projectDir/subworkflows/local/prepare_input/main"

include { BUILD_DATABASES as BUILD_HIFI_DATABASES } from "$projectDir/subworkflows/local/build_databases/main"
include { BUILD_DATABASES as BUILD_HIC_DATABASES  } from "$projectDir/subworkflows/local/build_databases/main"

include { GENOME_PROPERTIES } from "$projectDir/subworkflows/local/genome_properties/main"
include { COMPARE_LIBRARIES } from "$projectDir/subworkflows/local/compare_libraries/main"
include { SCREEN_READS      } from "$projectDir/subworkflows/local/screen_read_contamination/main"

include { HIFIASM  } from "$projectDir/modules/nf-core/hifiasm/main"
include { GFASTATS } from "$projectDir/modules/nf-core/gfastats/main"

include { PURGE_DUPLICATES } from "$projectDir/subworkflows/local/purge_dups/main"

include { COMPARE_ASSEMBLIES } from "$projectDir/subworkflows/local/compare_assemblies/main"
include { EVALUATE_ASSEMBLY  } from "$projectDir/subworkflows/local/evaluate_assembly/main"
include { ALIGN_RNASEQ       } from "$projectDir/subworkflows/local/align_rnaseq/main"

workflow {

    // Define constants
    def workflow_permitted_stages = [
        'inspect',      // 01 - Read inspection
        'preprocess',   // 02 - Read preprocessing
        'assemble',     // 03 - Assembly
        'purge',        // 04 - Duplicate purging
        'polish',       // 05 - Error polishing
        'screen',       // 06 - Contamination screening
        'scaffold',     // 07 - Scaffolding
        'curate'        // 08 - Rapid curation
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
    
    ch_assemblies = PREPARE_INPUT.out.assemblies
    // Assemble
    if ( 'assemble' in workflow_steps ) {
        // Run assemblers
        // Run in basic mode atm
        // Need to include a build ID here. 
        HIFIASM(
            PREPARE_INPUT.out.hifi,
            [], // paternal k-mers
            [], // maternal k-mers
            [], // Hi-C r1
            []  // Hi-C r2
        )
        GFASTATS( 
            HIFIASM.out.primary_contigs.mix( HIFIASM.out.alternate_contigs ),
            "fasta", // output format
            "",      // genome size
            "",      // target
            [],      // AGP file
            [],      // include bed
            [],      // exclude bed
            []       // SAK instructions
        )
        ch_hifiasm_out = GFASTATS.out.assembly.groupTuple()
            .join( HIFIASM.out.primary_contigs )
            .join( HIFIASM.out.alternate_contigs )
            .map { meta, fasta, pri_gfa, alt_gfa -> 
                [ meta, 
                    [ 
                        id: 'hifiasm', 
                        pri_fasta: fasta[0],  // These may not be grouped in primary, alternate order
                        alt_fasta: fasta[1],
                        pri_gfa: pri_gfa,
                        alt_gfa: alt_gfa
                    ] 
                ]
            }
        ch_assemblies = ch_assemblies.mix( ch_hifiasm_out )
        // Find mitochondria
    }

    // Purge duplicates
    if ( 'purge' in workflow_steps ) {
        ch_topurge = PREPARE_INPUT.out.hifi
                .map { meta, reads -> [ meta.findAll { ! (it.key in [ 'single_end' ]) }, reads ] }
                .combine( ch_assemblies, by:0 )
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
    }

    // Polish
    if ( 'polish' in workflow_steps ) {
        // Run assemblers
    }

    // Contamination screen
    if ( 'screen' in workflow_steps ) {
        // Kraken2
        // Blobtoolkit
        // FCS-Genome
    }

    // Scaffold
    if ( 'scaffold' in workflow_steps ) {
        // Run scaffolder
    }

    // Curate
    if ( 'curate' in workflow_steps ) {
        // Run assemblers
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

    }

    // Align RNAseq
    if( 'alignRNA' in workflow_steps ) {
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
