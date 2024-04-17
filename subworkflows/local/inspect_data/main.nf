include { GENOME_PROPERTIES } from "$projectDir/subworkflows/local/genome_properties/main"
include { COMPARE_LIBRARIES } from "$projectDir/subworkflows/local/compare_libraries/main"
include { SCREEN_READS      } from "$projectDir/subworkflows/local/screen_read_contamination/main"

workflow INSPECT_DATA {
    take:
    hifi_reads     // [ meta, hifi ]
    hifi_histogram // [ meta, hifi_hist, ktab ]
    hic_histogram  // [ meta, hic_hist, ktab ]

    main:
    ch_versions = Channel.empty()
    // QC Steps
    GENOME_PROPERTIES ( hifi_histogram )
    COMPARE_LIBRARIES ( hifi_histogram.join( hic_histogram ) )
    SCREEN_READS (
        PREPARE_INPUT.out.hifi,
        // TODO:: Allow custom database ala nf-core/genomeassembler.
        file( params.mash.screen_db, checkIfExists: true )
    )

    emit:
    versions = ch_versions
}