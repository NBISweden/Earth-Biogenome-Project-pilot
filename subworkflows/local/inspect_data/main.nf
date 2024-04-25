include { combineByMetaKeys } from "$projectDir/modules/local/functions"
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
        hifi_reads,
        // TODO:: Allow custom database ala nf-core/genomeassembler.
        file( params.mash.screen_db, checkIfExists: true )
    )

    ch_hifi_with_kmer_cov = combineByMetaKeys(
        hifi_reads,
        GENOME_PROPERTIES.out.kmer_cov,
        keySet: ['id','sample'],
        meta: 'lhs'
    )
    .map { meta, reads, kmer_cov -> [ meta + [ kmercov: kmer_cov ], reads ] }

    emit:
    hifi     = ch_hifi_with_kmer_cov
    versions = ch_versions
}