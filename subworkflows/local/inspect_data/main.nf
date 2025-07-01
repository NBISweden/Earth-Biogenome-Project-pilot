include { combineByMetaKeys } from "$projectDir/modules/local/functions"
include { GENOME_PROPERTIES } from "$projectDir/subworkflows/local/genome_properties/main"
include { COMPARE_LIBRARIES } from "$projectDir/subworkflows/local/compare_libraries/main"

workflow INSPECT_DATA {
    take:
    hifi_reads     // [ meta, hifi ]
    hifi_histogram // [ meta, hifi_hist, ktab ]
    hic_histogram  // [ meta, hic_hist, ktab ]

    main:
    // QC Steps
    GENOME_PROPERTIES ( hifi_histogram )
    COMPARE_LIBRARIES ( hifi_histogram.join( hic_histogram ) )
    ch_versions =  GENOME_PROPERTIES.out.versions.mix(
        COMPARE_LIBRARIES.out.versions)

    ch_hifi_with_kmer_cov = combineByMetaKeys(
        hifi_reads,
        GENOME_PROPERTIES.out.kmer_cov,
        keySet: ['id','sample'],
        meta: 'lhs'
    )
    .map { meta, reads, kmer_cov -> [ meta + [ kmercov: kmer_cov ], reads ] }

    GENOME_PROPERTIES.out.logs
        .mix( COMPARE_LIBRARIES.out.logs )
        .set { logs }

    emit:
    hifi = ch_hifi_with_kmer_cov
    logs
    versions = ch_versions
}
