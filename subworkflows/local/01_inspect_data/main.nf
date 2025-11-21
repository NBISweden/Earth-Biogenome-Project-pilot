include { combineByMetaKeys                 } from "../../../modules/local/functions"
include { SEQKIT_STATS as SEQKIT_HIFI_STATS } from "../../../modules/nf-core/seqkit/stats/main"
include { SEQKIT_STATS as SEQKIT_HIC_STATS  } from "../../../modules/nf-core/seqkit/stats/main"
include { FASTQC                            } from "../../../modules/nf-core/fastqc/main"
include { GENOME_PROPERTIES                 } from "./genome_properties"
include { COMPARE_LIBRARIES                 } from "./compare_libraries"

workflow INSPECT_DATA {
    take:
    hifi_reads              // [ meta, hifi ]
    hic_reads               // [ meta, hic ]
    hifi_histogram          // [ meta, hifi_hist, ktab ]
    hic_histogram           // [ meta, hic_hist, ktab ]
    ch_smudgeplot_threshold // Int.

    main:
    // Quantify data
    SEQKIT_HIFI_STATS( hifi_reads )
    SEQKIT_HIC_STATS( hic_reads )
    FASTQC( hic_reads )

    // Generate K-mer histogram
    GENOME_PROPERTIES ( hifi_histogram, ch_smudgeplot_threshold )
    COMPARE_LIBRARIES ( hifi_histogram.join( hic_histogram ) )
    ch_versions =  GENOME_PROPERTIES.out.versions.mix(
        COMPARE_LIBRARIES.out.versions,
        SEQKIT_HIFI_STATS.out.versions.first(),
        SEQKIT_HIC_STATS.out.versions.first(),
        FASTQC.out.versions.first()
    )

    ch_hifi_with_kmer_cov = combineByMetaKeys(
        hifi_reads,
        GENOME_PROPERTIES.out.kmer_cov,
        keySet: ['id','sample'],
        meta: 'lhs'
    )
    .map { meta, reads, kmer_cov -> [ meta + [ kmercov: kmer_cov ], reads ] }

    GENOME_PROPERTIES.out.logs
        .mix(
            COMPARE_LIBRARIES.out.logs,
            SEQKIT_HIFI_STATS.out.stats.map { _meta, stats -> stats },
            SEQKIT_HIC_STATS.out.stats.map { _meta, stats -> stats},
            FASTQC.out.zip.map { _meta, zip -> zip }
        )
        .set { logs }

    emit:
    hifi = ch_hifi_with_kmer_cov
    logs
    versions = ch_versions
}
