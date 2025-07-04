include { FASTP           } from '../../../modules/nf-core/fastp/main'
include { SAMTOOLS_IMPORT } from '../../../modules/nf-core/samtools/import/main'
include { SAMTOOLS_INDEX  } from '../../../modules/nf-core/samtools/index/main'

workflow CONVERT_FASTQ_CRAM {
    take:
    ch_fastq        // Channel: [ meta, fastq ]
    bool_trim_fastq // Boolean: Trim with Fastp

    main:
    // Filter
    FASTP(
        ch_fastq.filter { bool_trim_fastq },
        [],
        false,
        false,
        false,
    )

    // Convert to cram
    SAMTOOLS_IMPORT(ch_fastq.filter { !bool_trim_fastq }.mix(FASTP.out.reads))

    // Index Cram
    SAMTOOLS_INDEX(SAMTOOLS_IMPORT.out.cram)

    // Versions
    FASTP.out.versions.first().mix(
        SAMTOOLS_IMPORT.out.versions.first(),
        SAMTOOLS_INDEX.out.versions.first(),
    ).set { versions }

    emit:
    fastq = bool_trim_fastq ? FASTP.out.reads : ch_fastq
    cram  = SAMTOOLS_IMPORT.out.cram
    crai  = SAMTOOLS_INDEX.out.crai
    versions
}
