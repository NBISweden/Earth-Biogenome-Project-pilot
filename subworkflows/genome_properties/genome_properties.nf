#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FASTQ      } from "../../modules/nf-core/modules/samtools/fastq/main"

include { KMC_HIST            } from "../../modules/local/kmc/kmc_hist/kmc_hist"
include { KMC_DUMP            } from "../../modules/local/kmc/kmc_dump/kmc_dump"
include { GENOMESCOPE         } from "../../modules/local/genomescope/genomescope"
include { SMUDGEPLOT_CUTOFF   } from "../../modules/local/smudgeplot/smudgeplot_cutoff/smudgeplot_cutoff"
include { SMUDGEPLOT_HETKMERS } from "../../modules/local/smudgeplot/smudgeplot_hetkmers/smudgeplot_hetkmers"
include { SMUDGEPLOT_PLOT     } from "../../modules/local/smudgeplot/smudgeplot_plot/smudgeplot_plot"

include { PREPARE_INPUT       } from "../../subworkflows/prepare_input/prepare_input"

workflow {
    PREPARE_INPUT ( params.input )
    GENOME_PROPERTIES ( PREPARE_INPUT.out.hifi )
}

workflow GENOME_PROPERTIES {

    take:
    reads_ch  // input type: [ [ id: 'sample_name' ], [ file('path/to/reads') ] ]

    /* Genome properties workflow:
        - Estimate genome depth of coverage from reads
        - Generate k-mer histogram
        - Smudgeplot
    */
    main:
    reads_ch.transpose()   // Transform to [ [ id: 'sample_name'], file('/path/to/read')  ]
        .branch { meta, filename ->
            bam_ch: filename.toString().endsWith(".bam")
            fastq_ch: true // assume everything else is fastq
        }.set { input }
    // Convert BAMS to FASTQ
    SAMTOOLS_FASTQ( input.bam_ch )

    // Generate GenomeScope Profile
    KMC_HIST ( input.fastq_ch
        .mix( SAMTOOLS_FASTQ.out.fastq )
        .groupTuple()
    )
    GENOMESCOPE ( KMC_HIST.out.histogram )

    // Generate Smudgeplot
    SMUDGEPLOT_CUTOFF ( KMC_HIST.out.histogram )
    KMC_DUMP ( KMC_HIST.out.count_db.join( SMUDGEPLOT_CUTOFF.out.bounds ) )
    SMUDGEPLOT_HETKMERS ( KMC_DUMP.out.histogram )
    SMUDGEPLOT_PLOT ( SMUDGEPLOT_HETKMERS.out.coverage_tsv )

    versions_ch = KMC_HIST.out.versions.first().mix(
        KMC_DUMP.out.versions.first(),
        GENOMESCOPE.out.versions.first(),
        SMUDGEPLOT_CUTOFF.out.versions.first(),
        SMUDGEPLOT_HETKMERS.out.versions.first(),
        SMUDGEPLOT_PLOT.out.versions.first()
    )

    emit:
    versions = versions_ch

}
