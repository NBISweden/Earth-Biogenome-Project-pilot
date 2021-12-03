#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KMC_HIST            } from "../../modules/kmc/kmc_hist/kmc_hist"
include { KMC_DUMP            } from "../../modules/kmc/kmc_dump/kmc_dump"
include { GENOMESCOPE         } from "../../modules/genomescope/genomescope"
include { SMUDGEPLOT_CUTOFF   } from "../../modules/smudgeplot/smudgeplot_cutoff/smudgeplot_cutoff"
include { SMUDGEPLOT_HETKMERS } from "../../modules/smudgeplot/smudgeplot_hetkmers/smudgeplot_hetkmers"
include { SMUDGEPLOT_PLOT     } from "../../modules/smudgeplot/smudgeplot_plot/smudgeplot_plot"

workflow {
    GENOME_PROPERTIES ( Channel.fromPath( params.reads ) )
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
    // Generate GenomeScope Profile
    KMC_HIST ( reads_ch )
    GENOMESCOPE ( KMC_HIST.out.histogram )

    // Generate Smudgeplot
    SMUDGEPLOT_CUTOFF ( KMC_HIST.out.count_db )
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

}
