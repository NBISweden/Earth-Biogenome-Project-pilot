#! /usr/bin/env nextflow
include { FASTK_HISTEX } from "../../../../modules/nf-core/fastk/histex/main"
include { GENESCOPEFK  } from "../../../../modules/nf-core/genescopefk/main"

include { MERQURYFK_PLOIDYPLOT } from "../../../../modules/nf-core/merquryfk/ploidyplot/main"
include { MERQURYFK_KATGC      } from "../../../../modules/nf-core/merquryfk/katgc/main"

workflow GENOME_PROPERTIES {

    take:
    fastk_hist_ktab   // [ meta, fastk_hist, fastk_ktab ]

    /* Genome properties workflow:
        - Estimate genome depth of coverage from reads
        - Generate k-mer histogram
        - Smudgeplot
    */
    main:
    // Generate GenomeScope Profile
    FASTK_HISTEX ( fastk_hist_ktab.map { meta, hist, _ktab -> [ meta, hist ] } )
    GENESCOPEFK ( FASTK_HISTEX.out.hist )

    // Generate Smudgeplot
    MERQURYFK_PLOIDYPLOT ( fastk_hist_ktab )

    // Generage GC plot
    MERQURYFK_KATGC ( fastk_hist_ktab )

    GENESCOPEFK.out.linear_plot.mix (
        GENESCOPEFK.out.log_plot,
        GENESCOPEFK.out.transformed_linear_plot,
        GENESCOPEFK.out.transformed_log_plot,
        MERQURYFK_PLOIDYPLOT.out.stacked_ploidy_plot_png,
        MERQURYFK_KATGC.out.stacked_gc_plot_png
    )
    .map { _meta, img -> img }
    .set { logs }

    FASTK_HISTEX.out.versions.first().mix(
        GENESCOPEFK.out.versions.first(),
        MERQURYFK_PLOIDYPLOT.out.versions.first(),
        MERQURYFK_KATGC.out.versions.first()
    ).set { versions }

    emit:
    kmer_cov = GENESCOPEFK.out.kmer_cov
    logs
    versions
}
