#! /usr/bin/env nextflow
include { FASTK_HISTEX } from "$projectDir/modules/nf-core/fastk/histex/main"
include { GENESCOPEFK  } from "$projectDir/modules/nf-core/genescopefk/main"

include { MERQURYFK_PLOIDYPLOT } from "$projectDir/modules/nf-core/merquryfk/ploidyplot/main"
include { MERQURYFK_KATGC      } from "$projectDir/modules/nf-core/merquryfk/katgc/main"

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
    FASTK_HISTEX ( fastk_hist_ktab.map { meta, hist, ktab -> [ meta, hist ] } )
    GENESCOPEFK ( FASTK_HISTEX.out.hist )

    // Print warning if genome size estimate is outside predicted range.
    GENESCOPEFK.out.summary
        .view { meta, summary ->
            def genome_size_estimates = summary.readLines().find { it.startsWith("Genome Haploid Length") } =~ /[0-9,]+/
            if ( genome_size_estimates.size() == 1 ) {
                if ( meta.sample.genome_size < 0.9 * genome_size_estimates[0] ||
                    meta.sample.genome_size > 1.1 * genome_size_estimates[0] ) {
                        log.warn "GeneScopeFK genome size estimate differs from GOAT estimate"
                }
            } else {
                // Min and Max estimate
                if ( meta.sample.genome_size < genome_size_estimates[0] ||
                    meta.sample.genome_size > genome_size_estimates[1] ) {
                        log.warn "GeneScopeFK genome size estimate differs from GOAT estimate"
                }
            }
        }

    // Generate Smudgeplot
    MERQURYFK_PLOIDYPLOT ( fastk_hist_ktab )

    // Generage GC plot
    MERQURYFK_KATGC ( fastk_hist_ktab )

    versions_ch = FASTK_HISTEX.out.versions.first().mix(
            GENESCOPEFK.out.versions.first(),
            MERQURYFK_PLOIDYPLOT.out.versions.first(),
            MERQURYFK_KATGC.out.versions.first()
        )

    emit:
    kmer_cov = GENESCOPEFK.out.kmer_cov
    versions = versions_ch
}
