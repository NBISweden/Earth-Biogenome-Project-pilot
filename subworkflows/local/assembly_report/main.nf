include { TOL_SEARCH                 } from "$projectDir/modules/local/tol/search"
include { REPORT_DTOL                } from "$projectDir/modules/local/report/dtol"
include { REPORT_GENOMETRAITS        } from "$projectDir/modules/local/report/genometraits"
include { REPORT_SOFTWAREVERSIONS    } from "$projectDir/modules/local/report/softwareversions"
// include { QUARTO as QUARTO_DTOL      } from "$projectDir/modules/local/quarto"
// include { QUARTO as QUARTO_GENESCOPE } from "$projectDir/modules/local/quarto"
include { QUARTO } from "$projectDir/modules/local/quarto"
// include { MULTIQC as MULTIQC_FULL    } from "$projectDir/modules/nf-core/multiqc/main"
// include { MULTIQC as MULTIQC_SUMMARY } from "$projectDir/modules/nf-core/multiqc/main"
include { MULTIQC } from "$projectDir/modules/nf-core/multiqc/main"

workflow ASSEMBLY_REPORT {
    take:
    meta
    logs
    quarto_files // [ meta, notebook, files ]
    versions

    main:
    // DTOL table
    REPORT_DTOL( TOL_SEARCH( meta.map{ meta -> meta.sample.taxid } ).json )

    // Genome traits table
        // Expected vs Observed
        // Haploid Size   // GOAT vs Assembly
        // Haploid Number // GOAT vs HiC
        // Ploidy         // GOAT vs HiC
        // Sample Sex     // GOAT vs HiC
    REPORT_GENOMETRAITS( meta )

    // MultiQC panels from Quarto
    QUARTO( quarto_files )
    // Data profile
    // Input
    // GenomeScope jsons // plots?

    // Quality metrics // MQC summary stats
        // GC // Quast?
        // NG50
        // QV
        // Kmer completeness
        // Busco

    // Contact map

    // Merqury

    // Blobtools

    def mqc_files = logs.mix(
        REPORT_DTOL.out.tsv,
        REPORT_GENOMETRAITS.out.tsv,
        REPORT_SOFTWAREVERSIONS(
            versions
                .mix( QUARTO.out.versions.first() )
                .collect()
        ).yml,
        QUARTO.out.html
    )
    MULTIQC(
        mqc_files.collect().dump(tag:'MultiQC'),
        file("$projectDir/configs/multiqc_summary_report_config.yml", checkIfExists: true),
        params.multiqc.summary_report_extra_config ? file(params.multiqc.summary_report_extra_config, checkIfExists: true) : [],
        []
    )

}