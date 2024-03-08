include { TOL_SEARCH                 } from "$projectDir/modules/local/tol/search"
include { QUARTO as QUARTO_DTOL      } from "$projectDir/modules/local/quarto"
include { QUARTO as QUARTO_GENESCOPE } from "$projectDir/modules/local/quarto"
include { MULTIQC as MULTIQC_FULL    } from "$projectDir/modules/nf-core/multiqc/main"
include { MULTIQC as MULTIQC_SUMMARY } from "$projectDir/modules/nf-core/multiqc/main"

workflow ASSEMBLY_REPORT {
    take:
    meta
    logs
    versions

    main:
    mqc_files = logs.mix(
        TOL_SEARCH( meta.map{ meta -> meta.sample.taxId } ).json
            .map { json ->
                def data = [
                    tolId: json.species[0].tolIds[0].tolId,
                    species: json.species[0].scientificName,
                    class: json.species[0].taxaClass,
                    order: json.species[0].order
                ].collect { key, value -> "$key\t$value" }.join('\n')
            }
            .collectFile(name:'DToL.tsv', storeDir:"${params.outdir}/reporting"),
        meta.map { meta ->
                """\
                Genome traits\tExpected\tObserved
                Haploid Size\t${meta.sample.genome_size}\tunknown
                Haploid Number\t${meta.sample.haploid_number}\tunknown
                Ploidy\t${meta.sample.ploidy}\tunknown
                """.stripIndent()
            }
        .collectFile(name:'Genome_traits.tsv', storeDir:"${params.outdir}/reporting")
    )


    // Also genome traits // Expected vs Observed
        // Haploid Size   // Assembly vs GOAT
        // Haploid Number //
        // Ploidy         // Ploidyplot vs GOAT
        // Sample Sex     // Chr/ assignment? vs GOAT?

    // Data profile
    // GenomeScope jsons // plots?
    //

    // Quality metrics // MQC summary stats
        // GC // Quast?
        // NG50
        // QV
        // Kmer completeness
        // Busco

    // Contact map

    // Merqury

    // Blobtools

    MULTIQC_SUMMARY(
        mqc_files,
        file("$projectDir/configs/multiqc_summary_report_config.yml", checkIfExists: true),
        params.multiqc.summary_report_extra_config ? file(params.multiqc.summary_report_extra_config, checkIfExists: true) : [],
        [] // BGE logo?
    )

    // MULTIQC_FULL(
    //     mqc_files,
    //     file("$projectDir/configs/multiqc_full_report_config.yml", checkIfExists: true),
    //     params.multiqc.full_report_extra_config ? file(params.multiqc.full_report_extra_config, checkIfExists: true) : [],
    //     [] // NBIS logo?
    // )
}