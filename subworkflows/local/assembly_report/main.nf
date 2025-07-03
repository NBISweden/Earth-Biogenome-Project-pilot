include { TOL_SEARCH                 } from "$projectDir/modules/local/tol/search"
include { REPORT_DTOL                } from "$projectDir/modules/local/report/dtol"
include { REPORT_GENOMETRAITS        } from "$projectDir/modules/local/report/genometraits"
include { REPORT_SOFTWAREVERSIONS    } from "$projectDir/modules/local/report/softwareversions"
include { QUARTO_NOTEBOOK            } from "$projectDir/modules/local/quarto/notebook/main.nf"

workflow ASSEMBLY_REPORT {
    take:
    notebook        // Channel: [ meta:Map, notebook:Path ]
    logs            // Channel: Path
    versions        // Channel: Path
    executed_steps  // Object: Map

    main:
    // DTOL table
    REPORT_DTOL( TOL_SEARCH( notebook.map{ meta, notebook -> meta.sample.tax_id } ).json )

    // Genome traits table
        // Expected vs Observed
        // Haploid Size   // GOAT vs Assembly
        // Haploid Number // GOAT vs HiC
        // Ploidy         // GOAT vs HiC
        // Sample Sex     // GOAT vs HiC
    REPORT_GENOMETRAITS( notebook.map{ meta, notebook -> meta } )
    REPORT_SOFTWAREVERSIONS( versions.toSortedList().dump(tag:'versions', pretty: true) )
    def mqc_files = logs.mix(
        REPORT_DTOL.out.tsv,
        REPORT_GENOMETRAITS.out.tsv,
        REPORT_SOFTWAREVERSIONS.out.yml,
    )

    QUARTO_NOTEBOOK(
        notebook.collect(),
        mqc_files.collect().dump(tag:'MultiQC', pretty: true),
        Channel.value(executed_steps.collect{ k, v -> "$k: ${v}" }.join('\n')).collectFile(),
    )

    emit:
    report = QUARTO_NOTEBOOK.out.html

}