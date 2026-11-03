nextflow.enable.types = true

include { TOL_SEARCH          } from "../../../modules/local/tol/search"
include { REPORT_GENOMETRAITS } from "../../../modules/local/report/genometraits"
include { QUARTO_NOTEBOOK     } from "../../../modules/local/quarto/notebook/main.nf"
include { WRITE_FILE          } from "../../../modules/local/report/write_file"

workflow ASSEMBLY_REPORT {
    take:
    metadata     : Channel<Record>   // The sample
    notebook     : Path              // Quarto notebook
    aux_files    : List<Path>        // Notebooks it includes
    logs         : Channel<Path>
    versions     : Value<Path>
    report_params: Map<String,Boolean>

    main:
    // DTOL table
    dtol = TOL_SEARCH(
        metadata.map { m -> record(id: m.id, taxid: "${m.sample.tax_id}") }
    )

    // Genome traits table
        // Expected vs Observed
        // Haploid Size   // GOAT vs Assembly
        // Haploid Number // GOAT vs HiC
        // Ploidy         // GOAT vs HiC
        // Sample Sex     // GOAT vs HiC
    traits = REPORT_GENOMETRAITS(
        metadata.map { m ->
            record(id: m.id, sample_name: m.sample.name, genome_size: "${m.sample.genome_size}", haploid_number: "${m.sample.haploid_number}", ploidy: "${m.sample.ploidy}")
        }
    )

    mqc_files = logs
        .mix(dtol.map { d -> d.tsv })
        .mix(traits.map { t -> t.tsv })
        .collect()

    params_file = WRITE_FILE(
        channel.of( record(id: 'report_params', filename: 'params.yml', content: report_params.keySet().toList().collect { k -> "${k}: ${report_params[k]}" }.join('\n')) )
    )

    report = QUARTO_NOTEBOOK(
        metadata
            .map { m -> record(id: m.id, meta: m) }
            .combine(mqc: mqc_files)
            .combine(versions: versions)
            .combine(params_file: params_file.map { p -> p.file }.collect().map { l -> l.toList()[0] })
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    notebook: notebook,
                    aux_files: aux_files,
                    log_files: x.mqc.toList() + [ x.versions ],
                    report_params: x.params_file
                )
            }
    )

    emit:
    report: Channel<Record> = report
    dtol  : Channel<Record> = dtol
}
