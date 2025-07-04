process REPORT_DTOL {
    tag "${tol_search_json.species}"
    label 'process_single'

    input:
    val tol_search_json

    output:
    path "DToL.tsv", emit: tsv

    when:
    task.ext.when == null || task.ext.when

    exec:
    def tol_table = []
    tol_table << "ToL ID\tSpecies\tClass\tOrder"
    tol_table << [
        tol_search_json.tol_id,
        tol_search_json.species,
        tol_search_json.class,
        tol_search_json.order,
    ].join("\t")
    file("${task.workDir}/DToL.tsv").text = tol_table.join('\n')
}
