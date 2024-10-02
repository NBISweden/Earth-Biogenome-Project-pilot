process REPORT_DTOL {
    tag "${tol_search_json.species[0].scientificName}"
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
        tol_search_json.species[0].tolIds[0]?.tolId?: "No ToL ID",
        tol_search_json.species[0].scientificName,
        tol_search_json.species[0].taxaClass,
        tol_search_json.species[0].order
    ].join("\t")
    file("$task.workDir/DToL.tsv").text = tol_table.join('\n')
}
