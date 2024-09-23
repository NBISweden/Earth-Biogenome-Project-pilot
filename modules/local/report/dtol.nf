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
    def tol_table = [
        tolId: tol_search_json.species[0].tolIds[0]?.tolId?: "No ToL ID",
        species: tol_search_json.species[0].scientificName,
        class: tol_search_json.species[0].taxaClass,
        order: tol_search_json.species[0].order
    ].collect { key, value -> "$key\t$value" }.join('\n')
    file("$task.workDir/DToL.tsv").text = tol_table
}
