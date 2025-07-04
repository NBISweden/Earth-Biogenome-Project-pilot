process TOL_SEARCH {
    tag "Taxid: $taxid"
    label 'process_single'

    input:
    val taxid

    when:
    task.ext.when == null || task.ext.when

    exec:
    def args = task.ext.args ?: ''
    def response = new URL("https://id.tol.sanger.ac.uk/api/v2/species?taxonomyId=$taxid").text
    def lazy_json = new groovy.json.JsonSlurper().parseText(response)
    json = [
        tol_id:  lazy_json.species[0]['tolIds'][0]?.tolId?: "No ToL ID",
        species: lazy_json.species[0]['scientificName'],
        class:   lazy_json.species[0]['taxaClass'],
        order:   lazy_json.species[0]['order'],
    ]

    output:
    val json, emit: json
}
