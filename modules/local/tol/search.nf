import groovy.json.JsonSlurper

process TOL_SEARCH {
    tag "$taxid"
    label 'process_single'

    input:
    val taxid

    output:
    val json, emit: json

    when:
    task.ext.when == null || task.ext.when

    exec:
    def args = task.ext.args ?: ''
    def response = new URL("https://id.tol.sanger.ac.uk/api/v2/species?taxonomyId=$taxid").text
    json = new JsonSlurper().parseText(response)
}
