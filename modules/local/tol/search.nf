process TOL_SEARCH {
    tag "Taxid: $taxid"
    label 'process_single'

    input:
    val taxid

    output:
    val json, emit: json

    exec:
    println "TOL_SEARCH: Beginning"
    def args = task.ext.args ?: ''
    def response = new URL("https://id.tol.sanger.ac.uk/api/v2/species?taxonomyId=$taxid").text
    json = new java.util.concurrent.ConcurrentHashMap(new groovy.json.JsonSlurper().parseText(response)) // Otherwise returns a LazyMap which causes caching problems.
    println "TOL_SEARCH: End"
}
