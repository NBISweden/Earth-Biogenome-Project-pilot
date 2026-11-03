nextflow.enable.types = true

process TOL_SEARCH {
    tag "Taxid: ${taxid}"
    label 'process_single'

    input:
    record(
        id: String,
        taxid: String
    )

    output:
    record(
        id: id,
        json_file: file("${taxid}_tol_info.json"),
        tsv: file("DToL.tsv")
    )

    exec:
    def response = file("https://id.tol.sanger.ac.uk/api/v2/species?taxonomyId=${taxid}").text
    def lazy_json = new groovy.json.JsonSlurper().parseText(response) as Map
    def species = (lazy_json.species as List)[0] as Map
    def tol_ids = species['tolIds'] as List
    def tol_id = tol_ids ? (tol_ids[0] as Map)['tolId'] : "No ToL ID"
    def tol_table = [
        "ToL ID\tSpecies\tClass\tOrder",
        [ tol_id ?: "No ToL ID", species['scientificName'], species['taxaClass'], species['order'] ].join("\t")
    ]
    def jsonBuilder = new groovy.json.JsonBuilder(lazy_json)
    file("${task.workDir}/${taxid}_tol_info.json").text = jsonBuilder.toPrettyString()
    file("${task.workDir}/DToL.tsv").text = tol_table.join('\n')
}
