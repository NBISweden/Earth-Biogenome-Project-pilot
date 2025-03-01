process ENA_TAXQUERY {
    // directives
    tag { species_name }

    input:
    val species_name

    exec:
    def encoded_species_name = URLEncoder.encode(species_name)
    def ena_response = new URL("https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/${encoded_species_name}").text
    def lazy_json = new groovy.json.JsonSlurper().parseText(ena_response)
    json = [
        taxId: lazy_json[0]['taxId'],
        scientificName: lazy_json[0]['scientificName'],
        geneticCode: lazy_json[0]['geneticCode'],
        mitochondrialGeneticCode: lazy_json[0]['mitochondrialGeneticCode'],
        lineage: lazy_json[0]['lineage'],
    ]

    output:
    val json, emit: json
}
