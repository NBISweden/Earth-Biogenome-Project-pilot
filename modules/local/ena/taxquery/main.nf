process ENA_TAXQUERY {
    // directives
    tag { species_name }

    input:
    val species_name

    exec:
    def encoded_species_name = URLEncoder.encode(species_name)
    def ena_response = new URL("https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/${encoded_species_name}").text
    json = new groovy.json.JsonSlurper().parseText(ena_response) as HashMap

    output:
    val json, emit: json
}
