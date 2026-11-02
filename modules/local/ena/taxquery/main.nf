nextflow.enable.types = true

// Taxonomy of a species from the ENA taxonomy service. An exec process, so the result is cached
process ENA_TAXQUERY {
    tag "${species_name}"

    input:
    record(
        id: String,
        species_name: String
    )

    output:
    record(
        id: id,
        json: json
    )

    exec:
    def encoded_species_name = URLEncoder.encode(species_name, 'UTF-8')
    def ena_response = file("https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/${encoded_species_name}").text
    def entries = new groovy.json.JsonSlurper().parseText(ena_response) as List
    assert entries : "Species '${species_name}' not found in the ENA taxonomy"
    def entry = entries[0] as Map
    json = [
        taxId                   : entry['taxId'],
        scientificName          : entry['scientificName'],
        geneticCode             : entry['geneticCode'],
        mitochondrialGeneticCode: entry['mitochondrialGeneticCode'],
        lineage                 : entry['lineage'],
    ]
}
