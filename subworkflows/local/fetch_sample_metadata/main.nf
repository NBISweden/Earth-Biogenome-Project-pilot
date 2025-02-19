include { ENA_TAXQUERY      } from "../../../modules/local/ena/taxquery/main"
include { GOAT_TAXONSEARCH  } from "../../../modules/nf-core/goat/taxonsearch/main"

workflow FETCH_SAMPLE_METADATA {
    take:
    ch_sample // Map: [ sample: [ name: species, ... ] ]

    main:
    // Metadata to retrieve
    def sample_keys_ena = [ 'tax_id', 'genetic_code', 'mito_code', 'domain' ]
    def sample_keys_goat = [ 'genome_size', 'haploid_number', 'ploidy' ]
    // Prepare channel
    ch_with_id = ch_sample.map { yaml -> yaml + [id: yaml.sample.name.replace(" ", "_") ] }
    // Get species metadata from ENA
    ENA_TAXQUERY(
        ch_with_id
            .filter { yaml -> !yaml.sample.keySet().containsAll( sample_keys_ena ) } // Skip process if all is user defined
            .map { yaml -> yaml.sample.name }
    )
    // ch_with_ena = [ sample: [ name: ... ] ] x [ taxid: ... , geneticCode: ..., ...  ]
    ch_with_ena = ch_with_id
        .combine( ENA_TAXQUERY.out.json.ifEmpty( [[]] ) )
        .map { input, ena ->
            ena ? input.deepMerge(
                [
                    sample: [
                        tax_id: input.sample.tax_id ?: ena.taxId,
                        genetic_code: input.sample.genetic_code ?: ena.geneticCode,
                        mito_code: input.sample.mito_code ?: ena.mitochondrialGeneticCode,
                        domain: input.sample.domain ?: ena.lineage.tokenize(';').head(),
                    ]
                ]
            ) : input
        }
    GOAT_TAXONSEARCH(
        ch_with_ena
            .filter { meta -> meta.sample.domain == "Eukaryota" } // Eukaryotes only
            .filter { meta -> !meta.sample.keySet().containsAll( sample_keys_goat ) || !params.busco.lineages } // Skip if all is user defined
            .map { meta -> tuple( meta, meta.sample.name, [] ) }
    )
    // ch_metadata = [ sample: [ name: ..., tax_id: ... ] ] x [ [ sample: [ name: ... ] ], tsv ]
    ch_metadata = ch_with_ena
        .combine( GOAT_TAXONSEARCH.out.taxonsearch.ifEmpty( [[],[]] ) )
        .map { input, _goat_meta, goat_tsv ->
            def updated_metadata = input
            if( goat_tsv ){
                def species = goat_tsv.splitCsv( sep:"\t", header: true ).find { tsv -> tsv.scientific_name == input.sample.name }
                assert species != null : "GOAT_TAXONSEARCH failed to retrieve species information"
                def busco_lineages = goat_tsv.splitCsv( sep:"\t", header: true )
                    .findAll { row -> row.odb10_lineage }
                    .collect { row -> row.odb10_lineage }
                    .join(',')
                updated_metadata = updated_metadata.deepMerge(
                    [
                        sample: [
                            genome_size: input.sample.genome_size ?: species.genome_size,
                            haploid_number: input.sample.haploid_number ?: species.haploid_number,
                            ploidy: input.sample.ploidy ?: species.ploidy,
                        ],
                        settings: [ busco: [ lineages: params.busco.lineages?: busco_lineages ] ]
                    ]
                )
            }
            updated_metadata
        }

    emit:
    metadata = ch_metadata
}