include { ENA_TAXQUERY      } from "../../../modules/local/ena/taxquery/main"
include { GOAT_TAXONSEARCH  } from "../../../modules/nf-core/goat/taxonsearch/main"

workflow FETCH_SAMPLE_METADATA {
    take:
    ch_sample // Map: [ sample: [ name: species, ... ] ]

    main:
    // Prepare channel
    ch_with_id = ch_sample.map { yaml -> yaml +
        [
            id: yaml.sample.name.replace(" ", "_"),
            sample: new SampleInfo(
                sampleName: yaml.sample.name,
                taxId: yaml.sample.tax_id,
                geneticCode: yaml.sample.genetic_code,
                mitoCode: yaml.sample.mito_code,
                domain: yaml.sample.domain,
                genomeSize: yaml.sample.genome_size,
                haploidNumber: yaml.sample.haploid_number,
                ploidy: yaml.sample.ploidy,
            )
        ]
    }
    // Get species metadata from ENA
    ENA_TAXQUERY(
        ch_with_id
            .filter { meta -> !meta.sample.hasComponents(['taxId', 'geneticCode', 'mitoCode', 'domain']) } // Skip process if all is user defined
            .map { yaml -> yaml.sample.sampleName() }
    )
    // ch_with_ena = [ sample: [ name: ... ] ] x [ taxid: ... , geneticCode: ..., ...  ]
    ch_with_ena = ch_with_id
        .combine( ENA_TAXQUERY.out.json.ifEmpty( [[]] ) )
        .map { input, ena ->
            ena ? input.deepMerge(
                [
                    sample: input.sample.copyWith(
                            taxId: input.sample.taxId() ?: ena.taxId.toInteger(),
                            geneticCode: input.sample.geneticCode() ?: ena.geneticCode.toInteger(),
                            mitoCode: input.sample.mitoCode() ?: ena.mitochondrialGeneticCode.toInteger(),
                            domain: input.sample.domain() ?: ena.lineage.tokenize(';').head(),
                    )
                ]
            ) : input
        }
    // def sample_keys_goat = [ 'genome_size', 'haploid_number', 'ploidy' ]
    GOAT_TAXONSEARCH(
        ch_with_ena
            .filter { meta -> meta.sample.isEukaryota() } // Eukaryotes only
            .filter { meta -> !meta.sample.hasComponents([ 'genomeSize', 'haploidNumber', 'ploidy' ]) || !params.busco.lineages } // Skip if all is user defined
            .map { meta -> tuple( meta, meta.sample.sampleName(), [] ) }
    )
    // ch_metadata = [ sample: [ name: ..., tax_id: ... ] ] x [ [ sample: [ name: ... ] ], tsv ]
    ch_metadata = ch_with_ena
        .combine( GOAT_TAXONSEARCH.out.taxonsearch.ifEmpty( [[],[]] ) )
        .map { input, _goat_meta, goat_tsv ->
            def updated_metadata = input
            if( goat_tsv ){
                def species = goat_tsv.splitCsv( sep:"\t", header: true ).find { tsv -> tsv.scientific_name == input.sample.sampleName() }
                assert species != null : "GOAT_TAXONSEARCH failed to retrieve species information"
                def busco_lineages = goat_tsv.splitCsv( sep:"\t", header: true )
                    .findAll { row -> row.odb10_lineage }
                    .collect { row -> row.odb10_lineage }
                    .join(',')
                updated_metadata = updated_metadata.deepMerge(
                    [
                        sample: input.sample.copyWith(
                                genomeSize: input.sample.genomeSize() ?: species.genome_size.toInteger(),
                                haploidNumber: input.sample.haploidNumber() ?: species.haploid_number.toInteger(),
                                ploidy: input.sample.ploidy() ?: species.ploidy.toInteger(),
                        ),
                        settings: [ busco: [ lineages: params.busco.lineages?: busco_lineages ] ]
                    ]
                )
            } else {
                updated_metadata = updated_metadata.deepMerge(
                    [
                        settings: [ busco: [ lineages: params.busco.lineages ] ]
                    ]
                )
            }
            updated_metadata
        }

    emit:
    metadata = ch_metadata
}