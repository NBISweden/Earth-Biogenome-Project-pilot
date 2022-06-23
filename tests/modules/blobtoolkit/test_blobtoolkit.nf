#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLOBTOOLKIT } from '../../../../modules/local/blobtoolkit/blobtoolkit'

workflow TEST_BLOBTOOLKIT {

    input = [
        [ id: 'test' ],
        [
            file( params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true ),
        ],
        file( params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true )
    ]
    BLOBTOOLKIT(
        input,
        'eukaryota_odb10',  // busco lineages
        params.busco_lineages_path,
        params.uniprot_path
        params.ncbi_nt_path
        params.ncbi_taxonomy_path
    )
}