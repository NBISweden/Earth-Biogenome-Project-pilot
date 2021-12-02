#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INSPECTOR } from '../../../modules/inspector/inspector'

workflow TEST_INSPECTOR {

    input = [
        [ id: 'test' ],
        [
            file( params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true )
        ],
        file( params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true )
    ]
    // Test: No reference file
    INSPECTOR( input, [] )
}
