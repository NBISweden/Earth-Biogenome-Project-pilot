#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KAT_COMP } from '../../../../modules/kat/kat_comp/kat_comp'

workflow TEST_KAT_COMP {

    input = [
        [ id: 'test' ],
        [
            file( params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true ),
            file( params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true )
        ],
        file( params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true )
    ]
    KAT_COMP( input )
}