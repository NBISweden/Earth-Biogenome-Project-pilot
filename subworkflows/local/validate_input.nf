#! /usr/bin/env nextflow

/*
Notes:
    - All fastq reads assumed to be hifi for now
    - Assume single set of reads; can groupTuple later if needed
    - Can multimap channels if different output channels are needed
    - Input type currently CSV:
        // Expected format CSV format:
        // ['build_id', 'assembly_path', 'fastq']

*/

nextflow.enable.dsl = 2

workflow {

    VALIDATE_INPUT( params.assembly_csv )

}

workflow VALIDATE_INPUT {

    take:
    asm_csv  // string: path to csv file

    main:
    Channel.fromPath( asm_csv, checkIfExists: true )
        .splitCsv( header: ['build_id', 'assembly_path', 'fastq'], skip: 1 )
        .map { asm -> [ [ id: asm.build_id ], [ file(asm.fastq) ], file(asm.assembly_path) ] }
        .set { assembly_ch }

    emit:
    assemblies = assembly_ch

}
