#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

// include { VALIDATE_INPUT } from "$projectDir/subworkflows/local/validate_input"
include { PREPARE_INPUT  } from "$projectDir/subworkflows/prepare_input/prepare_input"

// include { BUSCO       } from "$projectDir/modules/busco"
// include { BLOBTOOLKIT } from "$projectDir/modules/blobtoolkit"
include { QUAST          } from "$projectDir/modules/nf-core/modules/quast/main"
include { INSPECTOR      } from "$projectDir/modules/local/inspector/inspector"

workflow {

    PREPARE_INPUT ( params.input )
    ASSEMBLY_VALIDATION(
        PREPARE_INPUT.out.assemblies,
        PREPARE_INPUT.out.hifi,
        params.reference ? file( params.reference, checkIfExists:true ) : []
    )
}

workflow ASSEMBLY_VALIDATION {

     take:
     assembly_ch   // input type: [ [ id: 'sample_name' ], [ id:'assemblerX_build1', path:'/path/to/assembly' ] ]
     reads_ch      // input type: [ [ id: 'sample_name' ], [ file('path/to/reads') ] ]
     reference_ch  // optional: file( reference_genome ) for comparison

    /* Assembly validation workflow:
        - Contamination check ( BLOBTOOLKIT )
        - K-mer spectra check
        - Coverage check ( BLOBTOOLKIT )
        - Gene space check ( BUSCO )
        - Mis assembly signal check
    */
    main:
    // BUSCO( assembly )
    // BLOBTOOLKIT( assembly )
    QUAST (
        assembly_ch.map { sample, assembly -> assembly.path }
            .collect(),
        reference_ch,
        [], // gff
        reference_ch, // true / false to use reference_ch
        []
    )
    INSPECTOR ( reads_ch.join( assembly_ch ), reference_ch )

}
