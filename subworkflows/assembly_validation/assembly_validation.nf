#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VALIDATE_INPUT } from "$projectDir/subworkflows/local/validate_input"

// include { BUSCO       } from "$projectDir/modules/busco"
// include { BLOBTOOLKIT } from "$projectDir/modules/blobtoolkit"
include { INSPECTOR      } from "$projectDir/modules/inspector/inspector"

workflow {

    VALIDATE_INPUT( params.assembly_csv )
    ASSEMBLY_VALIDATION(
        VALIDATE_INPUT().out.assemblies,
        params.reference ? file( params.reference, checkIfExists:true ) : []
    )
}

workflow ASSEMBLY_VALIDATION {

     take:
     assembly_ch   // input type: [ [ id: 'sample_name' ], [ file('path/to/assembly') ] ]
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
    INSPECTOR( assembly_ch, reference_ch )

}
