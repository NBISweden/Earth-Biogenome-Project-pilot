#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO       } from "$projectDir/modules/busco"
include { BLOBTOOLKIT } from "$projectDir/modules/blobtoolkit"

workflow {
    ASSEMBLY_VALIDATION( Channel.fromPath( params.assembly ) )
}

workflow ASSEMBLY_VALIDATION {

     take:
     assembly_ch

    /* Assembly validation workflow:
        - Contamination check ( BLOBTOOLKIT )
        - K-mer spectra check
        - Coverage check ( BLOBTOOLKIT )
        - Gene space check ( BUSCO )
        - Mis assembly signal check
    */
    main:
    BUSCO( assembly )
    BLOBTOOLKIT( assembly )

}
