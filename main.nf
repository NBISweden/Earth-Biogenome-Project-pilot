#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPARE_INPUT  } from "$projectDir/subworkflows/prepare_input/prepare_input"

include { ASSEMBLY_VALIDATION } from "$projectDir/subworkflows/assembly_validation/assembly_validation"

workflow {

    // The primary workflow for the Earth Biogenome Project Pilot

}

workflow VALIDATE_ASSEMBLIES {

    PREPARE_INPUT ( params.input )
    ASSEMBLY_VALIDATION(
        PREPARE_INPUT.out.assemblies,
        PREPARE_INPUT.out.hifi,
        params.reference ? file( params.reference, checkIfExists: true ) : []
    )

}

workflow.onComplete {
    log.info("""\
    Thank you for using our workflow.

    Results are located in the folder: $params.outdir
    """)
}