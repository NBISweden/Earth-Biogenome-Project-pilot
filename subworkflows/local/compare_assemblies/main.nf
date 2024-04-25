include { getPrimaryAssembly } from "$projectDir/modules/local/functions"
include { QUAST              } from "$projectDir/modules/nf-core/quast/main"

workflow COMPARE_ASSEMBLIES {

    take:
    assembly_ch        // input type: [ [ id: 'sample_name' ], [ id:'assemblerX_build1', primary_asm_path: '/path/to/primary_asm', alternate_asm_path: '/path/to/alternate_asm' ] ]

    main:
    QUAST (
        getPrimaryAssembly( assembly_ch )
            .groupTuple(),
        params.reference ? file( params.reference, checkIfExists: true ) : [],
        []              // No GFF
    )
    versions_ch = QUAST.out.versions

    emit:
    versions = versions_ch

}