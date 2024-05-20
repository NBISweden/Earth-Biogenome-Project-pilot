include { getPrimaryAssembly } from "$projectDir/modules/local/functions"
include { QUAST              } from "$projectDir/modules/nf-core/quast/main"

workflow COMPARE_ASSEMBLIES {

    take:
    assembly_ch        // input type: [ meta, AssemblyMap ]

    main:
    QUAST (
        getPrimaryAssembly( assembly_ch )
            .groupTuple(),
        params.reference ? file( params.reference, checkIfExists: true ) : [],
        []              // No GFF
    )
    versions_ch = QUAST.out.versions
    QUAST.out.tsv
        .map { it[1] } // Remove meta
        .set { logs }

    emit:
    versions = versions_ch
    logs

}