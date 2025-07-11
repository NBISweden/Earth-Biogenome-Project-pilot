include { getPrimaryAssembly } from "../../../modules/local/functions"
include { QUAST              } from "../../../modules/nf-core/quast/main"

workflow COMPARE_ASSEMBLIES {

    take:
    assembly_ch        // input type: [ meta, AssemblyMap ]

    main:
    QUAST (
        getPrimaryAssembly( assembly_ch )
            .map{ meta, assembly -> tuple(meta.subMap(['id','sample']), assembly) }
            .groupTuple(),
        params.reference ? file( params.reference, checkIfExists: true ) : [],
        []              // No GFF
    )
    versions_ch = QUAST.out.versions
    QUAST.out.tsv
        .map { _meta, tsv -> tsv }
        .set { logs }

    emit:
    versions = versions_ch
    logs

}