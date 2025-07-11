include { constructAssemblyRecord } from "../../../modules/local/functions"
include { FCSGX_FETCHDB           } from "../../../modules/local/fcsgx/fetchdb"
include { FCSGX_RUNGX             } from "../../../modules/local/fcsgx/rungx"
include { FCSGX_CLEAN             } from "../../../modules/local/fcsgx/clean"
include { deepMergeMaps           } from "../../../modules/local/functions"

workflow DECONTAMINATE {
    take:
    ch_assemblies // [ meta, assembly ]

    main:
    FCSGX_FETCHDB ( params.fcs.database ? Channel.empty() : Channel.fromPath( params.fcs.manifest, checkIfExists: true ) )
    ch_fcs_database = params.fcs.database ? Channel.fromPath( params.fcs.database, checkIfExists: true, type: 'dir' ) : FCSGX_FETCHDB.out.database
    ch_to_screen = ch_assemblies.flatMap { meta, assembly ->
            def updated_meta = deepMergeMaps(meta, [ assembly: [ stage: 'decontaminated', build: "${meta.assembly.assembler}-decontaminated-${meta.assembly.id}" ] ] )
            assembly.alt_fasta ? [
                [ updated_meta + [ haplotype: 0 ], updated_meta.sample.tax_id, assembly.pri_fasta ],
                [ updated_meta + [ haplotype: 1 ], updated_meta.sample.tax_id, assembly.alt_fasta ]
            ] : [
                [ updated_meta + [ haplotype: 0 ], updated_meta.sample.tax_id, assembly.pri_fasta ]
            ]
        }
    FCSGX_RUNGX( ch_to_screen, ch_fcs_database.collect(), params.fcs.ramdisk_path )
    FCSGX_CLEAN(
        ch_to_screen.join( FCSGX_RUNGX.out.fcs_gx_report, by: 0 )
            .map { meta, _taxid, asm, rpt -> [ meta, asm, rpt ] }
    )
    ch_cleaned_assemblies = constructAssemblyRecord(
        FCSGX_CLEAN.out.clean_fasta
            .map { meta, asm -> [ meta.subMap(meta.keySet()-['haplotype']), asm ] },
        true
    )

    FCSGX_FETCHDB.out.versions.mix(
        FCSGX_RUNGX.out.versions,
        FCSGX_CLEAN.out.versions,
    ).set { versions }

    emit:
    assemblies = ch_cleaned_assemblies
    versions
}