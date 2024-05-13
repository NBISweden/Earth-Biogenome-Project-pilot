include { combineByMetaKeys   } from "$projectDir/modules/local/functions"
include { getEachAssembly     } from "$projectDir/modules/local/functions"
include { getPrimaryAssembly  } from "$projectDir/modules/local/functions"
include { BUSCO               } from "$projectDir/modules/nf-core/busco/main"
include { MERQURYFK_MERQURYFK } from "$projectDir/modules/local/merquryfk/merquryfk"
// include { INSPECTOR           } from "$projectDir/modules/local/inspector/inspector"

workflow EVALUATE_ASSEMBLY {

    take:
    assembly_ch        // input type: [ meta, [ id:'assemblerX_build1', pri_fasta: '/path/to/primary_asm', alt_fasta: '/path/to/alternate_asm' ] ]
    fastk_db           // input type: [ meta, [ 'path/to/reads.hist' ], [ '/path/to/reads.ktab' ] ]

    main:

    // Kmer consistency check
    MERQURYFK_MERQURYFK (
        combineByMetaKeys (
            fastk_db,
            getEachAssembly(assembly_ch),
            keySet: ['id','sample'],
            meta: 'rhs'
        )
    )

    // Evaluate core gene space coverage
    busco_input = getPrimaryAssembly(assembly_ch)
        .flatMap { meta, asm ->
            if ( meta.settings?.busco?.lineages ) {
                // Use lineages from params.busco.lineages/GOAT.
                meta.settings.busco.lineages.tokenize(',').collect{ [ meta, asm, it ] }
            } else {
                // auto-detect.
                [ [ meta, asm, 'auto' ] ]
            }
        }
        .dump(tag: 'BUSCO')
        .multiMap { meta, asm, line ->
            assembly_ch: [ meta, asm ]
            lineage_ch: line
        }
    BUSCO (
        busco_input.assembly_ch,
        'genome',
        busco_input.lineage_ch,
        params.busco.lineages_db_path ? file( params.busco.lineages_db_path, checkIfExists: true ) : [],
        []
    )

    // TODO: Add Merqury stats to logs
    BUSCO.out.batch_summary
        .set { logs }
    MERQURYFK_MERQURYFK.out.versions.mix(
        BUSCO.out.versions
    ).set { versions }

    emit:
    logs
    versions
}