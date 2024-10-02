include { combineByMetaKeys   } from "$projectDir/modules/local/functions"
include { getEachAssembly     } from "$projectDir/modules/local/functions"
include { getPrimaryAssembly  } from "$projectDir/modules/local/functions"
include { BUSCO               } from "$projectDir/modules/nf-core/busco/main"
include { MERQURYFK_MERQURYFK } from "$projectDir/modules/local/merquryfk/merquryfk"
include { MERQURY             } from "$projectDir/modules/nf-core/merqury/main"
// include { INSPECTOR           } from "$projectDir/modules/local/inspector/inspector"

workflow EVALUATE_ASSEMBLY {

    take:
    assembly_ch        // input type: [ meta, [ id:'assemblerX_build1', pri_fasta: '/path/to/primary_asm', alt_fasta: '/path/to/alternate_asm' ] ]
    fastk_db           // input type: [ meta, [ 'path/to/reads.hist' ], [ '/path/to/reads.ktab' ] ]
    meryl_db           // input type: [ meta, 'path/to/reads.union.meryldb' ]

    main:

    // Kmer consistency check
    MERQURYFK_MERQURYFK (
        combineByMetaKeys (
            fastk_db,
            getEachAssembly(assembly_ch),
            keySet: ['id','sample'],
            meta: 'rhs'
        ) // [ meta, hist, ktab, assembly ]
    )

    MERQURY (
        combineByMetaKeys (
            meryl_db,
            getEachAssembly(assembly_ch),
            keySet: ['id','sample'],
            meta: 'rhs'
        ) // [ meta, meryldb, assembly ]
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

    BUSCO.out.short_summaries_txt
        .mix(
            MERQURYFK_MERQURYFK.out.stats,
            MERQURYFK_MERQURYFK.out.qv,
            // MERQURYFK_MERQURYFK.out.assembly_qv, // Contig names are missing in first column
            MERQURYFK_MERQURYFK.out.spectra_cn_st_png,
            MERQURYFK_MERQURYFK.out.spectra_asm_st_png,
            MERQURYFK_MERQURYFK.out.false_duplications,
            MERQURY.out.scaffold_qv,
            MERQURY.out.spectra_cn_st_png,
            MERQURY.out.spectra_asm_st_png,
        )
        .map { meta, file -> file }
        .set { logs }
    MERQURYFK_MERQURYFK.out.versions.first().mix(
        MERQURY.out.versions.first(),
        BUSCO.out.versions.first()
    ).set { versions }

    emit:
    logs
    versions
}