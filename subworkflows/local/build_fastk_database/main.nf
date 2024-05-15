include { FASTK_FASTK     } from "$projectDir/modules/nf-core/fastk/fastk/main"
include { FASTK_MERGE     } from "$projectDir/modules/nf-core/fastk/merge/main"

workflow BUILD_FASTK_DATABASE {

    take:
    fastx_data

    main:
    FASTK_FASTK ( fastx_data )
    fkdb_ch = FASTK_FASTK.out.hist.groupTuple()
        .join( FASTK_FASTK.out.ktab.groupTuple(), remainder: true )
        .join( FASTK_FASTK.out.prof.groupTuple(), remainder: true )
        .map { meta, hist, ktab, prof -> [meta.findAll { it.key != 'single_end' } , hist, ktab ? ktab.flatten() : [] , prof ? prof.flatten() : [] ] }
        .branch { meta, hist, ktab, prof ->
            single_hist: hist.size() == 1
            multi_hist : hist.size() > 1
        }
    FASTK_MERGE ( fkdb_ch.multi_hist )
    fk_single = fkdb_ch.single_hist.multiMap { meta, hist, ktab, prof ->
        hist: [ meta, hist ]
        ktab: [ meta, ktab ]
        prof: [ meta, prof ]
    }
    versions_ch = FASTK_FASTK.out.versions.first().mix( FASTK_MERGE.out.versions.first() )

    emit:
    fastk_histogram = fk_single.hist.mix(FASTK_MERGE.out.hist)
    fastk_ktab      = fk_single.ktab.mix(FASTK_MERGE.out.ktab)
    fastk_hist_ktab = fk_single.hist.mix(FASTK_MERGE.out.hist).join(fk_single.ktab.mix(FASTK_MERGE.out.ktab))
    versions        = versions_ch

}