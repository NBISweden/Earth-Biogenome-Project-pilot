include { FASTK_FASTK     } from "$projectDir/modules/nf-core/modules/fastk/fastk/main"
include { FASTK_MERGE     } from "$projectDir/modules/nf-core/modules/fastk/merge/main"
include { MERYL_COUNT     } from "$projectDir/modules/nf-core/modules/meryl/count/main"
include { MERYL_UNIONSUM  } from "$projectDir/modules/nf-core/modules/meryl/unionsum/main"
include { MERYL_HISTOGRAM } from "$projectDir/modules/nf-core/modules/meryl/histogram/main"

workflow BUILD_DATABASES {
    
    take:
    fastx_data

    main:
    // Build Meryl or Fastk databases (Containers -> FastK, Conda -> Meryl)
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
        
    MERYL_COUNT ( fastx_data )
    MERYL_UNIONSUM ( MERYL_COUNT.out.meryl_db )
    MERYL_HISTOGRAM ( MERYL_UNIONSUM.out.meryl_db )
    versions_ch = versions_ch.mix(
        MERYL_COUNT.out.versions.first(),
        MERYL_UNIONSUM.out.versions.first(),
        MERYL_HISTOGRAM.out.versions.first()
    )

    emit:
    fastk_histogram = fk_single.hist.mix(FASTK_MERGE.out.hist)
    fastk_ktab      = fk_single.ktab.mix(FASTK_MERGE.out.ktab)
    meryl_histogram = MERYL_HISTOGRAM.out.hist
    meryl_uniondb   = MERYL_UNIONSUM.out.meryl_db
    versions        = versions_ch

}