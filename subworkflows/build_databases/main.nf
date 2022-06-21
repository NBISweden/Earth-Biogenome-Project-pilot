include { FASTK_FASTK     } from "$projectDir/modules/nf-core/modules/fastk/fastk/main"
include { MERYL_COUNT     } from "$projectDir/modules/nf-core/modules/meryl/count/main"
include { MERYL_UNIONSUM  } from "$projectDir/modules/nf-core/modules/meryl/unionsum/main"
include { MERYL_HISTOGRAM } from "$projectDir/modules/nf-core/modules/meryl/histogram/main"

workflow BUILD_DATABASES {
    
    take:
    fastx_data

    main:
    // Build Meryl or Fastk databases (Containers -> FastK, Conda -> Meryl)
    FASTK_FASTK ( fastx_data )
    versions_ch = FASTK_FASTK.out.versions.first()
        
    MERYL_COUNT ( fastx_data )
    MERYL_UNIONSUM ( MERYL_COUNT.out.meryl_db )
    MERYL_HISTOGRAM ( MERYL_UNIONSUM.out.meryl_db )
    versions_ch = versions_ch.mix(
        MERYL_COUNT.out.versions.first(),
        MERYL_UNIONSUM.out.versions.first(),
        MERYL_HISTOGRAM.out.versions.first()
    )

    emit:
    fastk_histogram = FASTK_FASTK.out.hist
    fastk_ktab      = FASTK_FASTK.out.ktab
    meryl_histogram = MERYL_HISTOGRAM.out.hist
    meryl_uniondb   = MERYL_UNIONSUM.out.meryl_db
    versions        = versions_ch

}