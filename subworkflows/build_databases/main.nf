include { FASTK_FASTK     } from "$projectDir/modules/local/fastk/fastk"
include { FASTK_HISTEX    } from "$projectDir/modules/local/fastk/histex"
include { MERYL_COUNT     } from "$projectDir/modules/nf-core/modules/meryl/count/main"
include { MERYL_UNIONSUM  } from "$projectDir/modules/nf-core/modules/meryl/unionsum/main"
include { MERYL_HISTOGRAM } from "$projectDir/modules/nf-core/modules/meryl/histogram/main"

workflow BUILD_DATABASES {
    
    take:
    hifi_data

    main:
    // Build Meryl or Fastk databases (Containers -> FastK, Conda -> Meryl)
    FASTK_FASTK ( hifi_data )
    FASTK_HISTEX ( FASTK_FASTK.out.hist )
    versions_ch = FASTK_FASTK.out.versions.first()
        .mix( FASTK_HISTEX.out.versions.first() )
        
    MERYL_COUNT ( hifi_data )
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
    fastk_histex    = FASTK_HISTEX.out.hist
    meryl_histogram = MERYL_HISTOGRAM.out.hist
    meryl_uniondb   = MERYL_UNIONSUM.out.meryl_db
    versions        = versions_ch

}