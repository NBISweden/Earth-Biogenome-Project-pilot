include { MERQURYFK_KATCOMP } from "$projectDir/modules/nf-core/modules/merquryfk/katcomp/main"

workflow COMPARE_LIBRARIES {

    take:
    libraryA_hist_ktab_libraryB_hist_ktab // [ meta, libraryA_hist, libraryA_ktab, libraryB_hist, libraryB_ktab ]

    main:
    MERQURYFK_KATCOMP ( libraryA_hist_ktab_libraryB_hist_ktab )
    versions_ch = MERQURYFK_KATCOMP.out.versions.first()

    emit:
    versions = versions_ch    

}