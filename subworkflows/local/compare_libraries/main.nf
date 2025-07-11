include { MERQURYFK_KATCOMP } from "../../../modules/nf-core/merquryfk/katcomp/main"

workflow COMPARE_LIBRARIES {

    take:
    libraryA_hist_ktab_libraryB_hist_ktab // [ meta, libraryA_hist, libraryA_ktab, libraryB_hist, libraryB_ktab ]

    main:
    MERQURYFK_KATCOMP ( libraryA_hist_ktab_libraryB_hist_ktab )
    versions_ch = MERQURYFK_KATCOMP.out.versions.first()

    emit:
    versions = versions_ch
    logs = MERQURYFK_KATCOMP.out.stacked_png.map{ meta, img -> img }
}