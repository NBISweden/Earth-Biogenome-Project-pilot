include { MASH_SCREEN } from "$projectDir/modules/nf-core/modules/mash/screen/main"

workflow SCREEN_READS {

    take:
    reads_ch
    mash_screen_db_sketch

    main:
    MASH_SCREEN ( reads_ch, mash_screen_db_sketch )

    // TODO: Visualise results

    emit:
    versions = MASH_SCREEN.out.versions.first()
}