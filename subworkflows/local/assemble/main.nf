include { ASSEMBLE_HIFI } from "$projectDir/subworkflows/local/assemble_hifi/main"

workflow ASSEMBLE {
    take:
    hifi_reads // [ meta, hifi_reads ]

    main:
    ch_versions = Channel.empty()

    // TODO: Make strategy check
    ASSEMBLE_HIFI( hifi_reads )
    ch_raw_assemblies = ASSEMBLE_HIFI.out.assemblies

    emit:
    raw_assemblies = ch_raw_assemblies
    versions = ch_versions
}