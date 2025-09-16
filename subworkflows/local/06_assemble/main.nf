include { ASSEMBLE_HIFI } from "./assemble_hifi/main"

workflow ASSEMBLE {
    take:
    hifi_reads // [ meta, hifi_reads ]

    main:
    ch_versions = Channel.empty()

    // TODO: Make strategy check
    ASSEMBLE_HIFI( hifi_reads )
    ch_raw_assemblies = ASSEMBLE_HIFI.out.assemblies
    ch_versions = ch_versions.mix( ASSEMBLE_HIFI.out.versions )

    emit:
    raw_assemblies = ch_raw_assemblies
    versions       = ch_versions
    logs           = ASSEMBLE_HIFI.out.logs
}