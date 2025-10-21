include { ASSEMBLE_HIFI       } from "./assemble_hifi"
include { ASSEMBLE_ORGANELLES } from "./assemble_organelles"

workflow ASSEMBLE {
    take:
    hifi_reads // [ meta, hifi_reads ]
    organelle_assembly_mode // contigs, reads, or none

    main:
    ch_versions = Channel.empty()

    // Nuclear assembly
    // TODO: Make strategy check
    ASSEMBLE_HIFI( hifi_reads )
    ch_raw_assemblies = ASSEMBLE_HIFI.out.assemblies
    ch_versions = ch_versions.mix( ASSEMBLE_HIFI.out.versions )

    // Organelle assembly
    if ( organelle_assembly_mode == 'reads' ) {
        ASSEMBLE_ORGANELLES ( hifi_reads, 'r' )
        ch_versions = ch_versions.mix(ASSEMBLE_ORGANELLES.out.versions)
    } else if ( organelle_assembly_mode == 'contigs' ){
        ASSEMBLE_ORGANELLES ( ch_raw_assemblies, 'c' )
        ch_versions = ch_versions.mix(ASSEMBLE_ORGANELLES.out.versions)
        // TODO: filter organelles from assemblies
    } // else organelle_assembly_mode == 'none'

    emit:
    // TODO emit organelle assemblies
    raw_assemblies = ch_raw_assemblies
    versions       = ch_versions
    logs           = ASSEMBLE_HIFI.out.logs
}