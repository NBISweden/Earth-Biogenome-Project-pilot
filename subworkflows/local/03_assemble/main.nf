include { ASSEMBLE_HIFI       } from "./assemble_hifi"
include { ASSEMBLE_ORGANELLES } from "./assemble_organelles"

workflow ASSEMBLE {
    take:
    hifi_reads // [ meta, hifi_reads ]
    nuclear_assembly_mode // true, false
    organelle_assembly_mode // contigs, reads, or none

    main:
    ch_raw_assemblies = Channel.empty()
    ch_logs = Channel.empty()
    ch_versions = Channel.empty()

    // Nuclear assembly
    if ( nuclear_assembly_mode ) {
        ASSEMBLE_HIFI( hifi_reads )
        ch_raw_assemblies = ASSEMBLE_HIFI.out.assemblies
        ch_logs = ASSEMBLE_HIFI.out.logs
        ch_versions = ch_versions.mix( ASSEMBLE_HIFI.out.versions )
    } // else nuclear_assembly_mode == false

    // Organelle assembly
    if ( organelle_assembly_mode == 'reads' ) {
        ASSEMBLE_ORGANELLES ( hifi_reads, 'r' )
        ch_versions = ch_versions.mix(ASSEMBLE_ORGANELLES.out.versions)
    } else if ( organelle_assembly_mode == 'contigs' && nuclear_assembly_mode ) {
        ASSEMBLE_ORGANELLES ( ch_raw_assemblies, 'c' )
        ch_versions = ch_versions.mix(ASSEMBLE_ORGANELLES.out.versions)
        // TODO: filter organelles from assemblies
    } else if (organelle_assembly_mode == 'contigs' && !nuclear_assembly_mode ) {
        error "ERROR: Organelle assembly with contigs also requires nuclear assembly."
    } // else organelle_assembly_mode == 'none'

    emit:
    // TODO emit organelle assemblies
    raw_assemblies = ch_raw_assemblies
    logs           = ch_logs
    versions       = ch_versions
}