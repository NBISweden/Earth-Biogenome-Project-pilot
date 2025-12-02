include { ASSEMBLE_HIFI       } from "./assemble_hifi"
include { ASSEMBLE_ORGANELLES } from "./assemble_organelles"

workflow ASSEMBLE {
    take:
    hifi_reads              // Channel: [ meta, hifi_reads ]
    oatkdb                  // Path: /path/to/oatkdb
    nuclear_assembly_mode   // Boolean: assemble chromosomes
    organelle_assembly_mode // Enum: [ contigs, reads, or none ]

    main:
    ch_raw_assemblies = channel.empty()
    ch_logs = channel.empty()
    ch_versions = channel.empty()

    // Nuclear assembly
    if ( nuclear_assembly_mode ) {
        ASSEMBLE_HIFI( hifi_reads )
        ch_raw_assemblies = ASSEMBLE_HIFI.out.assemblies
        ch_logs = ch_logs.mix(ASSEMBLE_HIFI.out.logs)
        ch_versions = ch_versions.mix( ASSEMBLE_HIFI.out.versions )
    } // else nuclear_assembly_mode == false

    // Organelle assembly
    if ( organelle_assembly_mode in ['contigs', 'reads'] ) {
        if ( organelle_assembly_mode == 'contigs' && !nuclear_assembly_mode ) {
            error "Organelle assembly in 'contigs' mode requires 'nuclear_assembly_mode = true'"
        }
        ASSEMBLE_ORGANELLES(
            hifi_reads,              // mitohifi + oatk: [ meta, reads ]
            ch_raw_assemblies,       // mitohifi + oatk: [ meta, assembly_map ]
            organelle_assembly_mode, // mitohifi + oatk: mode
            oatkdb                   // oatk: hmm database
        )
        ch_logs     = ch_logs.mix(ASSEMBLE_ORGANELLES.out.logs)
        ch_versions = ch_versions.mix(ASSEMBLE_ORGANELLES.out.versions)
        // TODO: filter organelle contigs from primary assembly
    } // else organelle_assembly_mode == 'none'

    emit:
    raw_assemblies = ch_raw_assemblies
    logs           = ch_logs
    versions       = ch_versions
}