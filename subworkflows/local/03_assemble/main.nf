include { ASSEMBLE_HIFI       } from "./assemble_hifi"
include { ASSEMBLE_ORGANELLES } from "./assemble_organelles"

workflow ASSEMBLE {
    take:
    hifi_reads              // [ meta, hifi_reads ]
    mito_hmm                // list: [ hmm_files ]
    plastid_hmm             // list: [ hmm_files ]
    nuclear_assembly_mode   // true, false
    organelle_assembly_mode // contigs, reads, or none

    main:
    ch_raw_assemblies = channel.empty()
    ch_logs = channel.empty()
    ch_versions = channel.empty()

    // Nuclear assembly
    if ( nuclear_assembly_mode ) {
        ASSEMBLE_HIFI( hifi_reads )
        ch_raw_assemblies = ASSEMBLE_HIFI.out.assemblies
        ch_logs = ASSEMBLE_HIFI.out.logs
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
            mito_hmm,                // oatk: user-given mito hmm files
            plastid_hmm              // oatk: user-given plastid hmm files
        )
        ch_versions = ch_versions.mix(ASSEMBLE_ORGANELLES.out.versions)
        // TODO: filter organelle contigs from primary assembly
    } // else organelle_assembly_mode == 'none'

    emit:
    // TODO emit organelle assemblies
    raw_assemblies = ch_raw_assemblies
    logs           = ch_logs
    versions       = ch_versions
}