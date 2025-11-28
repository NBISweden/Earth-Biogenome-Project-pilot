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
        ch_logs = ch_logs.mix(ASSEMBLE_HIFI.out.logs)
        ch_versions = ch_versions.mix( ASSEMBLE_HIFI.out.versions )
    } // else nuclear_assembly_mode == false

    // Organelle assembly
    if ( organelle_assembly_mode == 'contigs' && nuclear_assembly_mode ) {
        ASSEMBLE_ORGANELLES(
            hifi_reads,           // oatk: [ meta, reads ]
            ch_raw_assemblies,    // mitohifi: [ meta, assembly_map ]
            'c',                  // mitohifi: mode
            mito_hmm,             // oatk: mito hmm files
            plastid_hmm           // oatk: plastid hmm files
        )
        ch_logs     = ch_logs.mix(ASSEMBLE_ORGANELLES.out.logs)
        ch_versions = ch_versions.mix(ASSEMBLE_ORGANELLES.out.versions)
        // TODO: filter organelle contigs from primary assembly
    } else if ( organelle_assembly_mode == 'reads' ) {
        ASSEMBLE_ORGANELLES(
            hifi_reads,           // mitohifi & oatk: [ meta, reads ]
            channel.empty(),      // mitohifi: empty contigs channel
            'r',                  // mitohifi: mode
            mito_hmm,             // oatk: mito hmm files
            plastid_hmm           // oatk: plastid hmm files
        )
        ch_logs     = ch_logs.mix(ASSEMBLE_ORGANELLES.out.logs)
        ch_versions = ch_versions.mix(ASSEMBLE_ORGANELLES.out.versions)
    } else if ( organelle_assembly_mode == 'contigs' && !nuclear_assembly_mode ) {
        error "Organelle assembly in 'contigs' mode requires 'nuclear_assembly_mode = true'"
    } // else organelle_assembly_mode == 'none'

    emit:
    // TODO emit organelle assemblies
    raw_assemblies = ch_raw_assemblies
    logs           = ch_logs
    versions       = ch_versions
}