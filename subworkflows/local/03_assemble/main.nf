include { ASSEMBLE_HIFI       } from "./assemble_hifi"
include { ASSEMBLE_ORGANELLES } from "./assemble_organelles"

workflow ASSEMBLE {
    take:
    hifi_reads              // [ meta, hifi_reads ]
    organelle_assembly_mode // contigs, reads, or none
    mito_hmm                // list: [ hmm_files ]
    plastid_hmm             // list: [ hmm_files ]

    main:
    ch_versions = Channel.empty()

    // Nuclear assembly
    // TODO: Make strategy check
    ASSEMBLE_HIFI( hifi_reads )
    ch_raw_assemblies = ASSEMBLE_HIFI.out.assemblies
    ch_versions = ch_versions.mix( ASSEMBLE_HIFI.out.versions )

    // Organelle assembly
    if ( organelle_assembly_mode == 'contigs' ) {
        ASSEMBLE_ORGANELLES(
            Channel.empty(),      // mitohifi: empty reads channel
            ch_raw_assemblies,    // mitohifi: [ meta, assembly_map ]
            'c',                  // mitohifi: mode
            hifi_reads,           // oatk: [ meta, reads ]
            mito_hmm,             // oatk: mito hmm files
            plastid_hmm           // oatk: plastid hmm files
        )
        ch_versions = ch_versions.mix(ASSEMBLE_ORGANELLES.out.versions)
        // TODO: filter organelle contigs from primary assembly
    } else if ( organelle_assembly_mode == 'reads' ) {
        ASSEMBLE_ORGANELLES(
            hifi_reads,           // mitohifi: [ meta, reads ]
            Channel.empty(),      // mitohifi: empty contigs channel
            'r',                  // mitohifi: mode
            hifi_reads,           // oatk: [ meta, reads ]
            mito_hmm,             // oatk: mito hmm files
            plastid_hmm           // oatk: plastid hmm files
        )
        ch_versions = ch_versions.mix(ASSEMBLE_ORGANELLES.out.versions)
    } // else organelle_assembly_mode == 'none'

    emit:
    // TODO emit organelle assemblies
    raw_assemblies = ch_raw_assemblies
    versions       = ch_versions
    logs           = ASSEMBLE_HIFI.out.logs
}