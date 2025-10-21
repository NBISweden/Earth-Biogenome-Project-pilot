include { MITOHIFI_FINDMITOREFERENCE                     } from "../../../modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_CONTIGS } from "../../../modules/nf-core/mitohifi/mitohifi/main"
include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_READS   } from "../../../modules/nf-core/mitohifi/mitohifi/main"

workflow ASSEMBLE_ORGANELLES {
    take:
    ch_assembly_data // Channel: tuple(meta, assembly / reads )
    ch_assembly_mode // String: enum('c','r')

    main:
    ch_versions = Channel.empty()
    // Find mitochondria
    // TODO: Need to check options to mitohifi modules.
    MITOHIFI_FINDMITOREFERENCE( ch_assembly_data.map { meta, _assembly_data -> [ meta, meta.sample.name ] }.unique() )

    // Prepare mitohifi_ch based on mode (contigs or reads)
    mitohifi_ch = ch_assembly_data
        .combine(
            MITOHIFI_FINDMITOREFERENCE.out.fasta
                .join(MITOHIFI_FINDMITOREFERENCE.out.gb),
            by: 0
        )
        .multiMap { meta, assembly_data, mitofa, mitogb ->
            input: [ meta, ch_assembly_mode == "c" ? assembly_data.pri_fasta : assembly_data ]
            reference: mitofa
            genbank: mitogb
            mito_code: meta.sample.mito_code
        }
    // Run module alias to allow different ext configs for contigs vs. reads
    if ( ch_assembly_mode == 'c' ) {
        MITOHIFI_MITOHIFI_CONTIGS(
            mitohifi_ch.input,
            mitohifi_ch.reference,
            mitohifi_ch.genbank,
            ch_assembly_mode,
            mitohifi_ch.mito_code,
        )
        ch_versions = ch_versions.mix(
            MITOHIFI_FINDMITOREFERENCE.out.versions.first(),
            MITOHIFI_MITOHIFI_CONTIGS.out.versions.first()
        )
    } else if ( ch_assembly_mode == 'r' ) {
        MITOHIFI_MITOHIFI_READS(
            mitohifi_ch.input,
            mitohifi_ch.reference,
            mitohifi_ch.genbank,
            ch_assembly_mode,
            mitohifi_ch.mito_code,
        )
        ch_versions = ch_versions.mix(
            MITOHIFI_FINDMITOREFERENCE.out.versions.first(),
            MITOHIFI_MITOHIFI_READS.out.versions.first()
        )
    }

    emit:
    // TODO: emit filtered assembly, or contig list to purge. Also emit from 03_assemble.
    versions = ch_versions
}