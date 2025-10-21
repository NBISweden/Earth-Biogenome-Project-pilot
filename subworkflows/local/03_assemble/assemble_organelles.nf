include { MITOHIFI_FINDMITOREFERENCE                     } from "../../../modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_CONTIGS } from "../../../modules/nf-core/mitohifi/mitohifi/main"
include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_READS   } from "../../../modules/nf-core/mitohifi/mitohifi/main"

workflow ASSEMBLE_ORGANELLES {
    take:
    ch_input
    ch_input_mode

    main:
    ch_versions = Channel.empty()
    // Find mitochondria
    // TODO: Need to check options to mitohifi modules.
    MITOHIFI_FINDMITOREFERENCE( ch_input.map { meta, _input -> [ meta, meta.sample.name ] }.unique() )

    // Validate input mode and run MITOHIFI accordingly on contigs or reads
    if (["c", "r"].contains(ch_input_mode)) {
        if (ch_input_mode == "c") {
            mitohifi_ch = ch_input
                .combine(
                    MITOHIFI_FINDMITOREFERENCE.out.fasta
                        .join(MITOHIFI_FINDMITOREFERENCE.out.gb),
                    by: 0
                )
                .multiMap { meta, assembly, mitofa, mitogb ->
                    input: [ meta, assembly.pri_fasta ]
                    reference: mitofa
                    genbank: mitogb
                    mito_code: meta.sample.mito_code
                }
            MITOHIFI_MITOHIFI_CONTIGS(
                mitohifi_ch.input,
                mitohifi_ch.reference,
                mitohifi_ch.genbank,
                ch_input_mode,
                mitohifi_ch.mito_code,
            )
            ch_versions = ch_versions.mix(
                MITOHIFI_FINDMITOREFERENCE.out.versions.first(),
                MITOHIFI_MITOHIFI_CONTIGS.out.versions.first()
            )
        } else if (ch_input_mode == "r") {
            mitohifi_ch = ch_input
                .combine(
                    MITOHIFI_FINDMITOREFERENCE.out.fasta
                        .join(MITOHIFI_FINDMITOREFERENCE.out.gb),
                    by: 0
                )
                .multiMap { meta, reads, mitofa, mitogb ->
                    input: [ meta, reads ]
                    reference: mitofa
                    genbank: mitogb
                    mito_code: meta.sample.mito_code
                }
            MITOHIFI_MITOHIFI_READS(
                mitohifi_ch.input,
                mitohifi_ch.reference,
                mitohifi_ch.genbank,
                ch_input_mode,
                mitohifi_ch.mito_code,
            )
            ch_versions = ch_versions.mix(
                MITOHIFI_FINDMITOREFERENCE.out.versions.first(),
                MITOHIFI_MITOHIFI_READS.out.versions.first()
            )
        }
    } else {
        error "Invalid input mode: ${ch_input_mode}"
    }

    emit:
    // TODO: emit filtered assembly, or contig list to purge. Also emit from 03_assemble.
    versions = ch_versions
}