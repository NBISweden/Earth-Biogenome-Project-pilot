include { MITOHIFI_FINDMITOREFERENCE } from "../../../modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI          } from "../../../modules/nf-core/mitohifi/mitohifi/main"

workflow ASSEMBLE_ORGANELLES {
    take:
    ch_raw_assemblies

    main:
    ch_versions = Channel.empty()
    // Find mitochondria
    // TODO: Need to check options to mitohifi modules.
    MITOHIFI_FINDMITOREFERENCE( ch_raw_assemblies.map { meta, _assemblies -> [ meta, meta.sample.name ] }.unique() )
    mitohifi_ch = ch_raw_assemblies
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
    MITOHIFI_MITOHIFI(
        mitohifi_ch.input,
        mitohifi_ch.reference,
        mitohifi_ch.genbank,
        'c',
        mitohifi_ch.mito_code,
    )
    ch_versions = ch_versions.mix(
        MITOHIFI_FINDMITOREFERENCE.out.versions.first(),
        MITOHIFI_MITOHIFI.out.versions.first()
    )

    emit:
    // TODO: emit filtered assembly, or contig list to purge
    versions = ch_versions
}