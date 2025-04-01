include { MITOHIFI_FINDMITOREFERENCE } from "$projectDir/modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI          } from "$projectDir/modules/nf-core/mitohifi/mitohifi/main"

workflow ASSEMBLE_ORGANELLES {
    take:
    ch_raw_assemblies

    main:
    ch_versions = Channel.empty()
    // Find mitochondria
    // TODO: Need to check options to mitohifi modules.
    MITOHIFI_FINDMITOREFERENCE( ch_raw_assemblies.map { meta, assemblies -> [ meta, meta.sample.sampleName() ] }.unique() )
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
            mito_code: meta.sample.mitoCode()
        }
    MITOHIFI_MITOHIFI(
        mitohifi_ch.input,
        mitohifi_ch.reference,
        mitohifi_ch.genbank,
        'c',
        mitohifi_ch.mito_code,
    )

    emit:
    // TODO: emit filtered assembly, or contig list to purge
    versions = ch_versions
}