include { MITOHIFI_FINDMITOREFERENCE } from "../../../modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI          } from "../../../modules/nf-core/mitohifi/mitohifi/main"

workflow ASSEMBLE_ORGANELLES {
    take:
    ch_reads         // Channel: [ meta, reads_path ] - only populated if mode = 'r'
    ch_assemblies    // Channel: [ meta, assembly_map ] - only populated if mode = 'c'
    ch_assembly_mode // String: enum('c','r')

    main:
    ch_versions = Channel.empty()

    // Mix input data channels for mitohifi, only one will be populated
    ch_input_data = ch_reads.mix( ch_assemblies )

    // Attempt mitohifi workflow
    MITOHIFI_FINDMITOREFERENCE( ch_input_data.map { meta, _data -> [ meta, meta.sample.name ] }.unique() )

    // Prepare mitohifi input
    mitohifi_ch = ch_input_data
        .combine(
            MITOHIFI_FINDMITOREFERENCE.out.fasta
                .join(MITOHIFI_FINDMITOREFERENCE.out.gb),
            by: 0
        )
        .multiMap { meta, data, mitofa, mitogb ->
            input: [ meta, ch_assembly_mode == "c" ? data.pri_fasta : data ] // If c, use contigs. Else r, use reads.
            reference: mitofa
            genbank: mitogb
            mito_code: meta.sample.mito_code
        }

    // Run mitochondrial assembly
    MITOHIFI_MITOHIFI(
        mitohifi_ch.input,
        mitohifi_ch.reference,
        mitohifi_ch.genbank,
        ch_assembly_mode,
        mitohifi_ch.mito_code,
    )
    ch_versions = ch_versions.mix(
        MITOHIFI_FINDMITOREFERENCE.out.versions.first(),
        MITOHIFI_MITOHIFI.out.versions.first()
    )

    emit:
    // TODO: emit filtered assembly, or contig list to purge. Also emit from 03_assemble.
    versions = ch_versions
}