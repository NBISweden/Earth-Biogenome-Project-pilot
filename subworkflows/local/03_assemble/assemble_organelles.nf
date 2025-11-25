include { MITOHIFI_FINDMITOREFERENCE } from "../../../modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI          } from "../../../modules/nf-core/mitohifi/mitohifi/main"
include { OATK                       } from "../../../modules/nf-core/oatk/main"

workflow ASSEMBLE_ORGANELLES {
    take:
    ch_reads         // Channel: [ meta, reads_path ]
    ch_assemblies    // Channel: [ meta, assembly_map ]
    ch_assembly_mode // String: enum('contigs','reads')
    ch_mito_hmm      // list: [ hmm_files ]
    ch_plastid_hmm   // list: [ hmm_files ]

    main:
    ch_versions = channel.empty()

    ch_input_data = ch_assembly_mode == 'contigs' ? ch_assemblies : ch_reads

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
            input: [ meta, ch_assembly_mode == "contigs" ? data.pri_fasta : data ] // contigs ? data = Channel<Map> : data = Channel<Path>.
            reference: mitofa
            genbank: mitogb
            mito_code: meta.sample.mito_code
        }

    // Run mitochondrial assembly
    MITOHIFI_MITOHIFI(
        mitohifi_ch.input,
        mitohifi_ch.reference,
        mitohifi_ch.genbank,
        ch_assembly_mode.charAt(0), // mode: 'c'/'r'
        mitohifi_ch.mito_code,
    )

    // Run Oatk assembly
    OATK(
        ch_input_data.map { meta, data -> tuple(meta, ch_assembly_mode == "contigs" ? data.pri_gfa : data) }, // contigs ? data = Channel<Map> : data = Channel<Path>.
        ch_mito_hmm
            .map { _meta, hmm_files ->
                hmm_files
            },
        ch_plastid_hmm
            .map { _meta, hmm_files ->
                hmm_files
            }
            .ifEmpty( [ [], [], [], [], [] ] )
    )

    // Versions
    ch_versions = ch_versions.mix(
        MITOHIFI_FINDMITOREFERENCE.out.versions.first(),
        MITOHIFI_MITOHIFI.out.versions.first(),
        OATK.out.versions.first()
    )

    emit:
    // TODO: emit filtered assembly, or contig list to purge. Also emit from 03_assemble.
    versions = ch_versions
}