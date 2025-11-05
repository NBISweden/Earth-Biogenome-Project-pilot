include { MITOHIFI_FINDMITOREFERENCE } from "../../../modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI          } from "../../../modules/nf-core/mitohifi/mitohifi/main"
include { OATK                       } from "../../../modules/nf-core/oatk/main"

workflow ASSEMBLE_ORGANELLES {
    take:
    ch_assembly_data // Channel: tuple(meta, assembly / reads )
    ch_assembly_mode // String: enum('c','r')
    ch_reads         // Channel: tuple(meta, reads )
    ch_mito_hmm      // list: [ hmm_files ]
    ch_plastid_hmm   // list: [ hmm_files ]

    main:
    ch_versions = Channel.empty()

    // Attempt MitoHiFi workflow
    MITOHIFI_FINDMITOREFERENCE( ch_assembly_data.map { meta, _assembly_data -> [ meta, meta.sample.name ] }.unique() )

    // Prepare mitohifi_ch based on input data type (mode c: contigs, mode r: reads)
    mitohifi_ch = ch_assembly_data
        .combine(
            MITOHIFI_FINDMITOREFERENCE.out.fasta
                .join(MITOHIFI_FINDMITOREFERENCE.out.gb),
            by: 0
        )
        .multiMap { meta, assembly_data, mitofa, mitogb ->
            input: [ meta, ch_assembly_mode == "c" ? assembly_data.pri_fasta : assembly_data ] // If c, use contigs. Else r, use reads.
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

    // Prepare channels for Oatk fallback strategy upon Mitohifi failure
    ch_oatk_input = ch_reads
        .map { meta, reads ->
            // Create join key (ch_reads meta won't match MITOHIFI.out meta in contig mode)
            def meta_key = [ meta.sample ]
            [ meta_key, meta, reads ]
        }
        .join(
            MITOHIFI_MITOHIFI.out.fasta
                .map { meta, mito_fasta ->
                    // Create matching join key
                    def meta_key = [ meta.sample ]
                    [ meta_key, mito_fasta ]
                },
            by: 0,
            remainder: true) // null remainder when there is no Mitohifi fasta output
        .filter { _meta_key, _meta, _reads, mito_fasta ->
            mito_fasta == null
         }
        .map { _meta_key, meta, reads, _mito_fasta ->
            [ meta, reads ]
        }

    // Run Oatk assembly if Mitohifi failed
    OATK(
        ch_oatk_input,
        ch_mito_hmm
            .map { _meta, hmm_files ->
            hmm_files
        },
        ch_plastid_hmm
            .map { _meta, hmm_files ->
            hmm_files
        }
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