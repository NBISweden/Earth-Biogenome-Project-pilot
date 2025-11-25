include { MITOHIFI_FINDMITOREFERENCE } from "../../../modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI          } from "../../../modules/nf-core/mitohifi/mitohifi/main"
include { OATK                       } from "../../../modules/nf-core/oatk/main"
include { DNADOTPLOT                 } from "../../../modules/local/dnadotplot/main"

workflow ASSEMBLE_ORGANELLES {
    take:
    ch_reads         // Channel: [ meta, reads_path ]
    ch_assemblies    // Channel: [ meta, assembly_map ] - only populated if mode = 'c'
    ch_assembly_mode // String: enum('c','r')
    ch_mito_hmm      // list: [ hmm_files ]
    ch_plastid_hmm   // list: [ hmm_files ]

    main:
    ch_versions = channel.empty()

    // In mode "r", keep reads (ch_assemblies is empty). In mode "c", empty the reads channel and keep assemblies.
    ch_input_data = ch_reads.filter { ch_assembly_mode == 'r' }.mix( ch_assemblies )

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
            input: [ meta, ch_assembly_mode == "c" ? data.pri_fasta : data ] // If "c", use contigs, data = (Channel<Map>]. Else "r", use reads, data = (Channel<Path>).
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

    // Run Oatk assembly
    OATK(
        ch_reads,
        ch_mito_hmm
            .map { _meta, hmm_files ->
                hmm_files
            }
            .ifEmpty( [ [], [], [], [], [] ] ),
        ch_plastid_hmm
            .map { _meta, hmm_files ->
                hmm_files
            }
            .ifEmpty( [ [], [], [], [], [] ] )
    )

    // Mitohifi dot plots: final mitogenome vs. each candidate
    ch_final_vs_candidates = MITOHIFI_MITOHIFI.out.fasta
        .join(MITOHIFI_MITOHIFI.out.all_candidates_fa,
            by: 0
        )
        .flatMap { meta, fasta, candidates ->
            def candidate_list = candidates instanceof List ? candidates : [candidates]
            candidate_list.collect { candidate ->
                [ meta, fasta, candidate ]
            }
        }
    DNADOTPLOT(ch_final_vs_candidates)

    // TODO OATK dot plots:

    // Versions
    ch_versions = ch_versions.mix(
        MITOHIFI_FINDMITOREFERENCE.out.versions.first(),
        MITOHIFI_MITOHIFI.out.versions.first(),
        OATK.out.versions.first(),
        DNADOTPLOT.out.versions.first()
    )

    emit:
    // TODO: emit filtered assembly, or contig list to purge. Also emit from 03_assemble.
    versions = ch_versions
}