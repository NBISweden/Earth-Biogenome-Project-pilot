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

    // Get mitochondrial reference files
    MITOHIFI_FINDMITOREFERENCE( ch_input_data.map { meta, _data -> [ meta, meta.sample.name ] }.unique() )

    // Run mitohifi assembly
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
    MITOHIFI_MITOHIFI(
        mitohifi_ch.input,
        mitohifi_ch.reference,
        mitohifi_ch.genbank,
        ch_assembly_mode,
        mitohifi_ch.mito_code,
    )

    // Run oatk assembly
    OATK(
        ch_reads,
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

    // Dot plots
    // Reference vs final mitohifi mitogenome channel
    ch_ref_vs_final = MITOHIFI_FINDMITOREFERENCE.out.fasta
        .join(
            MITOHIFI_MITOHIFI.out.fasta,
            by: 0
        )
        .map { meta, ref_fasta, final_fasta -> [ meta, ref_fasta, final_fasta, 'reference' ] }

    // Final mitohifi mitogenome vs. each mitohifi candidate channel
    ch_final_vs_mitohifi = MITOHIFI_MITOHIFI.out.fasta
        .join(
            MITOHIFI_MITOHIFI.out.all_candidates_fa,
            by: 0
        )
        .flatMap { meta, fasta, candidates ->
            def candidate_list = candidates instanceof List ? candidates : [candidates]
            candidate_list.collect { candidate ->
                [ meta, fasta, candidate, 'mitohifi' ]
            }
        }
    // Final mitohifi mitogenome vs. each oatk candidate channel
    ch_final_vs_oatk = MITOHIFI_MITOHIFI.out.fasta
        .map { meta, fasta -> [ [ meta.id, meta.sample], meta, fasta ] }
        .combine(
            OATK.out.mito_fasta
                .map { meta, fasta -> [ [ meta.id, meta.sample ], fasta ] }
                .splitFasta(
                    file: true,
                    by: 1
                ),
            by: 0
        )
        .map { _key, meta, fasta, candidate -> [ meta, fasta, candidate, 'oatk' ] }

    // Mix channels and plot
    ch_dotplot_inputs = ch_ref_vs_final.mix( ch_final_vs_mitohifi, ch_final_vs_oatk )
    DNADOTPLOT( ch_dotplot_inputs )

    // Logs
    MITOHIFI_MITOHIFI.out.contigs_annotations
        .mix (
            MITOHIFI_MITOHIFI.out.stats,
            DNADOTPLOT.out.svg
        )
        .map { _meta, log -> log }
        .set { logs }

    // Versions
    ch_versions = ch_versions.mix(
        MITOHIFI_FINDMITOREFERENCE.out.versions.first(),
        MITOHIFI_MITOHIFI.out.versions.first(),
        OATK.out.versions.first(),
        DNADOTPLOT.out.versions.first()
    )

    emit:
    // TODO: emit filtered assembly, or contig list to purge. Also emit from 03_assemble.
    logs
    versions = ch_versions
}