include { MITOHIFI_FINDMITOREFERENCE } from "../../../modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI          } from "../../../modules/nf-core/mitohifi/mitohifi/main"
include { OATK_SELECTHMM             } from "../../../modules/local/oatk/selecthmm/main"
include { OATK                       } from "../../../modules/nf-core/oatk/main"
include { DNADOTPLOT                 } from "../../../modules/local/dnadotplot/main"

workflow ASSEMBLE_ORGANELLES {
    take:
    ch_reads         // Channel: [ meta, reads_path ]
    ch_assemblies    // Channel: [ meta, assembly_map ]
    ch_assembly_mode // String: enum('contigs','reads')
    ch_oatkdb        // Path: /path/to/oatkdb

    main:
    ch_versions = channel.empty()

    ch_input_data = ch_assembly_mode == 'contigs' ? ch_assemblies : ch_reads
    ch_meta_data = ch_input_data.map { meta, _data -> [ meta, meta.sample.name ] }.unique()

    // Attempt mitohifi workflow
    MITOHIFI_FINDMITOREFERENCE( ch_meta_data )

    // Run mitohifi assembly
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
    MITOHIFI_MITOHIFI(
        mitohifi_ch.input,
        mitohifi_ch.reference,
        mitohifi_ch.genbank,
        ch_assembly_mode[0], // mode: 'c'/'r'
        mitohifi_ch.mito_code,
    )

    // Run Oatk assembly
    OATK_SELECTHMM ( ch_meta_data.map{ meta, species_name -> tuple(meta, species_name, meta.sample.lineage) }, ch_oatkdb )
    OATK(
        ch_input_data.map { meta, data -> tuple(meta, ch_assembly_mode == "contigs" ? data.pri_gfa : data) }, // contigs ? data = Channel<Map> : data = Channel<Path>.
        OATK_SELECTHMM.out.mito_hmm,
        OATK_SELECTHMM.out.pltd_hmm.ifEmpty( [ [], [], [], [], [] ] )
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
        OATK_SELECTHMM.out.versions.first(),
        OATK.out.versions.first(),
        DNADOTPLOT.out.versions.first()
    )

    emit:
    // TODO: emit filtered assembly, or contig list to purge. Also emit from 03_assemble.
    logs
    versions = ch_versions
}