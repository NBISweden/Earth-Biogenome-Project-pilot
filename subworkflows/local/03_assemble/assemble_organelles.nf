nextflow.enable.types = true

include { MITOHIFI_FINDMITOREFERENCE } from "../../../modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI          } from "../../../modules/nf-core/mitohifi/mitohifi/main"
include { OATK_SELECTHMM             } from "../../../modules/local/oatk/selecthmm/main"
include { OATK                       } from "../../../modules/nf-core/oatk/main"
include { DNADOTPLOT                 } from "../../../modules/local/dnadotplot/main"

record Reads {
    id   : String
    meta : Record
    reads: List<Path>
}

record Assembly {
    id       : String
    meta     : Record
    pri_fasta: Path
    alt_fasta: Path?
    pri_gfa  : Path?
    alt_gfa  : Path?
}

record OrganelleInput {
    id         : String
    meta       : Record
    mito_input : Path  // mitohifi: contigs fasta or reads
    oatk_input : Path? // oatk: contigs graph or reads
}

record DotplotInput {
    id        : String
    meta      : Record
    ref       : Path
    query     : Path
    experiment: String
}

def dotplotInput(id: String, meta: Record, ref: Path, query: Path, experiment: String) -> DotplotInput {
    return record(id: id, meta: meta, ref: ref, query: query, experiment: experiment)
}

workflow ASSEMBLE_ORGANELLES {
    take:
    reads        : Channel<Reads>
    assemblies   : Channel<Assembly>
    assembly_mode: String // enum('contigs','reads')
    oatkdb       : Path?  // /path/to/oatkdb. Oatk is skipped without it.

    main:
    input_data = assembly_mode == 'contigs'
        ? assemblies.map { a -> record(id: a.id, meta: a.meta, mito_input: a.pri_fasta, oatk_input: a.pri_gfa) }
        : reads.map { r -> record(id: r.id, meta: r.meta, mito_input: r.reads[0], oatk_input: r.reads[0]) }

    // Attempt mitohifi workflow
    reference = MITOHIFI_FINDMITOREFERENCE(
        input_data.map { x -> record(id: x.id, meta: x.meta, species: x.meta.sample.name) }
    ).filter { r -> r.fasta != null && r.gb != null }

    // Run mitohifi assembly
    mitohifi = MITOHIFI_MITOHIFI(
        input_data
            .map { x -> record(id: x.id, mito_input: x.mito_input) }
            .join(reference, by: 'id')
            .map { j ->
                record(
                    id: j.id,
                    meta: j.meta,
                    input: j.mito_input,
                    ref_fa: j.fasta,
                    ref_gb: j.gb,
                    input_mode: assembly_mode == 'contigs' ? 'c' : 'r',
                    mito_code: "${j.meta.sample.mito_code}",
                    output_prefix: assembly_mode == 'contigs' ? "${j.meta.id}-${j.meta.assembly.build}" : "${j.meta.id}"
                )
            }
    )
    mitohifi_final = mitohifi.filter { m -> m.fasta != null }

    // Run Oatk assembly
    selecthmm = OATK_SELECTHMM(
        input_data
            .filter { _x -> oatkdb != null }
            .map { x -> record(id: x.id, meta: x.meta, species: x.meta.sample.name, lineage: x.meta.sample.lineage, oatkdb: oatkdb) }
    )
    oatk = OATK(
        input_data
            .map { x -> record(id: x.id, meta: x.meta, reads: x.oatk_input, oatk_args: x.oatk_input.name.endsWith(".gfa") ? "-G" : "") }
            .join(
                selecthmm
                    .filter { s -> s.mito_hmm }
                    .map { s -> record(id: s.id, mito_hmm: s.mito_hmm.toList(), pltd_hmm: s.pltd_hmm.toList()) },
                by: 'id'
            )
    )
    oatk_final = oatk.filter { o -> o.mito_fasta != null }

    // Dot plots
    // Reference vs final mitohifi mitogenome
    ref_vs_final = reference
        .map { r -> record(id: r.id, ref: r.fasta) }
        .join(mitohifi_final, by: 'id')
        .map { j -> dotplotInput(j.id, j.meta, j.ref, j.fasta, 'reference') }

    // Final mitohifi mitogenome vs. each mitohifi candidate
    final_vs_mitohifi = mitohifi_final.flatMap { m ->
        m.all_candidates_fa.collect { candidate -> dotplotInput(m.id, m.meta, m.fasta, candidate, 'mitohifi') }.toList()
    }

    // Final mitohifi mitogenome vs. each oatk candidate
    final_vs_oatk = mitohifi_final
        .map { m -> record(id: m.id, meta: m.meta, ref: m.fasta) }
        .join(
            oatk_final.flatMap { o ->
                o.mito_fasta.splitFasta(file: true, by: 1).collect { candidate -> record(id: o.id, query: candidate) }.toList()
            },
            by: 'id'
        )
        .map { j -> dotplotInput(j.id, j.meta, j.ref, j.query, 'oatk') }

    // Mix channels and plot
    dotplots = DNADOTPLOT(
        ref_vs_final
            .mix(final_vs_mitohifi)
            .mix(final_vs_oatk)
            .map { d -> record(id: d.id, meta: d.meta, ref: d.ref, query: d.query, experiment: d.experiment, dnadotplot_args: "--revcompl --window 40") }
    )

    // Mitochondrial assemblies
    mito_assemblies = mitohifi_final
        .map { m -> record(id: m.id, meta: m.meta, fasta: m.fasta) }
        .mix(oatk_final.map { o -> record(id: o.id, meta: o.meta, fasta: o.mito_fasta) })

    // Logs
    logs = mitohifi
        .map { m -> m.contigs_annotations }
        .mix(mitohifi.map { m -> m.stats })
        .mix(dotplots.map { d -> d.svg })
        .filter { f -> f != null }

    emit:
    // TODO: emit filtered assembly, or contig list to purge. Also emit from 03_assemble.
    mito_assemblies: Channel<Record> = mito_assemblies
    logs           : Channel<Path>   = logs
    reference      : Channel<Record> = reference
    mitohifi       : Channel<Record> = mitohifi
    oatk           : Channel<Record> = oatk
    dotplots       : Channel<Record> = dotplots
}
