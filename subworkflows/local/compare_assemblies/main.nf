nextflow.enable.types = true

include { QUAST    } from "../../../modules/nf-core/quast/main"
include { Assembly } from "../../../types/records"

workflow COMPARE_ASSEMBLIES {

    take:
    assemblies: Channel<Assembly> // The assemblies of all stages. One QUAST run per sample
    reference : Path?             // Reference genome to compare to

    main:
    // Compare the primary assemblies of a sample across the stages, from raw to curated
    quast = QUAST(
        assemblies
            .map { a -> tuple(a.meta.id, a) }
            .groupBy()
            .map { sample_id, results ->
                def stage_order = [ raw: 1, decontaminated: 2, purged: 3, polished: 4, scaffolded: 5, curated: 6 ]
                def sorted = results.toList().toSorted { a, b ->
                    def stage_a = stage_order[a.meta.assembly.stage] ?: 99
                    def stage_b = stage_order[b.meta.assembly.stage] ?: 99
                    stage_a <=> stage_b ?: a.pri_fasta.name <=> b.pri_fasta.name
                }
                record(
                    id: sample_id,
                    meta: record(id: sample_id, sample: sorted[0].meta.sample),
                    consensus: sorted.collect { a -> a.pri_fasta }.toList(),
                    fasta: reference,
                    output_prefix: "${sample_id}_quast_report"
                )
            }
    )

    emit:
    logs : Channel<Path>   = quast.map { q -> q.tsv }
    quast: Channel<Record> = quast
}
