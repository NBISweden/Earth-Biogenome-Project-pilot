nextflow.enable.types = true

include { QUAST } from "../../../modules/nf-core/quast/main"

record Assembly {
    id       : String
    meta     : Record
    pri_fasta: Path
    alt_fasta: Path?
    pri_gfa  : Path?
    alt_gfa  : Path?
}

workflow COMPARE_ASSEMBLIES {

    take:
    assemblies: Channel<Assembly>
    reference : Path? // Reference genome to compare to

    main:
    quast = QUAST(
        assemblies
            .map { a -> tuple(a.meta.id, a) }
            .groupBy()
            .map { sample_id, results ->
                def sorted = results.toList().toSorted { a, b -> a.pri_fasta.name <=> b.pri_fasta.name }
                record(
                    id: sample_id,
                    meta: record(id: sample_id, sample: sorted[0].meta.sample),
                    consensus: sorted.collect { a -> a.pri_fasta }.toList(),
                    fasta: reference,
                    gff: null,
                    output_prefix: "${sample_id}_quast_report"
                )
            }
    )

    emit:
    logs: Channel<Path> = quast.map { q -> q.tsv }
}
