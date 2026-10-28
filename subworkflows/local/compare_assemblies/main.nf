nextflow.enable.types = true

include { QUAST    } from "../../../modules/nf-core/quast/main"
include { Assembly } from "../../../types/records"

workflow COMPARE_ASSEMBLIES {

    take:
    assemblies: Channel<Assembly>
    reference : Path? // Reference genome to compare to

    main:
    // Compare the assemblies of one stage of a sample
    quast = QUAST(
        assemblies
            .map { a -> tuple("${a.meta.id}:${a.meta.assembly.stage}", a) }
            .groupBy()
            .map { _key, results ->
                def sorted = results.toList().toSorted { a, b -> a.pri_fasta.name <=> b.pri_fasta.name }
                def sample_id = sorted[0].meta.id
                def stage = sorted[0].meta.assembly.stage
                record(
                    id: "${sample_id}-${stage}",
                    meta: record(id: sample_id, sample: sorted[0].meta.sample, assembly: [ stage: stage ]),
                    consensus: sorted.collect { a -> a.pri_fasta }.toList(),
                    fasta: reference,
                    gff: null,
                    output_prefix: "${sample_id}_${stage}_quast_report"
                )
            }
    )

    emit:
    logs : Channel<Path>   = quast.map { q -> q.tsv }
    quast: Channel<Record> = quast
}
