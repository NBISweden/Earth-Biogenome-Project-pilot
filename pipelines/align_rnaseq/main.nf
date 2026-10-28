#!/usr/bin/env nextflow

/*
 * Alignment of RNAseq reads to the primary assembly with STAR.
 *
 * Takes one wide assembly record per assembly build (see types/records.nf). The RNAseq reads are
 * read from `rnaseq`: one element per library, a read pair or a single file.
 */

nextflow.enable.types = true

include { STAR_ALIGN          } from '../../modules/nf-core/star/align/main'
include { STAR_GENOMEGENERATE } from '../../modules/nf-core/star/genomegenerate/main'
include { Assembly            } from '../../types/records'

params {
    assemblies: Channel<Assembly> // Assemblies to align the reads to
}

workflow {

    main:
    index = STAR_GENOMEGENERATE(
        params.assemblies.map { a -> record(id: a.id, meta: a.meta, fasta: a.pri_fasta, gtf: null, star_args: null) } // No valid annotations for draft genome
    )

    aligned = STAR_ALIGN(
        index
            .map { i -> record(id: i.id, sample_id: i.meta.id, meta: i.meta, index: i.index) }
            .join(
                params.assemblies.map { a -> record(id: a.id, rnaseq: a.rnaseq ?: []) },
                by: 'id'
            )
            .flatMap { x ->
                x.rnaseq.collect { reads ->
                    record(
                        id: x.id,
                        meta: x.meta,
                        single_end: reads.size() == 1,
                        reads: reads,
                        index: x.index,
                        gtf: null,                 // No gtf
                        star_ignore_sjdbgtf: true, // Ignore gtf
                        seq_platform: 'Illumina',  // Sequencing Platform
                        seq_center: '',            // Sequencing Center
                        star_args: '--readFilesCommand zcat',
                        output_prefix: x.meta.id
                    )
                }.toList()
            }
    )

    publish:
    aligned = aligned
}

output {
    aligned: Channel<Record> {
        path { r -> "star/${r.meta.assembly.build}" }
    }
}
