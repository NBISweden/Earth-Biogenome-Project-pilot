nextflow.enable.types = true

include { STAR_ALIGN          } from "../../../modules/nf-core/star/align/main"
include { STAR_GENOMEGENERATE } from "../../../modules/nf-core/star/genomegenerate/main"

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

workflow ALIGN_RNASEQ {

    take:
    reads      : Channel<Reads>    // RNAseq reads
    assemblies : Channel<Assembly> // TODO: Select assembly stage

    main:
    index = STAR_GENOMEGENERATE(
        assemblies.map { a -> record(id: a.id, meta: a.meta, fasta: a.pri_fasta, gtf: null, star_args: null) } // No valid annotations for draft genome
    )
    aligned = STAR_ALIGN(
        index
            .map { i -> record(id: i.id, sample_id: i.meta.id, meta: i.meta, index: i.index) }
            .join(reads.map { r -> record(sample_id: r.id, single_end: r.meta.single_end, reads: r.reads) }, by: 'sample_id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    single_end: x.single_end,
                    reads: x.reads,
                    index: x.index,
                    gtf: null,                 // No gtf
                    star_ignore_sjdbgtf: true, // Ignore gtf
                    seq_platform: 'Illumina',  // Sequencing Platform
                    seq_center: '',            // Sequencing Center
                    star_args: '--readFilesCommand zcat',
                    output_prefix: x.meta.id
                )
            }
    )

    emit:
    aligned: Channel<Record> = aligned
}
