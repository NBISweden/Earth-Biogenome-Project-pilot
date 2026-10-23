/*
 * Workflow based around the purge_dups tool
 * https://github.com/dfguan/purge_dups
 */

nextflow.enable.types = true

include { ITERATIVE_PURGEDUPS } from "../../../modules/local/purgedups/main"

record Reads {
    id   : String
    meta : Record
    reads: List<Path>
}

record KmerCoverage {
    id      : String
    kmer_cov: String?
}

record Assembly {
    id       : String  // meta.assembly.build
    meta     : Record  // sample meta data, including meta.assembly
    pri_fasta: Path
    alt_fasta: Path?
    pri_gfa  : Path?
    alt_gfa  : Path?
}

workflow PURGE_DUPLICATES {

    take:
    assemblies: Channel<Assembly>
    hifi      : Channel<Reads>         // One merged reads file per sample
    kmer_cov  : Channel<KmerCoverage>  // Purging waits for the k-mer coverage of the sample

    main:
    // Reads of samples with a k-mer coverage estimate
    hifi_reads = hifi
        .map { h -> record(id: h.id, reads: h.reads) }
        .join(kmer_cov.map { k -> record(id: k.id) }, by: 'id')

    // Combine reads with assemblies and prepare for alignment
    purged = ITERATIVE_PURGEDUPS(
        assemblies
            .map { a -> record(id: a.id, sample_id: a.meta.id, meta: a.meta, assemblies: a.alt_fasta ? [ a.pri_fasta, a.alt_fasta ] : [ a.pri_fasta ]) }
            .join(hifi_reads.map { h -> record(sample_id: h.id, reads: h.reads) }, by: 'sample_id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    reads: x.reads,
                    assemblies: x.assemblies,
                    minimap2_reads_args: "-x map-hifi",
                    minimap2_self_args: "-x asm5 -DP",
                    pbcstat_args: "",
                    calcuts_args: "",
                    hist_plot_args: "",
                    split_fa_args: "",
                    purge_dups_args: "-2",
                    get_seqs_args: "-e",
                    seqkit_args: "-w 100",
                    output_prefix: "${x.meta.id}_${x.meta.assembly.build}"
                )
            }
    )

    assemblies_out = purged.map { p ->
        def sorted = p.purged.toList().toSorted { a, b -> a.name <=> b.name }
        record(
            id: p.id,
            meta: p.meta,
            pri_fasta: sorted[0],
            alt_fasta: sorted.size() > 1 ? sorted[1] : null,
            pri_gfa: null,
            alt_gfa: null
        )
    }

    emit:
    assemblies: Channel<Assembly> = assemblies_out
    logs      : Channel<Path>     = purged.flatMap { p -> p.logs.toList() }
}
