nextflow.enable.types = true

include { FASTP           } from '../../../modules/nf-core/fastp/main'
include { SAMTOOLS_IMPORT } from '../../../modules/nf-core/samtools/import/main'
include { SAMTOOLS_INDEX  } from '../../../modules/nf-core/samtools/index/main'

include { Reads } from "../../../types/records"

// Longest common prefix of two file names, without trailing separators or a trailing _R
def commonPrefix(first: String, last: String) -> String {
    def minLength = Math.min(first.length(), last.length())
    def prefixLength = (0..<minLength).takeWhile { i -> first[i] == last[i] }.size()
    def prefix = first.substring(0, prefixLength).replaceAll(/[^a-zA-Z0-9]+$/, '')
    return prefix.endsWith("_R") ? prefix.substring(0, prefix.length() - 2) : prefix
}

workflow CONVERT_FASTQ_CRAM {
    take:
    samples        : Channel<Record> // With the Hi-C read pairs in `hic`
    bool_trim_fastq: Boolean         // Trim with Fastp
    hic_type       : String          // Library type, recorded in the CRAM read group

    main:
    // One record per Hi-C read pair, numbered
    ch_fastq = samples.flatMap { s ->
        ( s.hic ?: [] ).withIndex().collect { pair, index ->
            record(id: s.id, meta: s.meta + record(single_end: false, pair_id: index), reads: pair)
        }.toList()
    }

    // Filter
    trimmed = FASTP(
        ch_fastq
            .filter { _r -> bool_trim_fastq }
            .map { r ->
                record(
                    id: r.id,
                    meta: r.meta,
                    single_end: r.meta.single_end,
                    reads: r.reads,
                    adapter_fasta: null,
                    discard_trimmed_pass: false,
                    save_trimmed_fail: false,
                    save_merged: false,
                    fastp_args: [
                        "--disable_adapter_trimming",
                        "--disable_quality_filtering",
                        "--disable_trim_poly_g",
                        "--disable_length_filtering",
                        "--trim_front1=5",
                        "--trim_front2=5",
                    ].join(' '),
                    output_prefix: commonPrefix(r.reads.first().baseName, r.reads.last().baseName)
                )
            }
    )

    fastq = ch_fastq
        .filter { _r -> !bool_trim_fastq }
        .mix(trimmed.map { t -> record(id: t.id, meta: t.meta, reads: t.reads.toList().toSorted()) })

    // Convert to cram
    cram = SAMTOOLS_IMPORT(
        fastq.map { r ->
            def name = r.reads.first().simpleName
            record(
                id: r.id,
                meta: r.meta,
                single_end: r.meta.single_end,
                reads: r.reads,
                samtools_args: [
                    '--output-fmt cram',
                    "-r ID:${name}",             // Read group identifier
                    "-r CN:${hic_type}",         // ? Center Name | instead describes library type
                    "-r PU:${name}",             // ? Platform Unit - Identifier for the sequencer unit | instead describes sample
                    "-r SM:${r.meta.id}",        // Sample Name
                ].join(' '),
                output_prefix: name.endsWith('_1') ? name.substring(0, name.length() - 2) : name
            )
        }
    )

    // Index Cram
    crai = SAMTOOLS_INDEX(
        cram.map { c -> record(id: c.id, meta: c.meta, input: c.cram) }
    )

    // The converted read pairs of each sample, in the order of the input
    hic = fastq
        .map { r -> tuple(r.id, r) }
        .groupBy()
        .map { sample_id, results ->
            def sorted = results.toList().toSorted { a, b -> a.meta.pair_id <=> b.meta.pair_id }
            record(id: sample_id, hic_fastq: sorted.collect { r -> r.reads.toList() }.toList())
        }
        .mix(samples.filter { s -> !s.hic }.map { s -> record(id: s.id, hic_fastq: null as List<List<Path>>) })

    emit:
    fastq: Channel<Reads>  = fastq
    hic  = hic
    cram : Channel<Record> = cram
    crai : Channel<Record> = crai
}
