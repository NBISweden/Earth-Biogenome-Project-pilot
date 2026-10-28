nextflow.enable.types = true

include { SEQKIT_STATS as SEQKIT_HIFI_STATS } from "../../../modules/nf-core/seqkit/stats/main"
include { SEQKIT_STATS as SEQKIT_HIC_STATS  } from "../../../modules/nf-core/seqkit/stats/main"
include { FASTQC                            } from "../../../modules/nf-core/fastqc/main"
include { GENOME_PROPERTIES                 } from "./genome_properties"
include { COMPARE_LIBRARIES                 } from "./compare_libraries"

include { Reads; FastkDatabase } from "../../../types/records"

// File name without the .fastq/.fasta and optional .gz extensions
def readsBaseName(reads: List<Path>) -> String {
    def first = reads[0]
    return first.name.replaceAll(/\.(fastq|fq|fasta|fa)(\.gz)?$/, '')
}

workflow INSPECT_DATA {
    take:
    samples  : Channel<Record> // With the reads (hifi, hic) and the FastK databases (fastk_hifi_*, fastk_hic_*)
    kmer_size: Integer         // GenomeScope k-mer length

    main:
    // Reads and databases of each sample
    hifi_reads = samples.flatMap { s ->
        ( s.hifi ?: [] ).collect { f -> record(id: s.id, meta: s.meta + record(single_end: true), reads: [ f ]) }.toList()
    }
    hic_reads = samples.flatMap { s ->
        ( s.hic ?: [] ).withIndex().collect { pair, index ->
            record(id: s.id, meta: s.meta + record(single_end: false, pair_id: index), reads: pair)
        }.toList()
    }
    hifi_database = samples
        .filter { s -> s.fastk_hifi_hist != null }
        .map { s -> record(id: s.id, meta: s.meta, hist: s.fastk_hifi_hist, ktab: s.fastk_hifi_ktab, prof: s.fastk_hifi_prof) }
    hic_database = samples
        .filter { s -> s.fastk_hic_hist != null }
        .map { s -> record(id: s.id, meta: s.meta, hist: s.fastk_hic_hist, ktab: s.fastk_hic_ktab, prof: s.fastk_hic_prof) }

    // Quantify data
    hifi_stats = SEQKIT_HIFI_STATS(
        hifi_reads.map { r -> record(id: r.id, meta: r.meta, reads: r.reads, output_prefix: "${readsBaseName(r.reads)}_hifi") }
    )
    hic_stats = SEQKIT_HIC_STATS(
        hic_reads.map { r -> record(id: r.id, meta: r.meta, reads: r.reads, output_prefix: "${readsBaseName(r.reads)}_hic") }
    )
    fastqc = FASTQC(
        hic_reads.map { r -> record(id: r.id, meta: r.meta, reads: r.reads, output_prefix: readsBaseName(r.reads)) }
    )

    // Generate K-mer histogram
    genome_properties = GENOME_PROPERTIES( hifi_database, kmer_size )
    compare_libraries = COMPARE_LIBRARIES(
        hifi_database
            .map { d -> record(id: d.id, meta: d.meta, fastk1_hist: d.hist, fastk1_ktab: d.ktab) }
            .join(hic_database.map { d -> record(id: d.id, fastk2_hist: d.hist, fastk2_ktab: d.ktab) }, by: 'id')
    )

    logs = genome_properties.logs
        .mix(compare_libraries.logs)
        .mix(hifi_stats.map { s -> s.stats })
        .mix(hic_stats.map { s -> s.stats })
        .mix(fastqc.flatMap { f -> f.zip.toList() })

    emit:
    logs         : Channel<Path>   = logs
    seqkit_hifi  : Channel<Record> = hifi_stats
    seqkit_hic   : Channel<Record> = hic_stats
    fastqc       : Channel<Record> = fastqc
    genescopefk  : Channel<Record> = genome_properties.genescopefk
    smudgeplot   : Channel<Record> = genome_properties.smudgeplot
    katgc        : Channel<Record> = genome_properties.katgc
    katcomp      : Channel<Record> = compare_libraries.katcomp
}
