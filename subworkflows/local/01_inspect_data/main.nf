nextflow.enable.types = true

include { SEQKIT_STATS as SEQKIT_HIFI_STATS } from "../../../modules/nf-core/seqkit/stats/main"
include { SEQKIT_STATS as SEQKIT_HIC_STATS  } from "../../../modules/nf-core/seqkit/stats/main"
include { FASTQC                            } from "../../../modules/nf-core/fastqc/main"
include { GENOME_PROPERTIES                 } from "./genome_properties"
include { COMPARE_LIBRARIES                 } from "./compare_libraries"

record Reads {
    id   : String
    meta : Record
    reads: List<Path>
}

record FastkDatabase {
    id  : String
    meta: Record
    hist: Path
    ktab: List<Path>
    prof: List<Path>
}

// File name without the .fastq/.fasta and optional .gz extensions
def readsBaseName(reads: List<Path>) -> String {
    def first = reads[0]
    return first.name.replaceAll(/\.(fastq|fq|fasta|fa)(\.gz)?$/, '')
}

workflow INSPECT_DATA {
    take:
    hifi_reads    : Channel<Reads>
    hic_reads     : Channel<Reads>
    hifi_database : Channel<FastkDatabase>
    hic_database  : Channel<FastkDatabase>
    kmer_size     : Integer // GenomeScope k-mer length

    main:
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
    hifi_kmer_cov: Channel<Record> = genome_properties.kmer_cov
    logs         : Channel<Path>   = logs
    seqkit_hifi  : Channel<Record> = hifi_stats
    seqkit_hic   : Channel<Record> = hic_stats
    fastqc       : Channel<Record> = fastqc
    genescopefk  : Channel<Record> = genome_properties.genescopefk
    smudgeplot   : Channel<Record> = genome_properties.smudgeplot
    katgc        : Channel<Record> = genome_properties.katgc
    katcomp      : Channel<Record> = compare_libraries.katcomp
}
