nextflow.enable.types = true

include { HIFIASM         } from "../../../modules/nf-core/hifiasm/main"
include { GFATOOLS_GFA2FA } from "../../../modules/local/gfatools/gfa2fa"

record Reads {
    id   : String
    meta : Record
    reads: List<Path>
}

record Assembly {
    id       : String  // meta.assembly.build
    meta     : Record  // sample meta data, including meta.assembly
    pri_fasta: Path
    alt_fasta: Path?
    pri_gfa  : Path?
    alt_gfa  : Path?
}

record HifiasmInput {
    id          : String
    meta        : Record
    reads       : List<Path>
    hifiasm_args: String
    hic_read1   : List<Path>
    hic_read2   : List<Path>
}

def hifiasmInput(id: String, meta: Record, reads: List<Path>, hifiasm_args: String, hic_read1: List<Path>, hic_read2: List<Path>) -> HifiasmInput {
    return record(id: id, meta: meta, reads: reads, hifiasm_args: hifiasm_args, hic_read1: hic_read1, hic_read2: hic_read2)
}

// Add the hifiasm configuration and the assembly build to the sample meta data
def assemblyMeta(meta: Record, key: String, args: String) -> Record {
    return meta + record(
        settings: meta.settings + [ hifiasm: [ id: key, args: args ] ],
        assembly: [ assembler: 'hifiasm', stage: 'raw', id: key, build: "hifiasm-raw-${key}" ]
    )
}

workflow ASSEMBLE_HIFI {
    take:
    hifi_reads      : Channel<Reads>
    hic_reads       : Channel<Reads>
    hifiasm_configs : Map<String,String> // [ build_id: "hifiasm args" ]. Empty for a single default build.
    hifiasm_with_hic: Boolean
    use_phased      : Boolean

    main:
    // One assembly per hifiasm configuration
    builds = hifi_reads.flatMap { r ->
        def configs = hifiasm_configs ?: [ 'default': '' ]
        configs.keySet().collect { key ->
            record(
                id: "hifiasm-raw-${key}",
                sample_id: r.id,
                meta: assemblyMeta(r.meta, key, configs[key]),
                reads: r.reads,
                hifiasm_args: configs[key]
            )
        }.toList()
    }

    // Hi-C reads of a sample, in pair order
    hic_by_sample = hic_reads
        .map { r -> tuple(r.id, r) }
        .groupBy()
        .map { sample_id, pairs ->
            def sorted = pairs.toSorted { a, b -> a.meta.pair_id <=> b.meta.pair_id }
            record(
                sample_id: sample_id,
                hic_read1: sorted.collect { p -> p.reads.first() }.toList(),
                hic_read2: sorted.collect { p -> p.reads.last() }.toList()
            )
        }

    inputs = builds
        .filter { _b -> !hifiasm_with_hic }
        .map { b -> hifiasmInput(b.id, b.meta, b.reads, b.hifiasm_args, [], []) }
        .mix(
            builds
                .filter { _b -> hifiasm_with_hic }
                .join(hic_by_sample, by: 'sample_id')
                .map { b -> hifiasmInput(b.id, b.meta, b.reads, b.hifiasm_args, b.hic_read1, b.hic_read2) }
        )

    hifiasm = HIFIASM(
        inputs.map { i ->
            record(
                id: i.id,
                meta: i.meta,
                long_reads: i.reads,
                ul_reads: [],
                paternal_kmer_dump: null,
                maternal_kmer_dump: null,
                hic_read1: i.hic_read1,
                hic_read2: i.hic_read2,
                bin_files: [],
                hifiasm_args: i.hifiasm_args
            )
        }
    )

    // Convert the contig graphs to fasta
    fastas = GFATOOLS_GFA2FA(
        hifiasm.flatMap { h ->
            def gfas = use_phased ? [ h.hap1_contigs, h.hap2_contigs ] : [ h.primary_contigs ]
            gfas.findAll { g -> g != null }.collect { g -> record(id: h.id, meta: h.meta, gfa: g, output_prefix: g.baseName) }.toList()
        }
    )

    gfas = hifiasm.map { h ->
        record(
            id: h.id,
            pri_gfa: use_phased ? h.hap1_contigs : h.primary_contigs,
            alt_gfa: use_phased ? h.hap2_contigs : null
        )
    }

    assemblies = fastas
        .map { f -> tuple(f.id, f) }
        .groupBy()
        .map { id, results ->
            def sorted = results.toSorted { a, b -> a.fasta.name <=> b.fasta.name }
            record(
                id: id,
                meta: sorted[0].meta,
                pri_fasta: sorted[0].fasta,
                alt_fasta: sorted.size() > 1 ? sorted[1].fasta : null
            )
        }
        .join(gfas, by: 'id')

    emit:
    assemblies: Channel<Assembly> = assemblies
    logs      : Channel<Path>     = hifiasm.map { h -> h.log }
}
