nextflow.enable.types = true

include { HIFIASM         } from "../../modules/nf-core/hifiasm/main"
include { GFATOOLS_GFA2FA } from "../../modules/local/gfatools/gfa2fa"
include { Assembly; Sample } from "../../types/records"

// Add the hifiasm configuration and the assembly build to the sample meta data
def assemblyMeta(meta: Record, key: String, args: String) -> Record {
    return meta + record(
        settings: meta.settings + [ hifiasm: [ id: key, args: args ] ],
        assembly: [ assembler: 'hifiasm', stage: 'raw', id: key, build: "hifiasm-raw-${key}" ]
    )
}

workflow ASSEMBLE_HIFI {
    take:
    samples         : Channel<Sample>
    hifiasm_configs : Map<String,String> // [ build_id: "hifiasm args" ]. Empty for a single default build.
    hifiasm_with_hic: Boolean
    use_phased      : Boolean

    main:
    // One assembly per hifiasm configuration. The Hi-C pairs are used when requested.
    builds = samples.flatMap { s ->
        def configs = hifiasm_configs ?: [ 'default': '' ]
        def pairs = hifiasm_with_hic ? ( s.hic ?: [] ) : []
        configs.keySet().collect { key ->
            record(
                id: "hifiasm-raw-${key}",
                sample: s,
                meta: assemblyMeta(s.meta, key, configs[key]),
                reads: [ s.hifi_merged ],
                hifiasm_args: configs[key],
                hic_read1: pairs.collect { p -> p.first() }.toList(),
                hic_read2: pairs.collect { p -> p.last() }.toList()
            )
        }.toList()
    }

    hifiasm = HIFIASM(
        builds.map { b ->
            record(
                id: b.id,
                meta: b.meta,
                long_reads: b.reads,
                ul_reads: [],
                hic_read1: b.hic_read1,
                hic_read2: b.hic_read2,
                bin_files: [],
                hifiasm_args: b.hifiasm_args
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
            pri_gfa: use_phased ? h.hap1_contigs : h.primary_contigs
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
        .join(builds.map { b -> record(id: b.id, sample: b.sample) }, by: 'id')
        // The sample record, with the build in place of the sample id and meta
        .map { x -> x.sample + record(id: x.meta.assembly.build, meta: x.meta, pri_fasta: x.pri_fasta, alt_fasta: x.alt_fasta, pri_gfa: x.pri_gfa) }

    emit:
    assemblies: Channel<Assembly> = assemblies
    logs      : Channel<Path>     = hifiasm.map { h -> h.log }
    hifiasm   : Channel<Record>   = hifiasm
    gfa2fa    : Channel<Record>   = fastas
}
