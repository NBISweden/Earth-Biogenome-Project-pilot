#!/usr/bin/env nextflow

/*
 * Assembly QC: k-mer consistency (MerquryFK, Merqury), gene space completeness (BUSCO) and
 * contiguity (gfastats, QUAST) of one or more assemblies of a sample.
 *
 * Runs on its own, or is included by the stage pipelines and the EBP assembly workflow:
 *
 *   include { workflow as ASSEMBLY_QC_PURGED } from '../assembly_qc'
 *   qc = ASSEMBLY_QC_PURGED( record(assemblies: purged, busco: params.busco) )
 *
 * Usage:
 *   nextflow run pipelines/assembly_qc -params-file qc_parameters.yml -profile singularity
 *
 * See README.md for the parameter file.
 */

nextflow.enable.types = true

include { BUILD_FASTK_DATABASE as BUILD_FASTK_HIFI_DATABASE } from '../../subworkflows/local/build_fastk_database/main'
include { BUILD_MERYL_DATABASE as BUILD_MERYL_HIFI_DATABASE } from '../../subworkflows/local/build_meryl_database/main'
include { EVALUATE_ASSEMBLY                                 } from '../../subworkflows/local/evaluate_assembly/main'
include { COMPARE_ASSEMBLIES                                } from '../../subworkflows/local/compare_assemblies/main'
include { Assembly; Reads; FastkDatabase; MerylDatabase; KmerParams; BuscoParams; assemblyOf; assemblyWith } from '../../types/records'

// The sample, assembly and reads given in a params file
record Sample {
    name       : String
    genome_size: Number  // Haploid genome size in bases. Used to compute N*-statistics.
    ploidy     : Integer?
}

record AssemblyInput {
    assembler: String
    stage    : String // raw, decontaminated, purged, polished, scaffolded or curated
    id       : String
    pri_fasta: Path
    alt_fasta: Path?
}

record ReadsInput {
    reads: Path
}

params {
    // Assemblies of a calling pipeline. Wide records, with the HiFi k-mer databases (fastk_hifi_*, meryl_hifi_db)
    // or the HiFi reads (hifi) to build them from.
    assemblies: Channel<Assembly>?

    // Sample and data of a run on its own. Give these in a -params-file.
    sample    : Sample?
    assembly  : List<AssemblyInput>?
    hifi      : List<ReadsInput>? // HiFi reads to build the k-mer databases

    // Options
    busco     : BuscoParams = record(lineages: null, lineages_db_path: null, max_lineages: null)
    fastk     : KmerParams  = record(kmer_size: null)
    meryl     : KmerParams  = record(kmer_size: null)
    reference : Path?       // Reference genome to compare to with QUAST
    use_phased: Boolean     = false // Evaluate the primary and alternate assembly separately
}

// Assemblies from the params of a run on its own
def assembliesFromParams(sample: Sample, assemblies: List<AssemblyInput>, hifi: List<ReadsInput>, busco_lineages: String) -> List<Assembly> {
    def sample_id = sample.name.replace(' ', '_')
    def sample_meta = record(
        id      : sample_id,
        sample  : [ name: sample.name, genome_size: sample.genome_size, ploidy: sample.ploidy ],
        settings: [ busco: [ lineages: busco_lineages ] ]
    )
    return assemblies.collect { a ->
        def build = "${a.assembler}-${a.stage}-${a.id}"
        record(
            id             : build,
            meta           : sample_meta + record(assembly: [ assembler: a.assembler, stage: a.stage, id: a.id, build: build ]),
            pri_fasta      : a.pri_fasta,
            alt_fasta      : a.alt_fasta,
            pri_gfa        : null,
            alt_gfa        : null,
            hifi           : hifi.collect { h -> h.reads }.toList(),
            hifi_merged    : null,
            hic            : null,
            fastk_hifi_hist: null,
            fastk_hifi_ktab: null,
            fastk_hifi_prof: null,
            meryl_hifi_db  : null,
            kmer_cov       : null
        )
    }.toList()
}

// The HiFi reads of an assembly's sample, one record per file
def hifiReads(a: Record) -> List<Reads> {
    def sample_meta = record(id: a.meta.id, sample: a.meta.sample, settings: a.meta.settings, single_end: true)
    return ( a.hifi ?: [] ).collect { f -> record(id: a.meta.id, meta: sample_meta, reads: [ f ]) }.toList()
}

def withDatabases(a: Record, fastk: FastkDatabase, meryl: MerylDatabase) -> Assembly {
    return assemblyWith(a, [ fastk_hifi_hist: fastk.hist, fastk_hifi_ktab: fastk.ktab, fastk_hifi_prof: fastk.prof, meryl_hifi_db: meryl.meryl_db ])
}

workflow {

    main:
    if ( !params.assemblies && !( params.sample && params.assembly ) ) {
        error "Give assemblies from a calling pipeline, or the sample and assembly params."
    }
    if ( !params.busco.lineages && !params.assemblies ) {
        error "Missing parameter busco.lineages. This pipeline does not look up the BUSCO lineages of the sample."
    }

    from_params = channel.fromList(
        params.sample && params.assembly
            ? assembliesFromParams(params.sample, params.assembly, params.hifi ?: [], params.busco.lineages ?: '')
            : []
    )
    inputs = from_params.mix(( params.assemblies ?: channel.empty() ).map { a -> assemblyOf(a) })

    // K-mer databases for the assemblies that did not come with them
    with_databases = inputs.filter { a -> a.fastk_hifi_hist != null && a.meryl_hifi_db != null }
    without_databases = inputs.filter { a -> a.fastk_hifi_hist == null || a.meryl_hifi_db == null }
    reads = without_databases.flatMap { a -> hifiReads(a) }.unique()
    fastk = BUILD_FASTK_HIFI_DATABASE( reads, 'hifi', params.fastk.kmer_size ?: 31 )
    meryl = BUILD_MERYL_HIFI_DATABASE( reads, 'hifi', params.meryl.kmer_size ?: 31 )
    built = without_databases
        .map { a -> record(id: a.id, sample_id: a.meta.id, assembly: a) }
        .join(fastk.database.map { d -> record(sample_id: d.id, fastk: d) }, by: 'sample_id')
        .join(meryl.uniondb.map { m -> record(sample_id: m.id, meryl: m) }, by: 'sample_id')
        .map { x -> withDatabases(x.assembly, x.fastk, x.meryl) }
    assemblies = with_databases.mix(built)

    evaluation = EVALUATE_ASSEMBLY(
        assemblies,
        params.use_phased,
        params.busco.max_lineages ?: 3,
        params.busco.lineages_db_path
    )
    comparison = COMPARE_ASSEMBLIES( assemblies, params.reference )

    publish:
    merquryfk = evaluation.merquryfk
    merqury   = evaluation.merqury
    busco     = evaluation.busco
    gfastats  = evaluation.gfastats
    quast     = comparison.quast
    logs      = evaluation.logs.mix(comparison.logs)
}

output {
    merquryfk: Channel<Record> {
        path { r -> "merquryfk/${r.meta.assembly.build}" }
    }
    merqury: Channel<Record> {
        path { r -> "merqury/${r.meta.assembly.build}" }
    }
    busco: Channel<Record> {
        path { r ->
            r.short_summaries_txt >> ("busco/${r.meta.assembly.build}/")
            r.short_summaries_json >> ("busco/${r.meta.assembly.build}/")
        }
    }
    gfastats: Channel<Record> {
        path { r -> "gfastats/${r.meta.assembly.build}" }
    }
    quast: Channel<Record> {
        path 'quast'
    }
    logs: Channel<Path> {
        path 'logs'
    }
}
