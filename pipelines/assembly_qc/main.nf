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
include { Assembly; Reads; FastkDatabase; MerylDatabase; KmerParams; BuscoParams } from '../../types/records'

// The sample, assembly and reads given in a params file
record SampleInput {
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
    sample    : SampleInput?
    assembly  : List<AssemblyInput>?
    hifi      : List<ReadsInput>? // HiFi reads to build the k-mer databases

    // Options
    busco     : BuscoParams = record(lineages: null, lineages_db_path: null, max_lineages: null)
    fastk     : KmerParams  = record(kmer_size: null)
    meryl     : KmerParams  = record(kmer_size: null)
    use_phased: Boolean     = false // Evaluate the primary and alternate assembly separately
}

// Assemblies from the params of a run on its own
def assembliesFromParams(sample: SampleInput, assemblies: List<AssemblyInput>, hifi: List<ReadsInput>, busco_lineages: String) -> List<Assembly> {
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
            pri_gfa        : null as Path,
            hifi           : hifi.collect { h -> h.reads }.toList(),
            hifi_merged    : null as Path,
            hic            : null as List<List<Path>>,
            rnaseq         : null as List<List<Path>>,
            fastk_hifi_hist: null as Path,
            fastk_hifi_ktab: null as List<Path>,
            fastk_hifi_prof: null as List<Path>,
            meryl_hifi_db  : null as Path
        )
    }.toList()
}

// The HiFi reads of an assembly's sample, one record per file
def hifiReads(a: Record) -> List<Reads> {
    // Only what belongs to the sample: the settings of a build (for example its hifiasm options) would make
    // the same reads count as different ones, and the k-mer databases would be built once per build
    def sample_meta = record(id: a.meta.id, sample: a.meta.sample, settings: record(busco: a.meta.settings.busco), single_end: true)
    return ( a.hifi ?: [] ).collect { f -> record(id: a.meta.id, meta: sample_meta, reads: [ f ]) }.toList()
}

// Fails when the assembly has neither HiFi k-mer databases nor HiFi reads to build them from
def requireHifiSource(a: Record) -> Boolean {
    if ( a.fastk_hifi_hist == null && !a.hifi ) {
        error "Assembly ${a.id} has no HiFi k-mer database and no HiFi reads (hifi) to build one."
    }
    return true
}

def withDatabases(a: Assembly, fastk: FastkDatabase, meryl: MerylDatabase) {
    return a + record(fastk_hifi_hist: fastk.hist, fastk_hifi_ktab: fastk.ktab, fastk_hifi_prof: fastk.prof, meryl_hifi_db: meryl.meryl_db)
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
    // The HiFi reads are required for the k-mer databases of an assembly that did not come with them
    inputs = from_params
        .mix( params.assemblies ?: channel.empty() )
        .filter { a -> requireHifiSource(a) }

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
        params.busco.max_lineages ?: 2,
        params.busco.lineages_db_path
    )

    publish:
    merquryfk = evaluation.merquryfk
    merqury   = evaluation.merqury
    busco     = evaluation.busco
    gfastats  = evaluation.gfastats
    logs      = evaluation.logs
}

output {
    merquryfk: Channel<Record> {
        path { r -> "merquryfk/${r.meta.assembly.build}" }
    }
    merqury: Channel<Record> {
        path { r -> "merqury/${r.meta.assembly.build}" }
    }
    busco: Channel<Record> {
        path { r -> "busco/${r.meta.assembly.build}" }
    }
    gfastats: Channel<Record> {
        path { r -> "gfastats/${r.meta.assembly.build}" }
    }
    logs: Channel<Path> {
        path 'logs'
    }
}
