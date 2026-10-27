#!/usr/bin/env nextflow

/*
 * Assembly QC: k-mer consistency (MerquryFK, Merqury), gene space completeness (BUSCO) and
 * contiguity (gfastats, QUAST) of one or more assemblies of a sample.
 *
 * Runs on its own, or as the evaluation step of the EBP assembly workflow.
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

record BuscoParams {
    lineages        : String? // Comma separated BUSCO lineages
    max_lineages    : Integer?
    lineages_db_path: Path?   // Directory with downloaded BUSCO lineages. Downloads when not set.
}

record KmerParams {
    kmer_size: Integer?
}

params {
    // Sample and data. Give these in a -params-file.
    sample    : Sample
    assembly  : List<AssemblyInput>
    hifi      : List<ReadsInput> // HiFi reads to build the k-mer databases

    // Options
    busco     : BuscoParams = record(lineages: null, max_lineages: null, lineages_db_path: null)
    fastk     : KmerParams  = record(kmer_size: null)
    meryl     : KmerParams  = record(kmer_size: null)
    reference : Path?       // Reference genome to compare to with QUAST
    use_phased: Boolean     = false // Evaluate the primary and alternate assembly separately
}

workflow {

    main:
    if ( !params.busco.lineages ) {
        error "Missing parameter busco.lineages. This pipeline does not look up the BUSCO lineages of the sample."
    }

    def sample_id = params.sample.name.replace(' ', '_')
    def sample_meta = record(
        id      : sample_id,
        sample  : [ name: params.sample.name, genome_size: params.sample.genome_size, ploidy: params.sample.ploidy ],
        settings: [ busco: [ lineages: params.busco.lineages ] ]
    )

    // HiFi reads, one record per file
    reads = channel.fromList(
        params.hifi.collect { h -> record(id: sample_id, meta: sample_meta, reads: [ h.reads ]) }.toList()
    )

    // Assemblies, with the assembly build in the meta data
    assemblies = channel.fromList(
        params.assembly.collect { a ->
            def build = "${a.assembler}-${a.stage}-${a.id}"
            record(
                id       : build,
                meta     : sample_meta + record(assembly: [ assembler: a.assembler, stage: a.stage, id: a.id, build: build ]),
                pri_fasta: a.pri_fasta,
                alt_fasta: a.alt_fasta,
                pri_gfa  : null,
                alt_gfa  : null
            )
        }.toList()
    )

    // K-mer databases of the reads
    fastk = BUILD_FASTK_HIFI_DATABASE( reads, 'hifi', params.fastk.kmer_size ?: 31 )
    meryl = BUILD_MERYL_HIFI_DATABASE( reads, 'hifi', params.meryl.kmer_size ?: 31 )

    evaluation = EVALUATE_ASSEMBLY(
        assemblies,
        fastk.database,
        meryl.uniondb,
        params.use_phased,
        params.busco.max_lineages ?: 3,
        params.busco.lineages_db_path
    )
    comparison = COMPARE_ASSEMBLIES( assemblies, params.reference )

    publish:
    merquryfk = evaluation.merquryfk
    busco     = evaluation.busco
    gfastats  = evaluation.gfastats
    quast     = comparison.quast
}

output {
    merquryfk: Channel<Record> {
        path { r -> "merquryfk/${r.meta.assembly.build}" }
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
}
