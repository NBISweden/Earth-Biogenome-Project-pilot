#!/usr/bin/env nextflow

/*
 * Contamination screening with FCS-GX, followed by assembly QC.
 * https://github.com/ncbi/fcs
 *
 * Takes one wide assembly record per assembly build (see types/records.nf). The assemblies carry
 * the stage `decontaminated` in meta.assembly, which callers set. Returns the cleaned assemblies
 * with the same fields. Needs the taxid of the sample in meta.sample.tax_id.
 */

nextflow.enable.types = true

include { FCSGX_FETCHDB                          } from '../../modules/local/fcsgx/fetchdb'
include { FCSGX_RUNGX                            } from '../../modules/local/fcsgx/rungx'
include { FCSGX_CLEAN                            } from '../../modules/local/fcsgx/clean'
include { workflow as ASSEMBLY_QC_DECONTAMINATED } from '../assembly_qc'
include { Assembly; QcParams; defaultQcParams; FcsParams; withFastas } from '../../types/records'

params {
    assemblies: Channel<Assembly>  // Assemblies to screen
    provided  : Channel<Assembly>? // Assemblies of this stage that were supplied. Evaluated, not screened.

    fcs       : FcsParams   = record(database: null, manifest: null, ramdisk_path: null)
    fetch_database: Boolean = true // Fetch the FCS-GX database at the start unless fcs.database is given

    // Assembly QC
    qc        : QcParams    = defaultQcParams() // Options of the assembly QC. use_phased: also process and evaluate the alternate assembly
}

workflow {

    main:
    def fcs_database = params.fcs.database
    def fcs_manifest = params.fcs.manifest ?: "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest"
    def ramdisk_path = params.fcs.ramdisk_path != null ? params.fcs.ramdisk_path : '/dev/shm/gxdb'

    // Fetch the database at the start, while the earlier stages run, unless there is a local copy.
    // The fetch does not wait for the assemblies: it takes up to two days on a first run.
    fetched = FCSGX_FETCHDB(
        channel.of( record(id: 'fcsgx_database', manifest: fcs_manifest, output_prefix: null) )
            .filter { _x -> !fcs_database && params.fetch_database }
    )
    // The local database, or the fetched one
    database = channel.of(fcs_database)
        .filter { db -> db != null }
        .mix(fetched.map { f -> f.database })
        .collect()
        .map { dbs -> dbs.toList()[0] }

    // One screen per haplotype
    to_screen = params.assemblies.flatMap { a ->
        def taxid = "${a.meta.sample.tax_id}"
        def haplotypes = a.alt_fasta ? [ a.pri_fasta, a.alt_fasta ] : [ a.pri_fasta ]
        haplotypes.withIndex().collect { fasta, index ->
            record(id: "${a.id}-hap${index + 1}", meta: a.meta + record(haplotype: index + 1), taxid: taxid, fasta: fasta)
        }.toList()
    }

    screened = FCSGX_RUNGX(
        to_screen
            .combine(gxdb: database)
            .map { s -> record(id: s.id, meta: s.meta, taxid: s.taxid, assembly: s.fasta, gxdb: s.gxdb, ramdisk_path: ramdisk_path) }
    )

    cleaned = FCSGX_CLEAN(
        to_screen
            .join(screened.map { r -> record(id: r.id, action_report: r.fcs_gx_report) }, by: 'id')
            .map { s -> record(id: s.id, meta: s.meta, assembly: s.fasta, action_report: s.action_report) }
    )

    // Recombine the haplotypes of each assembly, primary first, into the wide record
    cleaned_files = cleaned
        .map { c -> tuple(c.meta.assembly.build, c) }
        .groupBy()
        .map { build, results ->
            def sorted = results.toSorted { a, b -> a.clean_fasta.name <=> b.clean_fasta.name }
            record(
                id: build,
                cleaned: sorted.collect { c -> c.clean_fasta }.toList()
            )
        }
    cleaned_assemblies = params.assemblies
        .join(cleaned_files, by: 'id')
        .map { a -> withFastas(a, a.cleaned) }
    assemblies = cleaned_assemblies.mix(params.provided ?: channel.empty())

    qc = ASSEMBLY_QC_DECONTAMINATED(
        record(assemblies: assemblies) + params.qc
    )

    publish:
    assemblies    = assemblies
    fcsgx_reports = screened
    fcsgx_clean   = cleaned
    merquryfk     = qc.merquryfk
    merqury       = qc.merqury
    busco         = qc.busco
    gfastats      = qc.gfastats
    logs          = screened.map { r -> r.fcs_gx_report }.mix(qc.logs)
}

output {
    assemblies: Channel<Assembly> {
        path { a ->
            a.pri_fasta >> ("assemblies/${a.id}/")
            a.alt_fasta >> ("assemblies/${a.id}/")
        }
    }
    fcsgx_reports: Channel<Record> {
        path { r -> "fcsgx/${r.meta.assembly.build}" }
    }
    fcsgx_clean: Channel<Record> {
        path { r -> "fcsgx/${r.meta.assembly.build}" }
    }
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
