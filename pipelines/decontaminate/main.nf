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
include { Assembly; KmerParams; BuscoParams; assemblyWith } from '../../types/records'

record FcsParams {
    database    : Path?   // Local copy of the FCS-GX database. Fetched with the manifest when not set.
    manifest    : String? // URL of the FCS-GX database manifest
    ramdisk_path: String? // Copy the database here before running FCS-GX
}

params {
    assemblies: Channel<Assembly>  // Assemblies to screen
    provided  : Channel<Assembly>? // Assemblies of this stage that were supplied. Evaluated, not screened.

    fcs       : FcsParams   = record(database: null, manifest: null, ramdisk_path: null)

    // Assembly QC
    busco     : BuscoParams = record(lineages: null, lineages_db_path: null, max_lineages: null)
    fastk     : KmerParams  = record(kmer_size: null)
    meryl     : KmerParams  = record(kmer_size: null)
    reference : Path?
    use_phased: Boolean     = false
}

workflow {

    main:
    def fcs_database = params.fcs.database
    def fcs_manifest = params.fcs.manifest ?: "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest"
    def ramdisk_path = params.fcs.ramdisk_path != null ? params.fcs.ramdisk_path : '/dev/shm/gxdb'

    // Fetch the database only when there are assemblies to screen and no local copy
    fetched = FCSGX_FETCHDB(
        params.assemblies
            .map { _a -> record(id: 'fcsgx_database', manifest: fcs_manifest, output_prefix: null) }
            .unique()
            .filter { _x -> !fcs_database }
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
        .map { a ->
            assemblyWith(
                a,
                [
                    pri_fasta: a.cleaned[0],
                    alt_fasta: a.cleaned.size() > 1 ? a.cleaned[1] : null,
                    pri_gfa: null,
                    alt_gfa: null
                ]
            )
        }
    assemblies = cleaned_assemblies.mix(params.provided ?: channel.empty())

    qc = ASSEMBLY_QC_DECONTAMINATED(
        record(
            assemblies: assemblies,
            busco: params.busco,
            fastk: params.fastk,
            meryl: params.meryl,
            reference: params.reference,
            use_phased: params.use_phased
        )
    )

    publish:
    assemblies    = assemblies
    fcsgx_reports = screened
    fcsgx_clean   = cleaned
    merquryfk     = qc.merquryfk
    merqury       = qc.merqury
    busco         = qc.busco
    gfastats      = qc.gfastats
    quast         = qc.quast
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
