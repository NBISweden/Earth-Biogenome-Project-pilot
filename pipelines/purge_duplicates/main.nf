#!/usr/bin/env nextflow

/*
 * Duplicate purging with purge_dups, followed by assembly QC.
 * https://github.com/dfguan/purge_dups
 *
 * Takes one wide assembly record per assembly build (see types/records.nf). The HiFi reads
 * are read from `hifi_merged`. Returns the purged assemblies with the same fields.
 */

nextflow.enable.types = true

include { ITERATIVE_PURGEDUPS            } from '../../modules/local/purgedups/main'
include { workflow as ASSEMBLY_QC_PURGED } from '../assembly_qc'
include { Assembly; QcParams; defaultQcParams; withFastas } from '../../types/records'

params {
    assemblies: Channel<Assembly>  // Assemblies to purge
    provided  : Channel<Assembly>? // Assemblies of this stage that were supplied. Evaluated, not purged.

    // Assembly QC
    qc        : QcParams    = defaultQcParams() // Options of the assembly QC. use_phased: also process and evaluate the alternate assembly
}

// The assembly files of a purge, primary first
def purgedFiles(purged: List<Path>) -> List<Path> {
    return purged.toSorted { a, b -> a.name <=> b.name }
}

workflow {

    main:
    purged = ITERATIVE_PURGEDUPS(
        params.assemblies.map { a ->
            record(
                id: a.id,
                meta: a.meta,
                reads: [ a.hifi_merged ],
                assemblies: a.alt_fasta ? [ a.pri_fasta, a.alt_fasta ] : [ a.pri_fasta ],
                minimap2_reads_args: "-x map-hifi",
                minimap2_self_args: "-x asm5 -DP",
                pbcstat_args: "",
                calcuts_args: "",
                hist_plot_args: "",
                split_fa_args: "",
                purge_dups_args: "-2",
                get_seqs_args: "-e",
                seqkit_args: "-w 100",
                output_prefix: "${a.meta.id}_${a.meta.assembly.build}"
            )
        }
    )

    // The wide record, with the purged assembly files
    purged_assemblies = params.assemblies
        .join(purged.map { p -> record(id: p.id, purged: purgedFiles(p.purged.toList())) }, by: 'id')
        .map { a -> withFastas(a, a.purged) }
    assemblies = purged_assemblies.mix(params.provided ?: channel.empty())

    qc = ASSEMBLY_QC_PURGED(
        record(assemblies: assemblies) + params.qc
    )

    publish:
    assemblies = assemblies
    purge_dups = purged
    merquryfk  = qc.merquryfk
    merqury    = qc.merqury
    busco      = qc.busco
    gfastats   = qc.gfastats
    logs       = purged.flatMap { p -> p.logs.toList() }.mix(qc.logs)
}

output {
    assemblies: Channel<Assembly> {
        path { a ->
            a.pri_fasta >> ("assemblies/${a.id}/")
            a.alt_fasta >> ("assemblies/${a.id}/")
        }
    }
    purge_dups: Channel<Record> {
        path { r -> "purge_dups/${r.meta.assembly.build}" }
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
