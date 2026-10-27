nextflow.enable.types = true

include { FCSGX_FETCHDB } from "../../../modules/local/fcsgx/fetchdb"
include { FCSGX_RUNGX   } from "../../../modules/local/fcsgx/rungx"
include { FCSGX_CLEAN   } from "../../../modules/local/fcsgx/clean"

record Assembly {
    id       : String  // meta.assembly.build
    meta     : Record  // sample meta data, including meta.assembly
    pri_fasta: Path
    alt_fasta: Path?
    pri_gfa  : Path?
    alt_gfa  : Path?
}

record HaplotypeScreen {
    id    : String
    meta  : Record
    taxid : String
    fasta : Path
}

def screenInput(id: String, meta: Record, taxid: String, fasta: Path) -> HaplotypeScreen {
    return record(id: id, meta: meta, taxid: taxid, fasta: fasta)
}

workflow DECONTAMINATE {
    take:
    assemblies   : Channel<Assembly>
    fcs_database : Path?    // Local copy of the FCS-GX database. Fetched with fcs_manifest when null.
    fcs_manifest : String   // URL of the FCS-GX database manifest
    ramdisk_path : String?  // Copy the database here before running FCS-GX

    main:
    fetched = FCSGX_FETCHDB(
        channel.of( record(id: 'fcsgx_database', manifest: fcs_manifest, output_prefix: null) )
            .filter { _x -> !fcs_database }
    )
    // The local database, or the fetched one
    database = channel.of(fcs_database)
        .filter { db -> db != null }
        .mix(fetched.map { f -> f.database })
        .collect()
        .map { dbs -> dbs.toList()[0] }

    // One screen per haplotype
    to_screen = assemblies.flatMap { a ->
        def meta = a.meta + record(
            assembly: a.meta.assembly + [ stage: 'decontaminated', build: "${a.meta.assembly.assembler}-decontaminated-${a.meta.assembly.id}" ]
        )
        def taxid = "${a.meta.sample.tax_id}"
        def haplotypes = a.alt_fasta ? [ a.pri_fasta, a.alt_fasta ] : [ a.pri_fasta ]
        haplotypes.withIndex().collect { fasta, index ->
            screenInput("${meta.assembly.build}-hap${index + 1}", meta + record(haplotype: index + 1), taxid, fasta)
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

    // Recombine the haplotypes of each assembly, primary first
    assemblies_out = cleaned
        .map { c -> tuple(c.meta.assembly.build, c) }
        .groupBy()
        .map { build, results ->
            def sorted = results.toSorted { a, b -> a.clean_fasta.name <=> b.clean_fasta.name }
            record(
                id: build,
                meta: sorted[0].meta,
                pri_fasta: sorted[0].clean_fasta,
                alt_fasta: sorted.size() > 1 ? sorted[1].clean_fasta : null,
                pri_gfa: null,
                alt_gfa: null
            )
        }

    emit:
    assemblies: Channel<Assembly> = assemblies_out
    logs      : Channel<Path>     = screened.map { r -> r.fcs_gx_report }
    screened  : Channel<Record>   = screened
    cleaned   : Channel<Record>   = cleaned
}
