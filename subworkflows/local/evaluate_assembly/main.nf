nextflow.enable.types = true

include { BUSCO_BUSCO as BUSCO       } from "../../../modules/nf-core/busco/busco/main"
include { MERQURYFK_MERQURYFK        } from "../../../modules/nf-core/merquryfk/merquryfk"
include { MERQURY_MERQURY as MERQURY } from "../../../modules/nf-core/merqury/merqury/main"
include { GFASTATS                   } from "../../../modules/nf-core/gfastats/main"

record Assembly {
    id       : String  // meta.assembly.build
    meta     : Record  // sample meta data, including meta.assembly
    pri_fasta: Path
    alt_fasta: Path?
    pri_gfa  : Path?
    alt_gfa  : Path?
}

record FastkDatabase {
    id  : String
    meta: Record
    hist: Path
    ktab: List<Path>
    prof: List<Path>
}

record MerylDatabase {
    id      : String
    meta    : Record
    meryl_db: Path
}

workflow EVALUATE_ASSEMBLY {

    take:
    assemblies            : Channel<Assembly>
    fastk_db              : Channel<FastkDatabase>
    meryl_db              : Channel<MerylDatabase>
    use_phased            : Boolean
    busco_max_lineages    : Integer
    busco_lineages_db_path: Path?

    main:
    // Kmer consistency check
    merquryfk = MERQURYFK_MERQURYFK(
        assemblies
            .map { a -> record(id: a.id, sample_id: a.meta.id, meta: a.meta, assembly: a.pri_fasta, haplotigs: a.alt_fasta) }
            .join(fastk_db.map { d -> record(sample_id: d.id, fastk_hist: d.hist, fastk_ktab: d.ktab) }, by: 'sample_id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    fastk_hist: x.fastk_hist,
                    fastk_ktab: x.fastk_ktab,
                    assembly: x.assembly,
                    haplotigs: x.haplotigs,
                    mathaptab: [],
                    pathaptab: [],
                    merquryfk_args: "-k",
                    output_prefix: "${x.meta.assembly.build}_merquryfk"
                )
            }
    )

    merqury = MERQURY(
        assemblies
            .map { a -> record(id: a.id, sample_id: a.meta.id, meta: a.meta, assembly: a.alt_fasta ? [ a.pri_fasta, a.alt_fasta ] : [ a.pri_fasta ]) }
            .join(meryl_db.map { m -> record(sample_id: m.id, meryl_db: m.meryl_db) }, by: 'sample_id')
            .map { x -> record(id: x.id, meta: x.meta, meryl_db: x.meryl_db, assembly: x.assembly, output_prefix: x.meta.id) }
    )

    // Evaluate core gene space coverage
    // Use lineages from params.busco.lineages/GOAT.
    busco = BUSCO(
        assemblies
            .flatMap { a ->
                def fastas = use_phased && a.alt_fasta ? [ a.pri_fasta, a.alt_fasta ] : [ a.pri_fasta ]
                def lineages = a.meta.settings.busco.lineages.tokenize(',').takeRight(busco_max_lineages)
                fastas.collectMany { fasta ->
                    lineages.collect { lineage ->
                        record(
                            id: "${a.meta.assembly.build}-${lineage}",
                            meta: a.meta,
                            fasta: fasta,
                            mode: 'genome',
                            lineage: lineage,
                            busco_lineages_path: busco_lineages_db_path,
                            config_file: null,
                            clean_intermediates: true,
                            busco_args: busco_lineages_db_path ? '--tar --offline' : '--tar',
                            output_prefix: "${a.meta.assembly.build}-${lineage}"
                        )
                    }.toList()
                }.toList()
            }
    )

    // Calculate contiguity stats
    gfastats = GFASTATS(
        assemblies
            .flatMap { a ->
                def fastas = use_phased && a.alt_fasta ? [ a.pri_fasta, a.alt_fasta ] : [ a.pri_fasta ]
                fastas.collect { fasta ->
                    record(
                        id: fasta.baseName,
                        meta: a.meta,
                        assembly: fasta,
                        out_fmt: null,
                        genome_size: "${a.meta.sample.genome_size}",
                        target: "",
                        agpfile: null,
                        include_bed: null,
                        exclude_bed: null,
                        instructions: null,
                        gfastats_args: "--nstar-report",
                        output_prefix: fasta.baseName
                    )
                }.toList()
            }
    )

    logs = busco.map { b -> b.short_summaries_txt }
        .mix(merquryfk.map { m -> m.stats })
        .mix(merquryfk.map { m -> m.qv })
        // Contig names are missing in first column of assembly_qv
        .mix(merquryfk.map { m -> m.full_spectra_cn_st ?: m.part_spectra_cn_st })
        .mix(merquryfk.map { m -> m.spectra_asm_st })
        .mix(merquryfk.map { m -> m.false_duplications })
        .mix(merqury.flatMap { m -> m.scaffold_qv.toList() })
        .mix(merqury.map { m -> m.spectra_cn_st_png })
        .mix(merqury.map { m -> m.spectra_asm_st_png })
        .mix(gfastats.map { g -> g.assembly_summary })
        .filter { f -> f != null }

    emit:
    logs     : Channel<Path>   = logs
    merquryfk: Channel<Record> = merquryfk
    merqury  : Channel<Record> = merqury
    busco    : Channel<Record> = busco
    gfastats : Channel<Record> = gfastats
}
