include { combineByMetaKeys    } from "../../../modules/local/functions"
include { getEachAssembly      } from "../../../modules/local/functions"
include { getPrimaryAssembly   } from "../../../modules/local/functions"
include { BUSCO_BUSCO as BUSCO } from "../../../modules/nf-core/busco/busco/main"
include { MERQURYFK_MERQURYFK  } from "../../../modules/nf-core/merquryfk/merquryfk"
include { MERQURY              } from "../../../modules/nf-core/merqury/main"
include { GFASTATS             } from "../../../modules/nf-core/gfastats/main"

workflow EVALUATE_ASSEMBLY {

    take:
    assembly_ch        // input type: [ meta, [ id:'assemblerX_build1', pri_fasta: '/path/to/primary_asm', alt_fasta: '/path/to/alternate_asm' ] ]
    fastk_db           // input type: [ meta, [ 'path/to/reads.hist' ], [ '/path/to/reads.ktab' ] ]
    meryl_db           // input type: [ meta, 'path/to/reads.union.meryl' ]

    main:

    // Kmer consistency check
    MERQURYFK_MERQURYFK (
        combineByMetaKeys (
            fastk_db,
            getEachAssembly(assembly_ch),
            keySet: ['id','sample'],
            meta: 'rhs'
        ) // TODO: Change the function so we don't need the map
        .map{ meta, hist, ktab, assembly ->
            assembly instanceof Path ?
            tuple(meta, hist, ktab, assembly, []) :
            tuple(meta, hist, ktab, assembly.head(), assembly.last())
        }, // [ meta, hist, ktab, assembly ]
        [],
        [],
    )

    MERQURY (
        combineByMetaKeys (
            meryl_db,
            getEachAssembly(assembly_ch),
            keySet: ['id','sample'],
            meta: 'rhs'
        ) // [ meta, meryldb, assembly ]
    )

    // Evaluate core gene space coverage
    busco_input = getPrimaryAssembly(assembly_ch)
        .flatMap { meta, asm ->
            if ( meta.settings?.busco?.lineages ) {
                // Use lineages from params.busco.lineages/GOAT.
                meta.settings.busco.lineages.tokenize(',').collect{ [ meta, asm, it ] }
            } else {
                // auto-detect.
                [ [ meta, asm, 'auto' ] ]
            }
        }
        .dump(tag: 'BUSCO', pretty: true)
        .multiMap { meta, asm, line ->
            assembly_ch: [ meta, asm ]
            lineage_ch: line
        }
    BUSCO (
        busco_input.assembly_ch,
        'genome',
        busco_input.lineage_ch,
        params.busco.lineages_db_path ? file( params.busco.lineages_db_path, checkIfExists: true ) : [],
        [],
        true
    )

    // Calculate contiguity stats
    fasta_ch = assembly_ch.multiMap { meta, assembly ->
        fasta: [ meta, assembly.pri_fasta ]
        genome_size: meta.sample.genome_size
    }
    GFASTATS(
        fasta_ch.fasta,
        [],                   // output format: none
        fasta_ch.genome_size, // genome size
        "",                   // target
        [[],[]],              // AGP file
        [[],[]],              // include bed
        [[],[]],              // exclude bed
        [[],[]]               // SAK instructions
    )

    BUSCO.out.short_summaries_txt
        .mix(
            MERQURYFK_MERQURYFK.out.stats,
            MERQURYFK_MERQURYFK.out.qv,
            // MERQURYFK_MERQURYFK.out.assembly_qv, // Contig names are missing in first column
            MERQURYFK_MERQURYFK.out.full_spectra_cn_st.join(
                MERQURYFK_MERQURYFK.out.part_spectra_cn_st, remainder: true
            ).map { meta, full, part -> full ? tuple(meta, full): tuple(meta, part) },
            MERQURYFK_MERQURYFK.out.spectra_asm_st,
            MERQURYFK_MERQURYFK.out.false_duplications,
            MERQURY.out.scaffold_qv,
            MERQURY.out.spectra_cn_st_png,
            MERQURY.out.spectra_asm_st_png,
            GFASTATS.out.assembly_summary,
        )
        .map { _meta, file -> file }
        .set { logs }
    MERQURYFK_MERQURYFK.out.versions.first().mix(
        GFASTATS.out.versions.first(),
        MERQURY.out.versions.first(),
        BUSCO.out.versions.first()
    ).set { versions }

    emit:
    logs
    versions
}