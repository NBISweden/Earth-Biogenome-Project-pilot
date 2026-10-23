nextflow.enable.types = true

include { FASTK_HISTEX    } from "../../../modules/nf-core/fastk/histex/main"
include { GENESCOPEFK     } from "../../../modules/nf-core/genescopefk/main"
include { SMUDGEPLOT      } from "../../../modules/local/smudgeplot/main"
include { MERQURYFK_KATGC } from "../../../modules/nf-core/merquryfk/katgc/main"

record FastkDatabase {
    id  : String
    meta: Record
    hist: Path
    ktab: List<Path>
    prof: List<Path>
}

workflow GENOME_PROPERTIES {

    take:
    fastk_database: Channel<FastkDatabase>
    kmer_size     : Integer // GenomeScope k-mer length

    /* Genome properties workflow:
        - Estimate genome depth of coverage from reads
        - Generate k-mer histogram
        - Smudgeplot
    */
    main:
    // Generate GenomeScope Profile
    histex = FASTK_HISTEX(
        fastk_database.map { d -> record(id: d.id, meta: d.meta, histogram: d.hist, histex_args: '-G', output_prefix: "${d.id}_histex") }
    )
    genescope = GENESCOPEFK(
        histex.map { h ->
            record(
                id: h.id,
                meta: h.meta,
                histogram: h.hist,
                genescopefk_args: "--kmer_length ${kmer_size} --ploidy ${h.meta.sample.ploidy ?: 1}"
            )
        }
    )

    // Generate Smudgeplot
    // the genomescopefk kmer value = mean value of the heterozygous peak
    // smudgeplot -L advice = avoid cutting off the monoploid peak, thus / 2
    smudgeplot = SMUDGEPLOT(
        fastk_database
            .map { d -> record(id: d.id, ktab: d.ktab) }
            .join(genescope.map { g -> record(id: g.id, meta: g.meta, kmer_cov: g.kmer_cov) }, by: 'id')
            .map { s ->
                record(
                    id: s.id,
                    meta: s.meta,
                    fastk_ktab: s.ktab,
                    smudgeplot_threshold: Math.round(s.kmer_cov.toDouble() / 2),
                    smudgeplot_args2: "--json_report --format png",
                    output_prefix: s.id
                )
            }
    )

    // Generage GC plot
    katgc = MERQURYFK_KATGC(
        fastk_database.map { d -> record(id: d.id, meta: d.meta, fastk_hist: d.hist, fastk_ktab: d.ktab, output_prefix: "${d.id}_katgc") }
    )

    logs = genescope
        .flatMap { g -> [ g.linear_plot, g.log_plot, g.transformed_linear_plot, g.transformed_log_plot ] }
        .mix(smudgeplot.map { s -> s.smudgeplot })
        .mix(katgc.map { k -> k.stacked_gc_plot_png })
        .filter { f -> f != null }

    emit:
    kmer_cov: Channel<Record> = genescope.map { g -> record(id: g.id, meta: g.meta, kmer_cov: g.kmer_cov) }
    logs    : Channel<Path>   = logs
}
