nextflow.enable.types = true

include { MERQURYFK_KATCOMP } from "../../../modules/nf-core/merquryfk/katcomp/main"

record LibraryComparison {
    id         : String
    meta       : Record
    fastk1_hist: Path
    fastk1_ktab: List<Path>
    fastk2_hist: Path
    fastk2_ktab: List<Path>
}

workflow COMPARE_LIBRARIES {

    take:
    comparisons: Channel<LibraryComparison> // Two k-mer databases per sample

    main:
    katcomp = MERQURYFK_KATCOMP(
        comparisons.map { c ->
            record(
                id: c.id,
                meta: c.meta,
                fastk1_hist: c.fastk1_hist,
                fastk1_ktab: c.fastk1_ktab,
                fastk2_hist: c.fastk2_hist,
                fastk2_ktab: c.fastk2_ktab,
                output_prefix: "${c.id}_katcomp"
            )
        }
    )

    emit:
    logs   : Channel<Path>   = katcomp.map { k -> k.stacked_png }.filter { f -> f != null }
    katcomp: Channel<Record> = katcomp
}
