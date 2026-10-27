nextflow.enable.types = true

include { FASTK_FASTK } from "../../../modules/nf-core/fastk/fastk/main"
include { FASTK_MERGE } from "../../../modules/nf-core/fastk/merge/main"

record Reads {
    id   : String
    meta : Record
    reads: List<Path>
}

record FastkDatabase {
    id  : String
    meta: Record
    hist: Path
    ktab: List<Path>
    prof: List<Path>
}

def fastkDatabase(id: String, meta: Record, hist: Path, ktab: List<Path>, prof: List<Path>) -> FastkDatabase {
    return record(id: id, meta: meta, hist: hist, ktab: ktab, prof: prof)
}

workflow BUILD_FASTK_DATABASE {

    take:
    fastx_data: Channel<Reads>
    library   : String  // 'hifi' or 'hic', used in file names
    kmer_size : Integer

    main:
    counted = FASTK_FASTK(
        fastx_data.map { r ->
            record(
                id: r.id,
                meta: r.meta,
                reads: r.reads,
                fastk_args: "-t1 -k${kmer_size} -v",
                output_prefix: "${r.id}_${r.reads[0].baseName}_${library}"
            )
        }
    )

    // Group the k-mer tables of all read files (or pairs) of a sample
    grouped = counted
        .map { r -> tuple(record(id: r.meta.id, sample: r.meta.sample, settings: r.meta.settings), r) }
        .groupBy()
        .map { key, results ->
            record(
                id: key.id,
                meta: key,
                hist: results.collect { x -> x.hist }.toSorted(),
                ktab: results.collectMany { x -> x.ktab.toList() }.toSorted(),
                prof: results.collectMany { x -> x.prof.toList() }.toSorted()
            )
        }

    single = grouped
        .filter { g -> g.hist.size() == 1 }
        .map { g -> fastkDatabase(g.id, g.meta, g.hist[0], g.ktab, g.prof) }

    merged = FASTK_MERGE(
        grouped
            .filter { g -> g.hist.size() > 1 }
            .map { g ->
                record(
                    id: g.id,
                    meta: g.meta,
                    hist: g.hist,
                    ktab: g.ktab,
                    prof: g.prof,
                    fastk_args: g.ktab ? "-ht" : "-h",
                    output_prefix: "${g.id}_merged_${library}"
                )
            }
    )

    database = single.mix(
        merged.map { m -> fastkDatabase(m.id, m.meta, m.hist, m.ktab.toList().toSorted(), m.prof.toList().toSorted()) }
    )

    emit:
    database    : Channel<FastkDatabase> = database
    fastk_runs  : Channel<Record>        = counted
    fastk_merged: Channel<Record>        = merged
}
