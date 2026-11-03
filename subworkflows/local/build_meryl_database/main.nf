nextflow.enable.types = true

include { MERYL_COUNT     } from "../../../modules/nf-core/meryl/count/main"
include { MERYL_UNIONSUM  } from "../../../modules/nf-core/meryl/unionsum/main"
include { MERYL_HISTOGRAM } from "../../../modules/nf-core/meryl/histogram/main"
include { Reads; MerylDatabase } from "../../../types/records"

workflow BUILD_MERYL_DATABASE {

    take:
    fastx_data: Channel<Reads>
    library   : String  // 'hifi' or 'hic', used in file names
    kmer_size : Integer

    main:
    counted = MERYL_COUNT(
        fastx_data.map { r -> record(id: r.id, meta: r.meta, reads: r.reads, kvalue: kmer_size) }
    )

    // Union the k-mer databases of all read files with the same meta
    union = MERYL_UNIONSUM(
        counted
            .flatMap { r -> r.meryl_db.collect { db -> tuple(r.meta, db) }.toList() }
            .groupBy()
            .map { meta, dbs ->
                record(id: meta.id, meta: meta, meryl_dbs: dbs.toSorted(), kvalue: kmer_size, output_prefix: "${meta.id}_${library}")
            }
    )

    histogram = MERYL_HISTOGRAM(
        union.map { u -> record(id: u.id, meta: u.meta, meryl_db: u.meryl_db, kvalue: kmer_size) }
    )

    emit:
    histogram: Channel<Record> = histogram
    uniondb  : Channel<MerylDatabase> = union.map { u -> record(id: u.id, meta: u.meta, meryl_db: u.meryl_db) }
}
