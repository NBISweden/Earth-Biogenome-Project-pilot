/*
 * Record types shared by the workflow and the stage pipelines. Import them with
 * `include { Assembly } from '<path>/types/records'`.
 *
 * Stages take one wide record per subject and return it with their outputs added
 * (`record + record(...)`). Fields a stage does not need are ignored, and optional
 * fields can be left out.
 */

nextflow.enable.types = true

// One read file, or one read pair
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

record MerylDatabase {
    id      : String
    meta    : Record
    meryl_db: Path
}

// One record per sample, used by the read stages
record Sample {
    id         : String
    meta       : Record           // id, sample, settings
    hifi       : List<Reads>?     // HiFi reads, one record per file
    hifi_merged: Reads?           // All HiFi reads of the sample in one file
    hic        : List<Reads>?     // Hi-C read pairs
    rnaseq     : List<Reads>?
    isoseq     : List<Reads>?
    fastk_hifi : FastkDatabase?
    fastk_hic  : FastkDatabase?
    meryl_hifi : MerylDatabase?
    kmer_cov   : String?          // HiFi k-mer coverage estimated by GenomeScope
}

// One record per assembly build, used by the assembly stages. Carries the sample data the stages need.
record Assembly {
    id         : String           // meta.assembly.build
    meta       : Record           // sample meta data, including meta.assembly
    pri_fasta  : Path
    alt_fasta  : Path?
    pri_gfa    : Path?
    alt_gfa    : Path?
    hifi       : List<Reads>?
    hifi_merged: Reads?
    hic        : List<Reads>?
    fastk_hifi : FastkDatabase?
    meryl_hifi : MerylDatabase?
    kmer_cov   : String?
}

// Tool options shared by the pipelines
record KmerParams {
    kmer_size: Integer?
}

record BuscoParams {
    lineages        : String?  // Comma separated BUSCO lineages, or 'auto' to retrieve them from GOAT
    lineages_db_path: Path?    // Directory with downloaded BUSCO lineages. Downloads when not set.
    max_lineages    : Integer? // Number of lineages to analyse, starting from the most specific
}
