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

// One record per sample, used by the read stages. Fields are flat, like those of Assembly.
record Sample {
    id             : String
    meta           : Record               // id, sample, settings
    hifi           : List<Path>?          // HiFi reads, one file per element
    hifi_merged    : Path?                // All HiFi reads of the sample in one file
    hic            : List<List<Path>>?    // Hi-C read pairs
    rnaseq         : List<List<Path>>?    // RNAseq reads: a pair, or a single file
    isoseq         : List<Path>?
    fastk_hifi_hist: Path?                // HiFi FastK database
    fastk_hifi_ktab: List<Path>?
    fastk_hifi_prof: List<Path>?
    meryl_hifi_db  : Path?                // HiFi Meryl database
    kmer_cov       : String?              // HiFi k-mer coverage estimated by GenomeScope
}

// One record per assembly build, used by the assembly stages. Carries the sample data the stages need.
// Fields are flat so a stage can run on its own from a YAML samplesheet.
record Assembly {
    id             : String               // meta.assembly.build
    meta           : Record               // sample meta data, including meta.assembly
    pri_fasta      : Path
    alt_fasta      : Path?
    pri_gfa        : Path?
    alt_gfa        : Path?
    hifi           : List<Path>?          // HiFi reads, one file per element
    hifi_merged    : Path?                // All HiFi reads of the sample in one file
    hic            : List<List<Path>>?    // Hi-C read pairs
    rnaseq         : List<List<Path>>?    // RNAseq reads: a pair, or a single file
    fastk_hifi_hist: Path?                // HiFi FastK database
    fastk_hifi_ktab: List<Path>?
    fastk_hifi_prof: List<Path>?
    meryl_hifi_db  : Path?                // HiFi Meryl database
    kmer_cov       : String?              // HiFi k-mer coverage estimated by GenomeScope
}

// A new Assembly from the fields of a record, with some of them replaced. Stages build the records they
// return with this. A record with extra fields cannot be passed to a function that takes an Assembly,
// and an Assembly instance cannot be extended with `+`.
def assemblyWith(a: Record, changes: Map) -> Assembly {
    return record(
        id             : ( changes.containsKey('id') ? changes.id : a.id ) as String,
        meta           : ( changes.containsKey('meta') ? changes.meta : a.meta ) as Record,
        pri_fasta      : ( changes.containsKey('pri_fasta') ? changes.pri_fasta : a.pri_fasta ) as Path,
        alt_fasta      : ( changes.containsKey('alt_fasta') ? changes.alt_fasta : a.alt_fasta ) as Path,
        pri_gfa        : ( changes.containsKey('pri_gfa') ? changes.pri_gfa : a.pri_gfa ) as Path,
        alt_gfa        : ( changes.containsKey('alt_gfa') ? changes.alt_gfa : a.alt_gfa ) as Path,
        hifi           : ( changes.containsKey('hifi') ? changes.hifi : a.hifi ) as List<Path>,
        hifi_merged    : ( changes.containsKey('hifi_merged') ? changes.hifi_merged : a.hifi_merged ) as Path,
        hic            : ( changes.containsKey('hic') ? changes.hic : a.hic ) as List<List<Path>>,
        rnaseq         : ( changes.containsKey('rnaseq') ? changes.rnaseq : a.rnaseq ) as List<List<Path>>,
        fastk_hifi_hist: ( changes.containsKey('fastk_hifi_hist') ? changes.fastk_hifi_hist : a.fastk_hifi_hist ) as Path,
        fastk_hifi_ktab: ( changes.containsKey('fastk_hifi_ktab') ? changes.fastk_hifi_ktab : a.fastk_hifi_ktab ) as List<Path>,
        fastk_hifi_prof: ( changes.containsKey('fastk_hifi_prof') ? changes.fastk_hifi_prof : a.fastk_hifi_prof ) as List<Path>,
        meryl_hifi_db  : ( changes.containsKey('meryl_hifi_db') ? changes.meryl_hifi_db : a.meryl_hifi_db ) as Path,
        kmer_cov       : ( changes.containsKey('kmer_cov') ? changes.kmer_cov : a.kmer_cov ) as String
    )
}

def assemblyOf(a: Record) -> Assembly {
    return assemblyWith(a, [:])
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
