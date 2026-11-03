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
    fastk_hic_hist : Path?                // Hi-C FastK database
    fastk_hic_ktab : List<Path>?
    fastk_hic_prof : List<Path>?
}

// One record per assembly build, used by the assembly stages. Carries the sample data the stages need.
// Fields are flat so a stage can run on its own from a YAML samplesheet.
record Assembly {
    id             : String               // meta.assembly.build
    meta           : Record               // sample meta data, including meta.assembly
    pri_fasta      : Path
    alt_fasta      : Path?
    pri_gfa        : Path?                // Contig graph of the raw assembly. The stages after it have none
    hifi           : List<Path>?          // HiFi reads, one file per element
    hifi_merged    : Path?                // All HiFi reads of the sample in one file
    hic            : List<List<Path>>?    // Hi-C read pairs
    rnaseq         : List<List<Path>>?    // RNAseq reads: a pair, or a single file
    fastk_hifi_hist: Path?                // HiFi FastK database
    fastk_hifi_ktab: List<Path>?
    fastk_hifi_prof: List<Path>?
    meryl_hifi_db  : Path?                // HiFi Meryl database
}

// The assembly of a build with the files that a stage made: the primary and, when there is a second one,
// the alternate. The contig graph belongs to the assembly before the stage, so it is dropped.
// `+` keeps the other fields of the record, also those that a stage added to join on.
def withFastas(a: Assembly, fastas: List<Path>) {
    return a + record(pri_fasta: fastas[0], alt_fasta: fastas.size() > 1 ? fastas[1] : null, pri_gfa: null as Path)
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

record FcsParams {
    database    : Path?   // Local copy of the FCS-GX database. Fetched with the manifest when not set.
    manifest    : String? // URL of the FCS-GX database manifest
    ramdisk_path: String? // Copy the database here before running FCS-GX
}

// The options of the assembly QC. The stage pipelines take them as one param and pass them on.
record QcParams {
    busco     : BuscoParams
    fastk     : KmerParams
    meryl     : KmerParams
    use_phased: Boolean // Evaluate (and process, in the stages that do) the primary and the alternate assembly
}

def defaultQcParams() {
    return record(
        busco: record(lineages: null as String, lineages_db_path: null as Path, max_lineages: null as Integer),
        fastk: record(kmer_size: null as Integer),
        meryl: record(kmer_size: null as Integer),
        use_phased: false
    )
}
