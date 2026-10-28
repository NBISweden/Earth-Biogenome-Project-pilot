#! /usr/bin/env nextflow

/*
 * Development: See docs/development to understand the workflow programming model and
 * how channel contents are structured.
 */

nextflow.enable.types = true

// Data preparation
include { PREPARE_INPUT                                     } from "./subworkflows/local/prepare_input/main"
include { BUILD_FASTK_DATABASE as BUILD_FASTK_HIFI_DATABASE } from "./subworkflows/local/build_fastk_database/main"
include { BUILD_FASTK_DATABASE as BUILD_FASTK_HIC_DATABASE  } from "./subworkflows/local/build_fastk_database/main"
include { BUILD_MERYL_DATABASE as BUILD_MERYL_HIFI_DATABASE } from "./subworkflows/local/build_meryl_database/main"
include { BUILD_MERYL_DATABASE as BUILD_MERYL_HIC_DATABASE  } from "./subworkflows/local/build_meryl_database/main"
include { CONVERT_FASTQ_CRAM                                } from "./subworkflows/local/02_convert_fastq_cram/main"
// Stages
include { INSPECT_DATA                                      } from "./subworkflows/local/01_inspect_data/main"
include { ASSEMBLE                                          } from "./subworkflows/local/03_assemble/main"
include { DECONTAMINATE                                     } from "./subworkflows/local/04_decontaminate/main"
include { PURGE_DUPLICATES                                  } from "./subworkflows/local/05_purge_dups/main"
include { DVPOLISH                                          } from "./subworkflows/local/06_polish/main"
include { SCAFFOLD                                          } from "./subworkflows/local/07_scaffold/main"
include { SCAFFOLD_CURATION                                 } from "./subworkflows/local/08_scaffold_curation/main"
include { ALIGN_RNASEQ                                      } from "./subworkflows/local/09_align_rnaseq/main"
// Evaluation and reports
include { COMPARE_ASSEMBLIES                                } from "./subworkflows/local/compare_assemblies/main"
include { EVALUATE_ASSEMBLY                                 } from "./subworkflows/local/evaluate_assembly/main"
include { REPORT_VERSIONS                                   } from "./subworkflows/local/version_report/main"
include { ASSEMBLY_REPORT                                   } from "./subworkflows/local/assembly_report/main"

record Reads {
    id   : String
    meta : Record
    reads: List<Path>
}

record Assembly {
    id       : String  // meta.assembly.build
    meta     : Record  // sample meta data, including meta.assembly
    pri_fasta: Path
    alt_fasta: Path?
    pri_gfa  : Path?
    alt_gfa  : Path?
}

record KmerCoverage {
    id      : String
    kmer_cov: String?
}

record KmerParams {
    kmer_size: Integer?
}

record BuscoParams {
    lineages        : String? // A comma separated set of lineages or 'auto' (default: Retrieve from GOAT.)
    lineages_db_path: Path?
    max_lineages    : Integer? // The number of busco lineages to analyse starting from the most specific. (default: 2)
}

record FcsParams {
    database    : Path?
    manifest    : String?
    ramdisk_path: String?
}

// Override by supplying a params.yml file to the -params-file option in Nextflow.
params {

    // Workflow inputs
    input    : Path
    reference: Path? // Reference genome to compare the assemblies to

    // Stages to run. Output subdirectory for each stage.
    steps                  : String  = 'inspect,preprocess,assemble,purge,polish,screen,scaffold,curate'
    nuclear_assembly_mode  : Boolean = true      // values: true, false
    organelle_assembly_mode: String  = 'contigs' // values: 'reads', 'contigs', 'none'
    use_phased             : Boolean = false
    stage_dirs             : Map<String,String> = [
        inspect   : '01_read_inspection',
        preprocess: '02_read_preprocessing',
        assembly  : '03_assembly',
        screen    : '04_contamination_screen',
        purge     : '05_duplicate_purging',
        polish    : '06_error_polishing',
        scaffold  : '07_scaffolding',
        curate    : '08_curation',
        alignRNA  : '09_RNA_alignment',
        report    : '10_report',
    ]

    // Tool specific
    fastk      : KmerParams = record(kmer_size: null) // default: 31
    meryl      : KmerParams = record(kmer_size: null) // default: 31
    genescopefk: KmerParams = record(kmer_size: null) // default: 31
    hifiasm    : Map<String,String>? // Expects a Map of options [ opts01: "--opt1 --opt2 val", opts02: "--opt2 val --opt3" ]
    hifiasm_with_hic: Boolean = false
    busco      : BuscoParams = record(lineages: null, lineages_db_path: null, max_lineages: null)
    oatkdb     : Path?    // Setup: git clone --depth 1 https://github.com/c-zhou/OatkDB.git . Then oatkdb = "/path/to/OatkDB/v20230921"
    chunk_size : String   = '100.MB'  // DVPOLISH
    fcs        : FcsParams = record(database: null, manifest: null, ramdisk_path: null)

    // Reports
    quarto_assembly_report          : Path   = file("${projectDir}/assets/notebooks/assembly_report.qmd")
    quarto_assembly_report_aux_files: String = "${projectDir}/assets/notebooks/_*.qmd"
    multiqc_assembly_report_config  : Path   = file("${projectDir}/configs/multiqc_assembly_report_config.yml")

    // Hi-C
    hic_type         : String  = 'arima_v2' // define HiC sequencing type [ "arima_v2", "arima_v1", "omni-c", "" ]
    hic_bam_sorting  : String  = 'by_name'  // define how to sort the final HiC bam file: [ by_coordinate, by_name ]
    cooler_bin_size  : Integer = 1000       // define 1000bp bin size for cooler cload intervals
    hic_map_sort_by  : String  = 'length'   // define how to sort the HiC-maps: [ length, name, nosort ]
    hic_map_qv       : Integer = 0          // define minimum map quality for HiC alignments
    hifi_coverage_cap: Integer = 100        // pacbio hifi coverage cap for pretext and higlass coverage track, e.g. cap repeats at <hifi_coverage_cap>
    telomer_motif    : String  = "TTAGGG"
}

// Change the stage of an assembly. The build is used in file names.
def setAssemblyStage(a: Assembly, stage: String) -> Assembly {
    def build = "${a.meta.assembly.assembler}-${stage}-${a.meta.assembly.id}"
    return record(
        id       : build,
        meta     : a.meta + record(assembly: a.meta.assembly + [ stage: stage, build: build ]),
        pri_fasta: a.pri_fasta,
        alt_fasta: a.alt_fasta,
        pri_gfa  : a.pri_gfa,
        alt_gfa  : a.alt_gfa
    )
}

// Size in bytes of a memory string such as '100.MB' or '2 GB'
def chunkSizeBytes(size: String) -> Long {
    def units = [ B: 1L, KB: 1024L, MB: 1024L * 1024L, GB: 1024L * 1024L * 1024L, TB: 1024L * 1024L * 1024L * 1024L ]
    def parts = size.trim().toUpperCase().replace('.', ' ').tokenize(' ')
    return parts[0].toLong() * ( parts.size() > 1 ? units[parts[1]] : 1L )
}

def kmerCoverage(id: String, kmer_cov: String?) -> KmerCoverage {
    return record(id: id, kmer_cov: kmer_cov)
}

// Directory of the stage that produced an assembly
def assemblyStageDir(stage_dirs: Map<String,String>, stage: String) -> String {
    def dirs = [ raw: 'assembly', decontaminated: 'screen', purged: 'purge', polished: 'polish', scaffolded: 'scaffold', curated: 'curate' ]
    return stage_dirs[dirs[stage]]
}

def workflowVersion() -> String {
    def version = workflow.manifest.version ? "v${workflow.manifest.version}" : ""
    def sha = workflow.commitId ? "-g${workflow.commitId.substring(0, 7)}" : ""
    return "Workflow:\n  ${workflow.manifest.name}: ${version}${sha}\n  Nextflow: ${nextflow.version}"
}

workflow {

    main:
    // Define constants
    def permitted_stages = [
        'inspect',      // 01 - Read inspection
        'preprocess',   // 02 - Read preprocessing
        'assemble',     // 03 - Assembly
        'screen',       // 04 - Contamination screening
        'purge',        // 05 - Duplicate purging
        'polish',       // 06 - Error polishing
        'scaffold',     // 07 - Scaffolding
        'curate',       // 08 - Rapid curation
        'alignRNA'      // 09 - Align RNAseq data
    ]

    // Check input
    def steps = params.steps.tokenize(",")
    if ( !steps.every { step -> step in permitted_stages } ) {
        error "Unrecognised workflow step in ${params.steps} ( ${permitted_stages} )"
    }
    def enable = record(
        inspect   : 'inspect' in steps,
        preprocess: 'preprocess' in steps,
        assemble  : 'assemble' in steps,
        screen    : 'screen' in steps,
        purge     : 'purge' in steps,
        polish    : 'polish' in steps,
        scaffold  : 'scaffold' in steps,
        curate    : 'curate' in steps,
        alignRNA  : 'alignRNA' in steps
    )

    // The primary workflow for the Earth Biogenome Project Pilot
    log.info("""
    Running NBIS Earth Biogenome Project Assembly workflow.
    """)

    def fastk_kmer_size = params.fastk.kmer_size ?: 31
    def meryl_kmer_size = params.meryl.kmer_size ?: 31
    def genescopefk_kmer_size = params.genescopefk.kmer_size ?: 31

    // Read in data
    input = PREPARE_INPUT( params.input, params.busco.lineages )

    // Create Cram and trim
    convert = CONVERT_FASTQ_CRAM( input.hic, params.hic_type.startsWith('arima'), params.hic_type )

    // Build necessary databases
    // TODO: Migrate back to Meryl. Genome inspection missing KATGC and SMUDGEPLOT for meryldb
    fastk_hifi = BUILD_FASTK_HIFI_DATABASE( input.hifi, 'hifi', fastk_kmer_size )
    fastk_hic  = BUILD_FASTK_HIC_DATABASE( convert.fastq, 'hic', fastk_kmer_size )
    meryl_hifi = BUILD_MERYL_HIFI_DATABASE( input.hifi, 'hifi', meryl_kmer_size )
    BUILD_MERYL_HIC_DATABASE( convert.fastq, 'hic', meryl_kmer_size )

    // Data inspection
    inspect = INSPECT_DATA(
        input.hifi.filter { _r -> enable.inspect },
        input.hic.filter { _r -> enable.inspect },
        fastk_hifi.database.filter { _d -> enable.inspect },
        fastk_hic.database.filter { _d -> enable.inspect },
        genescopefk_kmer_size
    )
    // Purging waits for the k-mer coverage of the sample when it is estimated
    kmer_cov = enable.inspect
        ? inspect.hifi_kmer_cov.map { k -> kmerCoverage(k.id, k.kmer_cov) }
        : input.hifi_merged.map { r -> kmerCoverage(r.id, null) }

    // Assemble
    assemble = ASSEMBLE(
        input.hifi_merged.filter { _r -> enable.assemble },
        input.hic.filter { _r -> enable.assemble },
        params.oatkdb,
        params.nuclear_assembly_mode,
        params.organelle_assembly_mode,
        params.hifiasm ?: [:],
        params.hifiasm_with_hic,
        params.use_phased
    )
    raw_assemblies = assemble.raw_assemblies.mix(input.assemblies.filter { a -> a.meta.assembly.stage == 'raw' })

    // Contamination screen
    to_screen = raw_assemblies.map { a -> setAssemblyStage(a, 'decontaminated') }
    screen = DECONTAMINATE(
        to_screen.filter { _a -> enable.screen },
        params.fcs.database,
        params.fcs.manifest ?: "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest",
        params.fcs.ramdisk_path != null ? params.fcs.ramdisk_path : '/dev/shm/gxdb'
    )
    cleaned_assemblies = to_screen.filter { _a -> !enable.screen }
        .mix(screen.assemblies.filter { _a -> enable.screen })
        .mix(input.assemblies.filter { a -> a.meta.assembly.stage == 'decontaminated' })

    // Purge duplicates
    to_purge = cleaned_assemblies.map { a -> setAssemblyStage(a, 'purged') }
    purge = PURGE_DUPLICATES(
        to_purge.filter { _a -> enable.purge },
        input.hifi_merged,
        kmer_cov
    )
    purged_assemblies = to_purge.filter { _a -> !enable.purge }
        .mix(purge.assemblies.filter { _a -> enable.purge })
        .mix(input.assemblies.filter { a -> a.meta.assembly.stage == 'purged' })

    // Polish
    to_polish = purged_assemblies.map { a -> setAssemblyStage(a, 'polished') }
    polish = DVPOLISH(
        to_polish.filter { _a -> enable.polish },
        input.hifi_merged,
        meryl_hifi.uniondb,
        chunkSizeBytes(params.chunk_size),
        params.use_phased
    )
    polished_assemblies = to_polish.filter { _a -> !enable.polish }
        .mix(polish.assemblies.filter { _a -> enable.polish })
        .mix(input.assemblies.filter { a -> a.meta.assembly.stage == 'polished' })

    // Scaffold
    to_scaffold = polished_assemblies.map { a -> setAssemblyStage(a, 'scaffolded') }
    scaffold = SCAFFOLD(
        to_scaffold.filter { _a -> enable.scaffold },
        convert.fastq,
        params.use_phased,
        params.hic_type,
        params.hic_bam_sorting
    )
    scaffolded_assemblies = to_scaffold.filter { _a -> !enable.scaffold }
        .mix(scaffold.assemblies.filter { _a -> enable.scaffold })
        .mix(input.assemblies.filter { a -> a.meta.assembly.stage == 'scaffolded' })

    // Curate
    to_curate = scaffolded_assemblies.map { a -> setAssemblyStage(a, 'curated') }
    curate = SCAFFOLD_CURATION(
        to_curate.filter { _a -> enable.curate },
        convert.fastq,
        input.hifi,
        params.hic_map_sort_by,
        params.hic_map_qv,
        params.cooler_bin_size,
        params.hifi_coverage_cap,
        params.telomer_motif
    )

    // Align RNAseq
    rnaseq = ALIGN_RNASEQ(
        input.rnaseq.filter { _r -> enable.alignRNA },
        input.assemblies.filter { _a -> enable.alignRNA }
    )

    // Evaluate assemblies
    evaluated_assemblies = ( enable.assemble ? assemble.raw_assemblies : channel.empty() )
        .mix(input.assemblies)
        .mix(enable.screen ? screen.assemblies : channel.empty())
        .mix(enable.purge ? purge.assemblies : channel.empty())
        .mix(enable.polish ? polish.assemblies : channel.empty())
        .mix(enable.scaffold ? scaffold.assemblies : channel.empty())
    compared = COMPARE_ASSEMBLIES( evaluated_assemblies, params.reference )
    evaluated = EVALUATE_ASSEMBLY(
        evaluated_assemblies.collect().flatMap { l -> l.toList() }, // Introduce bottleneck to delay evaluation until previous steps are performed
        fastk_hifi.database,
        meryl_hifi.uniondb,
        params.use_phased,
        params.busco.max_lineages ?: 2,
        params.busco.lineages_db_path
    )

    // Version reporting
    versions = REPORT_VERSIONS(
        channel.topic("versions") as Channel<Tuple>,
        workflowVersion()
    )

    // Assembly report
    report_logs = channel.of(params.multiqc_assembly_report_config)
        .mix(inspect.logs)
        .mix(assemble.logs)
        .mix(screen.logs)
        .mix(purge.logs)
        .mix(polish.logs)
        .mix(scaffold.logs)
        .mix(compared.logs)
        .mix(evaluated.logs)
    report = ASSEMBLY_REPORT(
        input.metadata,
        params.quarto_assembly_report,
        files(params.quarto_assembly_report_aux_files, checkIfExists: true).toList(),
        report_logs,
        versions,
        [
            diagnostics            : 'debug' in workflow.profile.tokenize(","),
            nuclear_assembly_mode  : params.nuclear_assembly_mode,
            organelle_assembly_mode: params.organelle_assembly_mode in ['contigs', 'reads'],
            inspect                : enable.inspect,
            preprocess             : enable.preprocess,
            assemble               : enable.assemble,
            screen                 : enable.screen,
            purge                  : enable.purge,
            polish                 : enable.polish,
            scaffold               : enable.scaffold,
            curate                 : enable.curate,
            alignRNA               : enable.alignRNA
        ]
    )

    publish:
    // 01 - Read inspection
    goat              = input.goat
    seqkit_hifi_stats = inspect.seqkit_hifi
    seqkit_hic_stats  = inspect.seqkit_hic
    fastqc_hic        = inspect.fastqc
    fastk_runs        = fastk_hifi.fastk_runs.mix(fastk_hic.fastk_runs)
    fastk_merged      = fastk_hifi.fastk_merged.mix(fastk_hic.fastk_merged)
    genescopefk       = inspect.genescopefk
    smudgeplot        = inspect.smudgeplot
    katgc             = inspect.katgc
    kat_comp          = inspect.katcomp
    dtol_search       = report.dtol
    // 02 - Read preprocessing
    converted_bam_to_fasta = input.converted
    hic_cram          = convert.cram
    hic_crai          = convert.crai
    // 03 - Assembly
    hifiasm           = assemble.hifiasm
    gfa2fa            = assemble.gfa2fa
    mitohifi_references = assemble.mito_reference
    mitohifi          = assemble.mitohifi
    oatk              = assemble.oatk
    dotplots          = assemble.dotplots
    // 04 - Contamination screen
    fcsgx_reports     = screen.screened
    fcsgx_clean       = screen.cleaned
    // 05 - Duplicate purging
    purge_dups        = purge.purged
    // 06 - Error polishing
    dvpolish          = polish.final_asm
    // 07 - Scaffolding
    pairtools         = scaffold.pairtools
    yahs              = scaffold.yahs
    merged_haps       = scaffold.merged_haplotypes
    // 08 - Curation
    higlass_mcool     = curate.zoomify
    higlass_gaps      = curate.gaps
    higlass_telomers  = curate.telomer_bigwig
    pretext           = curate.pretext_tracks
    // 09 - RNA alignment
    star              = rnaseq
    // Evaluation
    merquryfk         = evaluated.merquryfk
    merqury           = evaluated.merqury
    busco             = evaluated.busco
    gfastats          = evaluated.gfastats
    quast             = compared.quast
    // 10 - Report
    report            = report.report

    onComplete:
    if ( workflow.success ) {
        log.info("""
        Thank you for using the NBIS Earth Biogenome Project Assembly workflow.
        The workflow completed successfully.

        Results are located in the folder: ${workflow.outputDir}
        """)
    } else {
        log.info("""
        Thank you for using the NBIS Earth Biogenome Project Assembly workflow.
        The workflow completed unsuccessfully.

        Please read over the error message. If you are unable to solve it, please
        post an issue at https://github.com/NBISweden/Earth-Biogenome-Project-pilot/issues
        where we will do our best to help.
        """)
    }
}

output {
    // 01 - Read inspection
    goat: Channel<Record> {
        path "${params.stage_dirs.inspect}/goat"
    }
    seqkit_hifi_stats: Channel<Record> {
        path "${params.stage_dirs.inspect}/seqkit_hifi_stats"
    }
    seqkit_hic_stats: Channel<Record> {
        path "${params.stage_dirs.inspect}/seqkit_hic_stats"
    }
    fastqc_hic: Channel<Record> {
        path "${params.stage_dirs.inspect}/fastqc_hic"
    }
    fastk_runs: Channel<Record> {
        path "${params.stage_dirs.inspect}/fastk"
    }
    fastk_merged: Channel<Record> {
        path "${params.stage_dirs.inspect}/fastk"
    }
    genescopefk: Channel<Record> {
        path "${params.stage_dirs.inspect}/genescopefk"
    }
    smudgeplot: Channel<Record> {
        path "${params.stage_dirs.inspect}/smudgeplot"
    }
    katgc: Channel<Record> {
        path "${params.stage_dirs.inspect}/katgc"
    }
    kat_comp: Channel<Record> {
        path "${params.stage_dirs.inspect}/kat_comp"
    }
    dtol_search: Channel<Record> {
        path { r -> r.json_file >> ("${params.stage_dirs.inspect}/dtol_search/") }
    }
    // 02 - Read preprocessing
    converted_bam_to_fasta: Channel<Record> {
        path "${params.stage_dirs.preprocess}/converted_bam_to_fasta"
    }
    hic_cram: Channel<Record> {
        path "${params.stage_dirs.preprocess}/hi-c_cram"
    }
    hic_crai: Channel<Record> {
        path "${params.stage_dirs.preprocess}/hi-c_cram"
    }
    // 03 - Assembly
    hifiasm: Channel<Record> {
        path { r -> "${params.stage_dirs.assembly}/${r.meta.assembly.build}" }
    }
    gfa2fa: Channel<Record> {
        path { r -> "${params.stage_dirs.assembly}/${r.meta.assembly.build}" }
    }
    mitohifi_references: Channel<Record> {
        path "${params.stage_dirs.assembly}/organelle/mitohifi/references"
    }
    mitohifi: Channel<Record> {
        path { r -> "${params.stage_dirs.assembly}/organelle/mitohifi/${params.organelle_assembly_mode == 'contigs' ? r.meta.assembly.build : r.meta.id}" }
    }
    oatk: Channel<Record> {
        path "${params.stage_dirs.assembly}/organelle/oatk"
    }
    dotplots: Channel<Record> {
        path "${params.stage_dirs.assembly}/organelle/dotplots"
    }
    // 04 - Contamination screen
    fcsgx_reports: Channel<Record> {
        path "${params.stage_dirs.screen}/fcsgx"
    }
    fcsgx_clean: Channel<Record> {
        path "${params.stage_dirs.screen}/fcsgx"
    }
    // 05 - Duplicate purging
    purge_dups: Channel<Record> {
        path "${params.stage_dirs.purge}/purge_dups"
    }
    // 06 - Error polishing
    dvpolish: Channel<Record> {
        path { r -> "${params.stage_dirs.polish}/dvpolish/${r.meta.assembly.build}" }
    }
    // 07 - Scaffolding
    pairtools: Channel<Record> {
        path { r -> r.stat >> ("${params.stage_dirs.scaffold}/pairtools/${r.meta.assembly.build}/") }
    }
    yahs: Channel<Record> {
        path { r -> "${params.stage_dirs.scaffold}/yahs/${r.meta.assembly.build}" }
    }
    merged_haps: Channel<Record> {
        path { r -> "${params.stage_dirs.scaffold}/merged_haps/${r.meta.assembly.build}" }
    }
    // 08 - Curation
    higlass_mcool: Channel<Record> {
        path { r -> r.mcool >> ("${params.stage_dirs.curate}/higlass/${r.meta.assembly.build}/") }
    }
    higlass_gaps: Channel<Record> {
        path { r -> r.beddb >> ("${params.stage_dirs.curate}/higlass/${r.meta.assembly.build}/") }
    }
    higlass_telomers: Channel<Record> {
        path { r -> r.bw >> ("${params.stage_dirs.curate}/higlass/${r.meta.assembly.build}/") }
    }
    pretext: Channel<Record> {
        path { r -> r.hitile >> ("${params.stage_dirs.curate}/pretext/${r.meta.assembly.build}/") }
    }
    // 09 - RNA alignment
    star: Channel<Record> {
        path { r -> "${params.stage_dirs.alignRNA}/star/${r.meta.assembly.build}" }
    }
    // Evaluation
    merquryfk: Channel<Record> {
        path { r -> "${assemblyStageDir(params.stage_dirs, r.meta.assembly.stage)}/merquryfk/${r.meta.assembly.build}" }
    }
    merqury: Channel<Record> {
        path { r -> "${assemblyStageDir(params.stage_dirs, r.meta.assembly.stage)}/merqury/${r.meta.assembly.build}" }
    }
    busco: Channel<Record> {
        path { r ->
            r.short_summaries_txt >> ("${assemblyStageDir(params.stage_dirs, r.meta.assembly.stage)}/busco/${r.meta.assembly.build}/")
            r.short_summaries_json >> ("${assemblyStageDir(params.stage_dirs, r.meta.assembly.stage)}/busco/${r.meta.assembly.build}/")
        }
    }
    gfastats: Channel<Record> {
        path { r -> "${assemblyStageDir(params.stage_dirs, r.meta.assembly.stage)}/gfastats/${r.meta.assembly.build}" }
    }
    quast: Channel<Record> {
        path "${params.stage_dirs.report}/quast"
    }
    // 10 - Report
    report: Channel<Record> {
        path "${params.stage_dirs.report}"
    }
}
