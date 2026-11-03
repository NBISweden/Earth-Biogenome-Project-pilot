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
include { INSPECT_DATA                                      } from "./subworkflows/local/01_inspect_data/main"
// Stage pipelines. Each one runs the assembly QC on the assemblies it produces.
include { workflow as ASSEMBLE          } from "./pipelines/assemble"
include { workflow as DECONTAMINATE     } from "./pipelines/decontaminate"
include { workflow as PURGE_DUPLICATES  } from "./pipelines/purge_duplicates"
include { workflow as DVPOLISH          } from "./pipelines/polish"
include { workflow as SCAFFOLD          } from "./pipelines/scaffold"
include { workflow as SCAFFOLD_CURATION } from "./pipelines/curate"
include { workflow as ALIGN_RNASEQ      } from "./pipelines/align_rnaseq"
// Comparison of the assemblies across the stages
include { workflow as COMPARE_STAGES    } from "./pipelines/compare_stages"
// Reports
include { REPORT_VERSIONS               } from "./subworkflows/local/version_report/main"
include { ASSEMBLY_REPORT               } from "./subworkflows/local/assembly_report/main"
include { Assembly; Sample; KmerParams; BuscoParams; FcsParams } from "./types/records"

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
    chunk_size : MemoryUnit = 100.MB  // DVPOLISH
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
def setAssemblyStage(a: Assembly, stage: String) {
    def build = "${a.meta.assembly.assembler}-${stage}-${a.meta.assembly.id}"
    return a + record(id: build, meta: ( a.meta + record(assembly: a.meta.assembly + [ stage: stage, build: build ]) ) as Record)
}

// Folder of a sample or an assembly build in the folder of a stage
def subjectDir(stage_dir: String, meta: Record) -> String {
    return "${stage_dir}/${meta.assembly?.build ?: meta.id}"
}

// Folder of a stage. A stage_dirs param that names only some stages keeps the default folders of the others
def stageDir(stage_dirs: Map<String,String>, key: String) -> String {
    def defaults = [
        inspect: '01_read_inspection', preprocess: '02_read_preprocessing', assembly: '03_assembly',
        screen: '04_contamination_screen', purge: '05_duplicate_purging', polish: '06_error_polishing',
        scaffold: '07_scaffolding', curate: '08_curation', alignRNA: '09_RNA_alignment', report: '10_report'
    ]
    return stage_dirs[key] ?: defaults[key]
}

// Directory of the stage that produced an assembly
def assemblyStageDir(stage_dirs: Map<String,String>, stage: String) -> String {
    def dirs = [ raw: 'assembly', decontaminated: 'screen', purged: 'purge', polished: 'polish', scaffolded: 'scaffold', curated: 'curate' ]
    return stageDir(stage_dirs, dirs[stage])
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
    // Hi-C reads are optional. A run goes as far as it can without them, but fails at the start if a step that
    // will run needs them
    input = PREPARE_INPUT(
        params.input,
        params.busco.lineages,
        [
            assemble    : enable.assemble,
            assemble_hic: enable.assemble && params.hifiasm_with_hic,
            scaffold    : enable.scaffold,
            curate      : enable.curate
        ]
    )

    // Create Cram and trim
    convert = CONVERT_FASTQ_CRAM( input.samples, params.hic_type.startsWith('arima'), params.hic_type )

    // Build necessary databases
    // TODO: Migrate back to Meryl. Genome inspection missing KATGC and SMUDGEPLOT for meryldb
    hifi_reads = input.samples.flatMap { s ->
        ( s.hifi ?: [] ).collect { f -> record(id: s.id, meta: s.meta + record(single_end: true), reads: [ f ]) }.toList()
    }
    fastk_hifi = BUILD_FASTK_HIFI_DATABASE( hifi_reads, 'hifi', fastk_kmer_size )
    fastk_hic  = BUILD_FASTK_HIC_DATABASE( convert.fastq, 'hic', fastk_kmer_size )
    meryl_hifi = BUILD_MERYL_HIFI_DATABASE( hifi_reads, 'hifi', meryl_kmer_size )
    BUILD_MERYL_HIC_DATABASE( convert.fastq, 'hic', meryl_kmer_size )

    // One wide record per sample. Stages take the fields they need. The HiFi reads are required, so every
    // sample has its HiFi databases.
    with_hifi_databases = input.samples
        .join(fastk_hifi.database.map { d -> record(id: d.id, fastk_hifi_hist: d.hist, fastk_hifi_ktab: d.ktab, fastk_hifi_prof: d.prof) }, by: 'id')
        .join(meryl_hifi.uniondb.map { m -> record(id: m.id, meryl_hifi_db: m.meryl_db) }, by: 'id')
    // The stages work on the converted Hi-C reads
    samples_ready = with_hifi_databases
        .join(convert.hic, by: 'id')
        .map { s -> s + record(hic: s.hic_fastq) }

    // Data inspection. A sample without Hi-C has no Hi-C database.
    inspect_samples = with_hifi_databases.join(
        fastk_hic.database.map { d -> record(id: d.id, fastk_hic_hist: d.hist, fastk_hic_ktab: d.ktab, fastk_hic_prof: d.prof) },
        by: 'id',
        remainder: true
    )
    inspect = INSPECT_DATA(
        inspect_samples.filter { _s -> enable.inspect },
        genescopefk_kmer_size
    )

    // Assemblies supplied as input, with the data of their sample
    user_assemblies = input.assemblies.join(
        samples_ready.map { s ->
            record(
                sample_id: s.id, hifi: s.hifi, hifi_merged: s.hifi_merged, hic: s.hic, rnaseq: s.rnaseq,
                fastk_hifi_hist: s.fastk_hifi_hist, fastk_hifi_ktab: s.fastk_hifi_ktab, fastk_hifi_prof: s.fastk_hifi_prof,
                meryl_hifi_db: s.meryl_hifi_db
            )
        },
        by: 'sample_id'
    )

    // The options of the assembly QC, passed to every stage
    def qc = record(busco: params.busco, fastk: params.fastk, meryl: params.meryl, use_phased: params.use_phased)

    // Assemble
    assemble = ASSEMBLE(
        record(
            samples: samples_ready.filter { _s -> enable.assemble },
            provided: user_assemblies.filter { a -> a.meta.assembly.stage == 'raw' },
            oatkdb: params.oatkdb,
            nuclear_assembly_mode: params.nuclear_assembly_mode,
            // The mode check of the pipeline must only apply when the stage runs
            organelle_assembly_mode: enable.assemble ? params.organelle_assembly_mode : 'none',
            hifiasm: params.hifiasm,
            hifiasm_with_hic: params.hifiasm_with_hic,
            qc: qc
        )
    )

    // Contamination screen
    to_screen = assemble.assemblies.map { a -> setAssemblyStage(a, 'decontaminated') }
    screen = DECONTAMINATE(
        record(
            fetch_database: enable.screen,
            assemblies: to_screen.filter { _a -> enable.screen },
            provided: user_assemblies.filter { a -> a.meta.assembly.stage == 'decontaminated' },
            fcs: params.fcs,
            qc: qc
        )
    )
    cleaned_assemblies = to_screen.filter { _a -> !enable.screen }.mix(screen.assemblies)

    // Purge duplicates
    to_purge = cleaned_assemblies.map { a -> setAssemblyStage(a, 'purged') }
    purge = PURGE_DUPLICATES(
        record(
            assemblies: to_purge.filter { _a -> enable.purge },
            provided: user_assemblies.filter { a -> a.meta.assembly.stage == 'purged' },
            qc: qc
        )
    )
    purged_assemblies = to_purge.filter { _a -> !enable.purge }.mix(purge.assemblies)

    // Polish
    to_polish = purged_assemblies.map { a -> setAssemblyStage(a, 'polished') }
    polish = DVPOLISH(
        record(
            assemblies: to_polish.filter { _a -> enable.polish },
            provided: user_assemblies.filter { a -> a.meta.assembly.stage == 'polished' },
            chunk_size: params.chunk_size,
            qc: qc
        )
    )
    polished_assemblies = to_polish.filter { _a -> !enable.polish }.mix(polish.assemblies)

    // Scaffold
    to_scaffold = polished_assemblies.map { a -> setAssemblyStage(a, 'scaffolded') }
    scaffold = SCAFFOLD(
        record(
            assemblies: to_scaffold.filter { _a -> enable.scaffold },
            provided: user_assemblies.filter { a -> a.meta.assembly.stage == 'scaffolded' },
            qc: qc,
            hic_type: params.hic_type,
            hic_bam_sorting: params.hic_bam_sorting
        )
    )
    scaffolded_assemblies = to_scaffold.filter { _a -> !enable.scaffold }.mix(scaffold.assemblies)

    // Curate
    to_curate = scaffolded_assemblies.map { a -> setAssemblyStage(a, 'curated') }
    curate = SCAFFOLD_CURATION(
        record(
            assemblies: to_curate.filter { _a -> enable.curate },
            provided: user_assemblies.filter { a -> a.meta.assembly.stage == 'curated' },
            qc: qc,
            hic_map_sort_by: params.hic_map_sort_by,
            hic_map_qv: params.hic_map_qv,
            cooler_bin_size: params.cooler_bin_size,
            hifi_coverage_cap: params.hifi_coverage_cap,
            telomer_motif: params.telomer_motif
        )
    )

    // Compare the assemblies across the stages: one QUAST run per sample. It waits for the last stage.
    stage_comparison = COMPARE_STAGES(
        record(
            assemblies: assemble.assemblies
                .mix(screen.assemblies)
                .mix(purge.assemblies)
                .mix(polish.assemblies)
                .mix(scaffold.assemblies)
                .mix(user_assemblies.filter { a -> a.meta.assembly.stage == 'curated' }),
            reference: params.reference
        )
    )

    // Align RNAseq
    rnaseq = ALIGN_RNASEQ(
        record(assemblies: user_assemblies.filter { _a -> enable.alignRNA })
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
        .mix(curate.logs)
        .mix(stage_comparison.logs)
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
    fcsgx_reports     = screen.fcsgx_reports
    fcsgx_clean       = screen.fcsgx_clean
    // 05 - Duplicate purging
    purge_dups        = purge.purge_dups
    // 06 - Error polishing
    dvpolish          = polish.dvpolish
    dvpolish_merqury_input    = polish.merqury_input
    dvpolish_merqury_polished = polish.merqury_polished
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
    // Assembly QC of each stage
    merquryfk         = assemble.merquryfk.mix(screen.merquryfk).mix(purge.merquryfk).mix(polish.merquryfk).mix(scaffold.merquryfk).mix(curate.merquryfk)
    merqury           = assemble.merqury.mix(screen.merqury).mix(purge.merqury).mix(polish.merqury).mix(scaffold.merqury).mix(curate.merqury)
    busco             = assemble.busco.mix(screen.busco).mix(purge.busco).mix(polish.busco).mix(scaffold.busco).mix(curate.busco)
    gfastats          = assemble.gfastats.mix(screen.gfastats).mix(purge.gfastats).mix(polish.gfastats).mix(scaffold.gfastats).mix(curate.gfastats)
    // 10 - Report
    quast             = stage_comparison.quast
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
        path { r -> "${stageDir(params.stage_dirs, 'inspect')}/${r.id}/goat" }
    }
    seqkit_hifi_stats: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'inspect')}/${r.meta.id}/seqkit_hifi_stats" }
    }
    seqkit_hic_stats: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'inspect')}/${r.meta.id}/seqkit_hic_stats" }
    }
    fastqc_hic: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'inspect')}/${r.meta.id}/fastqc_hic" }
    }
    fastk_runs: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'inspect')}/${r.meta.id}/fastk" }
    }
    fastk_merged: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'inspect')}/${r.meta.id}/fastk" }
    }
    genescopefk: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'inspect')}/${r.meta.id}/genescopefk" }
    }
    smudgeplot: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'inspect')}/${r.meta.id}/smudgeplot" }
    }
    katgc: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'inspect')}/${r.meta.id}/katgc" }
    }
    kat_comp: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'inspect')}/${r.meta.id}/kat_comp" }
    }
    dtol_search: Channel<Record> {
        path { r -> r.json_file >> ("${stageDir(params.stage_dirs, 'inspect')}/${r.id}/dtol_search/") }
    }
    // 02 - Read preprocessing
    converted_bam_to_fasta: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'preprocess')}/${r.meta.id}/converted_bam_to_fasta" }
    }
    hic_cram: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'preprocess')}/${r.meta.id}/hic_cram" }
    }
    hic_crai: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'preprocess')}/${r.meta.id}/hic_cram" }
    }
    // 03 - Assembly
    hifiasm: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'assembly'), r.meta)}/hifiasm" }
    }
    gfa2fa: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'assembly'), r.meta)}/gfa2fa" }
    }
    mitohifi_references: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'assembly')}/${r.meta.id}/organelle/mitohifi/references" }
    }
    mitohifi: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'assembly')}/${params.organelle_assembly_mode == 'contigs' ? r.meta.assembly.build : r.meta.id}/organelle/mitohifi" }
    }
    oatk: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'assembly'), r.meta)}/organelle/oatk" }
    }
    dotplots: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'assembly'), r.meta)}/organelle/dotplots" }
    }
    // 04 - Contamination screen
    fcsgx_reports: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'screen'), r.meta)}/fcsgx" }
    }
    fcsgx_clean: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'screen'), r.meta)}/fcsgx" }
    }
    // 05 - Duplicate purging
    purge_dups: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'purge'), r.meta)}/purge_dups" }
    }
    // 06 - Error polishing
    dvpolish: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'polish'), r.meta)}/dvpolish" }
    }
    dvpolish_merqury_input: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'polish'), r.meta)}/dvpolish/merqury_input" }
    }
    dvpolish_merqury_polished: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'polish'), r.meta)}/dvpolish/merqury_polished" }
    }
    // 07 - Scaffolding
    pairtools: Channel<Record> {
        path { r -> r.stat >> ("${subjectDir(stageDir(params.stage_dirs, 'scaffold'), r.meta)}/pairtools/") }
    }
    yahs: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'scaffold'), r.meta)}/yahs" }
    }
    merged_haps: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'scaffold'), r.meta)}/merged_haps" }
    }
    // 08 - Curation
    higlass_mcool: Channel<Record> {
        path { r -> r.mcool >> ("${subjectDir(stageDir(params.stage_dirs, 'curate'), r.meta)}/higlass/") }
    }
    higlass_gaps: Channel<Record> {
        path { r -> r.beddb >> ("${subjectDir(stageDir(params.stage_dirs, 'curate'), r.meta)}/higlass/") }
    }
    higlass_telomers: Channel<Record> {
        path { r -> r.bw >> ("${subjectDir(stageDir(params.stage_dirs, 'curate'), r.meta)}/higlass/") }
    }
    pretext: Channel<Record> {
        path { r -> r.hitile >> ("${subjectDir(stageDir(params.stage_dirs, 'curate'), r.meta)}/pretext/") }
    }
    // 09 - RNA alignment
    star: Channel<Record> {
        path { r -> "${subjectDir(stageDir(params.stage_dirs, 'alignRNA'), r.meta)}/star" }
    }
    // Assembly QC, next to the assembly it evaluates
    merquryfk: Channel<Record> {
        path { r -> "${assemblyStageDir(params.stage_dirs, r.meta.assembly.stage)}/${r.meta.assembly.build}/qc/merquryfk" }
    }
    merqury: Channel<Record> {
        path { r -> "${assemblyStageDir(params.stage_dirs, r.meta.assembly.stage)}/${r.meta.assembly.build}/qc/merqury" }
    }
    busco: Channel<Record> {
        path { r -> "${assemblyStageDir(params.stage_dirs, r.meta.assembly.stage)}/${r.meta.assembly.build}/qc/busco" }
    }
    gfastats: Channel<Record> {
        path { r -> "${assemblyStageDir(params.stage_dirs, r.meta.assembly.stage)}/${r.meta.assembly.build}/qc/gfastats" }
    }
    // 10 - Report
    quast: Channel<Record> {
        path { r -> "${stageDir(params.stage_dirs, 'report')}/quast/${r.meta.id}" }
    }
    report: Channel<Record> {
        path "${stageDir(params.stage_dirs, 'report')}"
    }
}
