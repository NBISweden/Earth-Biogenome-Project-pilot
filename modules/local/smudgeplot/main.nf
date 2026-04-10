process SMUDGEPLOT {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container 'ghcr.io/nbisweden/smudgeplot:1.1'

    input:
    tuple val(meta), path(fastk_ktab)
    val smudgeplot_threshold

    output:
    tuple val(meta), path("*.sma")                       , emit: annotated_kmer_pairs_tsv
    tuple val(meta), path("*.smu")                       , emit: kmer_pairs_tsv
    tuple val(meta), path("*.smudge_report.tsv")         , emit: smudgeplot_report_tsv
    tuple val(meta), path("*_centralities.{png,pdf}")    , emit: centralities
    tuple val(meta), path("*_centralities.txt")          , emit: centralities_txt
    tuple val(meta), path("*_smudgeplot.{png,pdf}")      , emit: smudgeplot
    tuple val(meta), path("*_smudgeplot_log10.{png,pdf}"), emit: smudgeplot_log10
    tuple val(meta), path("*_smudgeplot_report.json")    , emit: smudgeplot_report_json, optional: true
    tuple val("${task.process}"), val('smudgeplot'), eval('smudgeplot --version 2>&1 | sed "s/.* v//"'), emit: versions_smudgeplot, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastk'), val('1.1'), emit: versions_fastk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SMUDGEPLOT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=tmp
    smudgeplot \\
        hetmers \\
        $args \\
        -t ${task.cpus} \\
        -o $prefix \\
        -L $smudgeplot_threshold \\
        ${fastk_ktab.find{ it.toString().endsWith(".ktab") }}
    smudgeplot \\
        all \\
        $args2 \\
        -o $prefix \\
        ${prefix}.smu
    """
}
