process SMUDGEPLOT {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    // TODO: smudgeplot has an error rate calc bug (https://github.com/KamilSJaron/smudgeplot/issues/227), update container once fixed
    container 'ghcr.io/nbisweden/smudgeplot:1.0'

    input:
    tuple val(meta), path(fastk_ktab)

    output:
    tuple val(meta), path("*.smudge_report.tsv")     , emit: smudgeplot_report_tsv
    tuple val(meta), path("*_centralities.txt")      , emit: centralities_txt
    tuple val(meta), path("*_annotated_smu.txt")     , emit: kmer_pairs_txt
    tuple val(meta), path("*_centralities.png")      , emit: centralities_png      , optional: true
    tuple val(meta), path("*_centralities.pdf")      , emit: centralities_pdf      , optional: true
    tuple val(meta), path("*_smudgeplot.png")        , emit: smudgeplot_png        , optional: true
    tuple val(meta), path("*_smudgeplot.pdf")        , emit: smudgeplot_pdf        , optional: true
    tuple val(meta), path("*_smudgeplot_log10.png")  , emit: smudgeplot_log10_png  , optional: true
    tuple val(meta), path("*_smudgeplot_log10.pdf")  , emit: smudgeplot_log10_pdf  , optional: true
    tuple val(meta), path("*_smudgeplot_report.json"), emit: smudgeplot_report_json, optional: true
    path "versions.yml"                              , emit: versions

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
    def FASTK_VERSION = '1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    export MPLCONFIGDIR=tmp
    smudgeplot \\
        hetmers \\
        $args \\
        -t ${task.cpus} \\
        -o $prefix \\
        ${fastk_ktab.find{ it.toString().endsWith(".ktab") }}
    smudgeplot \\
        all \\
        $args2 \\
        -o $prefix \\
        ${prefix}.smu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smudgeplot: \$(smudgeplot --version 2>&1 | sed 's/.* v//')
        fastk: $FASTK_VERSION
    END_VERSIONS
    """
}
