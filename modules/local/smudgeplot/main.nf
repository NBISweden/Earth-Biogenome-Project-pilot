nextflow.enable.types = true

process SMUDGEPLOT {
    tag "${id}"
    label 'process_medium'

    ext [
        args: { smudgeplot_args ?: '' },
        args2: { smudgeplot_args2 ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container 'ghcr.io/nbisweden/smudgeplot:1.1'

    input:
    record(
        id: String,
        meta: Record?,
        fastk_ktab: List<Path>,
        smudgeplot_threshold: Integer,
        smudgeplot_args: String?,
        smudgeplot_args2: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        annotated_kmer_pairs_tsv: file("*.sma"),
        kmer_pairs_tsv: file("*.smu"),
        smudgeplot_report_tsv: file("*.smudge_report.tsv"),
        centralities: file("*_centralities.{png,pdf}"),
        centralities_txt: file("*_centralities.txt"),
        smudgeplot: file("*_smudgeplot.{png,pdf}"),
        smudgeplot_log10: file("*_smudgeplot_log10.{png,pdf}"),
        smudgeplot_report_json: file("*_smudgeplot_report.json", optional: true)
    )

    topic:
    tuple(task.process, 'smudgeplot', eval('smudgeplot --version 2>&1 | sed "s/.* v//"')) >> 'versions'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'fastk', '1.1') >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SMUDGEPLOT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix
    """
    export MPLCONFIGDIR=tmp
    smudgeplot \\
        hetmers \\
        ${task.ext.args} \\
        -t ${task.cpus} \\
        -o $prefix \\
        -L $smudgeplot_threshold \\
        ${fastk_ktab.find{ f -> f.name.endsWith(".ktab") }}
    smudgeplot \\
        all \\
        ${task.ext.args2} \\
        -o $prefix \\
        ${prefix}.smu
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}.sma ${prefix}.smu ${prefix}.smudge_report.tsv
    touch ${prefix}_centralities.png ${prefix}_centralities.txt
    touch ${prefix}_smudgeplot.png ${prefix}_smudgeplot_log10.png
    touch ${prefix}_smudgeplot_report.json
    """
}
