nextflow.enable.types = true

process GENESCOPEFK {
    tag "${id}"
    label 'process_low'

    ext [
        args: { genescopefk_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    input:
    record(
        id: String,
        meta: Record?,
        histogram: Path,
        genescopefk_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        linear_plot: file("${task.ext.prefix}_linear_plot.png"),
        log_plot: file("${task.ext.prefix}_log_plot.png"),
        model: file("${task.ext.prefix}_model.txt"),
        summary: file("${task.ext.prefix}_summary.txt"),
        transformed_linear_plot: file("${task.ext.prefix}_transformed_linear_plot.png"),
        transformed_log_plot: file("${task.ext.prefix}_transformed_log_plot.png"),
        log: file("${task.ext.prefix}_genescopefk.log"),
        kmer_cov: env('KMERCOV')
    )

    topic:
    tuple(task.process, 'genescopefk', "380815c420f50171f9234a0fd1ff426b39829b91") >> 'versions'
    tuple(task.process, 'R', eval("R --version | sed '1!d; s/.*version //; s/ .*//'")) >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GENESCOPEFK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix
    """
    #! /usr/bin/env bash

    GeneScopeFK.R \\
        ${task.ext.args} \\
        --input $histogram \\
        --output . \\
        --name_prefix ${prefix} > ${prefix}_genescopefk.log

    KMERCOV=\$( awk '/^kmercov/ {printf "%.2f", \$2}' *_model.txt )
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch "${prefix}_linear_plot.png"
    touch "${prefix}_log_plot.png"
    echo "kmercov 20.00" > "${prefix}_model.txt"
    touch "${prefix}_summary.txt"
    touch "${prefix}_transformed_linear_plot.png"
    touch "${prefix}_transformed_log_plot.png"
    touch "${prefix}_genescopefk.log"

    KMERCOV=\$( awk '/^kmercov/ {printf "%.2f", \$2}' *_model.txt )
    """
}
