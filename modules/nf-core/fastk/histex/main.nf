nextflow.enable.types = true

process FASTK_HISTEX {
    tag "${id}"
    label 'process_low'

    ext [
        args: { histex_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    //conda "${moduleDir}/environment.yml"
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    input:
    record(
        id: String,
        meta: Record?,
        histogram: Path,
        histex_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        hist: file("*.hist")
    )

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'fastk', '1.2') >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FASTK_HISTEX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    Histex \\
        ${task.ext.args} \\
        $histogram \\
        > ${task.ext.prefix}.hist
    """

    stub:
    """
    touch ${task.ext.prefix}.hist
    """
}
