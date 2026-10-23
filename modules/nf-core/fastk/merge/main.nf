nextflow.enable.types = true

process FASTK_MERGE {
    tag "${id}"
    label 'process_medium'

    ext [
        args: { fastk_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    //conda "${moduleDir}/environment.yml"
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    input:
    record(
        id: String,
        meta: Record?,
        hist: List<Path>,
        ktab: List<Path>,
        prof: List<Path>,
        fastk_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        hist: file("*.hist"),
        ktab: files("*.ktab*", hidden: true, optional: true),
        prof: files("*.{prof,pidx}*", hidden: true, optional: true)
    )

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'fastk', '1.2') >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FASTK_MERGE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    Fastmerge \\
        ${task.ext.args} \\
        -T$task.cpus \\
        ${task.ext.prefix} \\
        ${hist.join(' ')}
    """

    stub:
    def args          = task.ext.args
    def prefix        = task.ext.prefix
    def touch_hist    = args.contains('-h') ? "touch ${prefix}.hist"                   : ''
    def touch_ktab    = args.contains('-t') ? "touch ${prefix}.ktab .${prefix}.ktab.1" : ''
    """
    ${touch_hist}
    ${touch_ktab}
    """
}
