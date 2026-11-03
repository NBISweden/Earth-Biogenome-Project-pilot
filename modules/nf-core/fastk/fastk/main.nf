nextflow.enable.types = true

process FASTK_FASTK {
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
        reads: List<Path>,
        fastk_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        hist: file("*.hist"),
        log: file("*.log"),
        ktab: files("*.ktab*", hidden: true, optional: true),
        prof: files("*.{prof,pidx}*", hidden: true, optional: true)
    )

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'fastk', '1.2') >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FASTK_FASTK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix
    """
    mkdir tmp_uncompressed

    FastK \\
        ${task.ext.args} \\
        -T$task.cpus \\
        -Ptmp_uncompressed \\
        -M${task.memory.toGiga()} \\
        -N${prefix}_fk \\
        ${reads.join(' ')} \\
        1>${prefix}.fastK.log 2>&1

    find . -name '*.ktab*' -exec chmod a+r {} \\;

    rm -rf tmp_uncompressed
    """

    stub:
    def args       = task.ext.args
    def prefix     = task.ext.prefix
    def touch_ktab = args.contains('-t') ? "touch ${prefix}.ktab .${prefix}.ktab.1" : ''
    def touch_prof = args.contains('-p') ? "touch ${prefix}.prof .${prefix}.pidx.1" : ''
    """
    touch ${prefix}.hist
    $touch_ktab
    $touch_prof

    echo \\
    "FastK \\
        $args \\
        -T$task.cpus \\
        -M${task.memory.toGiga()} \\
        -N${prefix}_fk \\
        ${reads.join(' ')}" 1>${prefix}.fastK.log 2>&1
    """
}
