process FASTK_MERGE {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    //conda "${moduleDir}/environment.yml"
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    input:
    tuple val(meta), path(hist), path(ktab), path(prof)

    output:
    tuple val(meta), path("*.hist")                      , emit: hist
    tuple val(meta), path("*.ktab*", hidden: true)       , emit: ktab, optional: true
    tuple val(meta), path("*.{prof,pidx}*", hidden: true), emit: prof, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastk'), val('1.2'), emit: versions_fastk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FASTK_MERGE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    Fastmerge \\
        $args \\
        -T$task.cpus \\
        ${prefix} \\
        $hist
    """

    stub:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def touch_hist    = args.contains('-h') ? "touch ${prefix}.hist"                   : ''
    def touch_ktab    = args.contains('-t') ? "touch ${prefix}.ktab .${prefix}.ktab.1" : ''
    """
    ${touch_hist}
    ${touch_ktab}
    """
}
