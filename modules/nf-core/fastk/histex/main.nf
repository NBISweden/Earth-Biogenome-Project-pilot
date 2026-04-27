process FASTK_HISTEX {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    //conda "${moduleDir}/environment.yml"
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("*.hist"), emit: hist
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastk'), val('1.2'), emit: versions_fastk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FASTK_HISTEX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    Histex \\
        $args \\
        $histogram \\
        > ${prefix}.hist
    """

    stub:
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hist
    """
}
