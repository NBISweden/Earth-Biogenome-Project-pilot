process FASTK_HISTEX {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container 'oras://community.wave.seqera.io/library/fastk:1.1.0--c033613137caa2b9'

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("*.hist"), emit: hist
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FASTK_HISTEX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    Histex \\
        $args \\
        $histogram \\
        > ${prefix}.hist

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
    END_VERSIONS
    """
}
