process FASTK_FASTK {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.hist")                      , emit: hist
    tuple val(meta), path("*.ktab*", hidden: true)       , emit: ktab, optional: true
    tuple val(meta), path("*.{prof,pidx}*", hidden: true), emit: prof, optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FASTK_FASTK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '0e24fb45b71c4e14382ae1e1bc063bf66ea4e112' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir tmp_uncompressed
    
    FastK \\
        $args \\
        -T$task.cpus \\
        -Ptmp_uncompressed \\
        -M${task.memory.toGiga()} \\
        -N${prefix}_fk \\
        $reads

    rm -rf tmp_uncompressed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
    END_VERSIONS
    """
}
