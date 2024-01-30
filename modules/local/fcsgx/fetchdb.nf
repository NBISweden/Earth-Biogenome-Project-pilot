process FCSGX_FETCHDB {
    tag "$manifest.name"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.0--h4ac6f70_3':
        'biocontainers/ncbi-fcs-gx:0.5.0--h4ac6f70_3' }"

    input:
    val manifest // URL of manifest. Cannot stage locally.

    output:
    path "$prefix"      , emit: database
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "gxdb_$manifest.baseName"
    def VERSION = '0.5.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    sync_files.py \\
        get \\
        --mft "${manifest.toUriString()}" \\
        --dir "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs_gx: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def VERSION = '0.5.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs_gx: $VERSION
    END_VERSIONS
    """
}
