process FCSGX_FETCHDB {
    tag "$manifest.name"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::ncbi-fcs-gx=0.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.4--h4ac6f70_1':
        'biocontainers/ncbi-fcs-gx:0.5.4--h4ac6f70_1' }"

    input:
    val manifest // URL of manifest. Cannot stage locally.

    output:
    path "$prefix"      , emit: database
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fcs_gx'), val('0.5.4'), emit: versions_fcsgx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "gxdb_$manifest.baseName"
    """
    sync_files.py \\
        get \\
        --mft "${manifest.toUriString()}" \\
        --dir "$prefix"
    """
}
