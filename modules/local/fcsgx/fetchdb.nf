nextflow.enable.types = true

process FCSGX_FETCHDB {
    tag "${id}"
    label 'process_single'

    ext [
        prefix: { output_prefix ?: "gxdb_${manifest.tokenize('/').last().replaceAll(/\.[^.]*$/, '')}" }
    ]

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::ncbi-fcs-gx=0.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.4--h4ac6f70_1':
        'biocontainers/ncbi-fcs-gx:0.5.4--h4ac6f70_1' }"

    input:
    record(
        id: String,
        manifest: String, // URL of manifest. Cannot stage locally.
        output_prefix: String?
    )

    output:
    record(
        id: id,
        database: file("${task.ext.prefix}")
    )

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'fcs_gx', '0.5.4') >> 'versions'

    script:
    """
    sync_files.py \\
        get \\
        --mft "${manifest}" \\
        --dir "${task.ext.prefix}"
    """

    stub:
    """
    mkdir ${task.ext.prefix}
    """
}
