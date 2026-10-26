nextflow.enable.types = true

process TIDK_PLOT {
    tag "${id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tidk:0.2.41--hdbdd923_0':
        'biocontainers/tidk:0.2.41--hdbdd923_0' }"

    ext [
        args: { tidk_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        tsv: Path,
        tidk_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        svg: file("*.svg")
    )

    topic:
    tuple(task.process, 'tidk', eval("tidk --version | sed 's/tidk //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    tidk \\
        plot \\
        --output $prefix \\
        $args \\
        --tsv "$tsv"
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}.svg
    """
}
