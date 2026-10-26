nextflow.enable.types = true

process COOLER_ZOOMIFY {
    tag "${id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.10.4--pyhdfd78af_0' :
        'biocontainers/cooler:0.10.4--pyhdfd78af_0' }"

    ext [
        args: { cooler_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        cool: Path,
        cooler_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        mcool: file("*.mcool")
    )

    topic:
    tuple(task.process, 'cooler', eval('cooler --version 2>&1 | sed "s/cooler, version //"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    cooler zoomify \\
        $args \\
        -n $task.cpus \\
        -o ${prefix}.mcool \\
        $cool
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}.mcool
    """
}
