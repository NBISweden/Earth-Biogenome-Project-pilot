nextflow.enable.types = true

process SEQKIT_STATS {
    tag "${id}"
    label 'process_low'

    ext [
        args: { seqkit_args ?: '--all' },
        prefix: { output_prefix ?: "${id}" }
    ]

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    record(
        id: String,
        meta: Record?,
        reads: List<Path>,
        seqkit_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        stats: file("*.tsv")
    )

    topic:
    tuple(task.process, 'seqkit', eval("seqkit version | sed 's/^.*v//'")) >> 'versions'

    script:
    """
    seqkit stats \\
        --tabular \\
        --threads ${task.cpus} \\
        ${task.ext.args} \\
        ${reads.join(' ')} > '${task.ext.prefix}.tsv'
    """

    stub:
    """
    touch ${task.ext.prefix}.tsv
    """
}
