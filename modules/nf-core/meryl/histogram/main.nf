nextflow.enable.types = true

process MERYL_HISTOGRAM {
    tag "${id}"
    label 'process_low'

    ext [
        args: { meryl_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meryl:1.4.1--h4ac6f70_0':
        'biocontainers/meryl:1.4.1--h4ac6f70_0' }"

    input:
    record(
        id: String,
        meta: Record?,
        meryl_db: Path,
        kvalue: Integer,
        meryl_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        hist: file("*.hist")
    )

    topic:
    tuple(task.process, 'meryl', eval("meryl --version |& sed 's/meryl //'")) >> 'versions'

    script:
    """
    meryl histogram \\
        k=$kvalue \\
        threads=$task.cpus \\
        memory=${task.memory.toGiga()} \\
        ${task.ext.args} \\
        $meryl_db > ${task.ext.prefix}.hist
    """

    stub:
    """
    touch ${task.ext.prefix}.hist
    """
}
