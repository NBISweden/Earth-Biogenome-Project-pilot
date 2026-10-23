nextflow.enable.types = true

process MERYL_COUNT {
    tag "${id}"
    label 'process_high'

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
        reads: List<Path>,
        kvalue: Integer,
        meryl_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        meryl_db: files("*.meryl")
    )

    topic:
    tuple(task.process, 'meryl', eval("meryl --version |& sed 's/meryl //'")) >> 'versions'

    script:
    def reduced_mem = (task.memory.toGiga() * 0.9).intValue()
    """
    for READ in ${reads.join(' ')}; do
        meryl count \\
            k=${kvalue} \\
            threads=${task.cpus} \\
            memory=${reduced_mem} \\
            ${task.ext.args} \\
            \$READ \\
            output ${task.ext.prefix}.\${READ%.f*}.meryl
    done
    """

    stub:
    """
    for READ in ${reads.join(' ')}; do
        touch ${task.ext.prefix}.\${READ%.f*}.meryl
    done
    """
}
