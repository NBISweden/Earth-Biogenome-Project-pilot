nextflow.enable.types = true

process BIOBAMBAM_BAMMARKDUPLICATES2 {
    tag "${id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biobambam:2.0.185--h85de650_1'
        : 'biocontainers/biobambam:2.0.185--h85de650_1'}"

    ext [
        args: { biobambam_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        bam: Path,
        biobambam_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        bam: file("*.bam"),
        metrics: file("*.metrics.txt")
    )

    topic:
    tuple(task.process, 'biobambam', eval("bammarkduplicates2 --version |& sed '1!d; s/.*version //; s/.\$//'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    bammarkduplicates2 \\
        ${args} \\
        I=${bam} \\
        O=${prefix}.bam \\
        M=${prefix}.metrics.txt \\
        tmpfile=${prefix} \\
        markthreads=${task.cpus}
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.metrics.txt
    """
}
