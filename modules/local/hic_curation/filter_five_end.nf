nextflow.enable.types = true

process FILTER_FIVE_END {
    tag "${id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/samtools_perl:e2ce3f7265547d1f"

    ext [
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        bam: Path,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        bam: file("*.bam")
    )

    topic:
    tuple(task.process, 'samtools', eval('samtools --version | sed "1!d; s/samtools //"')) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${id}"

    """
    samtools view -h -@ $task.cpus ${bam} | \\
    filter_five_end.pl | \\
    samtools view -@ $task.cpus -Sb - > ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}.bam
    """
}
