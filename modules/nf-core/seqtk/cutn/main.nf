nextflow.enable.types = true

process SEQTK_CUTN {
    tag "${id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'biocontainers/seqtk:1.4--he4a0461_1' }"

    ext [
        args: { seqtk_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        fasta: Path,
        seqtk_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        bed: file("*.bed")
    )

    topic:
    tuple(task.process, 'seqtk', eval("seqtk 2>&1 | sed -n 's/^Version: //p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"

    """
    seqtk \\
        cutN \\
        $args \\
        -g $fasta \\
        > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"

    """
    touch ${prefix}.bed
    """

}
