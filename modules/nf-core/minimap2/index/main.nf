nextflow.enable.types = true

process MINIMAP2_INDEX {
    tag "${id}"
    label 'process_low'

    // Note: the versions here need to match the versions used in minimap2/align
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.29--h577a1d6_0' :
        'biocontainers/minimap2:2.29--h577a1d6_0' }"

    ext [
        args: { minimap2_args ?: '' }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        fasta: Path,
        minimap2_args: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        index: file("*.mmi")
    )

    topic:
    tuple(task.process, "minimap2", eval("minimap2 --version")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${fasta.baseName}.mmi \\
        $args \\
        $fasta
    """

    stub:
    """
    touch ${fasta.baseName}.mmi
    """
}
