nextflow.enable.types = true

process DVPOLISH_PBMM2_INDEX {
    tag "${id}"
    label 'process_medium'

    ext [
        args: { pbmm2_args ?: '' }
    ]

    // Note: the versions here need to match the versions used in pbmm2/align
    conda 'bioconda::pbmm2=1.13.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmm2:1.13.1--h9ee0642_0' :
        'biocontainers/pbmm2:1.13.1--h9ee0642_0' }"

    input:
    record(
        id: String,
        meta: Record?,
        fasta: Path,
        pbmm2_args: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        index: file("*.mmi")
    )

    topic:
    tuple(task.process, 'pbmm2', eval("pbmm2 --version 2>&1 | sed '1!d; s/pbmm2 //'")) >> 'versions'

    script:
    """
    pbmm2 index\\
        -j $task.cpus \\
        ${task.ext.args} \\
        $fasta \\
        ${fasta.baseName}.mmi
    """

    stub:
    """
    touch ${fasta.baseName}.mmi
    """
}
