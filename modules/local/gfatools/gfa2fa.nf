nextflow.enable.types = true

process GFATOOLS_GFA2FA {
    tag "${id}"
    label 'process_single'

    ext [
        args: { gfatools_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    conda "bioconda::gfatools=0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfatools:0.5--he4a0461_4':
        'biocontainers/gfatools:0.5--he4a0461_4' }"

    input:
    record(
        id: String,
        meta: Record?,
        gfa: Path,
        gfatools_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        fasta: file("*.fasta.gz")
    )

    topic:
    tuple(task.process, 'gfatools', eval("gfatools version | sed '1!d; s/.* //'")) >> 'versions'

    script:
    """
    gfatools \\
        gfa2fa \\
        ${task.ext.args} \\
        $gfa \\
        | gzip -c > ${task.ext.prefix}.fasta.gz
    """

    stub:
    """
    echo "" | gzip > ${task.ext.prefix}.fasta.gz
    """
}
