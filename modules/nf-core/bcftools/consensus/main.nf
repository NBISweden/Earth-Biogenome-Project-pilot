nextflow.enable.types = true

process BCFTOOLS_CONSENSUS {
    tag "${id}"
    label 'process_medium'

    ext [
        args: { bcftools_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    record(
        id: String,
        meta: Record?,
        vcf: Path,
        tbi: Path,
        fasta: Path,
        mask: Path?,
        bcftools_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        fasta: file('*.fa')
    )

    topic:
    tuple(task.process, 'bcftools', eval("bcftools --version | sed '1!d; s/^.*bcftools //'")) >> 'versions'

    script:
    def masking = mask ? "-m ${mask}" : ""
    """
    cat ${fasta} \\
        | bcftools \\
            consensus \\
            ${vcf} \\
            ${task.ext.args} \\
            ${masking} \\
            > ${task.ext.prefix}.fa
    """

    stub:
    """
    touch ${task.ext.prefix}.fa
    """
}
