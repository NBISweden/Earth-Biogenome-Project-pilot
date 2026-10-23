nextflow.enable.types = true

process SAMTOOLS_INDEX {
    tag "${id}"
    label 'process_low'

    ext [
        args: { samtools_args ?: '' }
    ]

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    record(
        id: String,
        meta: Record?,
        input: Path,
        samtools_args: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        index: file("*.{bai,csi,crai}")
    )

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    """
    samtools \\
        index \\
        -@ ${task.cpus} \\
        ${task.ext.args} \\
        ${input}
    """

    stub:
    def extension = input.extension == 'cram'
        ? "crai"
        : task.ext.args.contains("-c") ? "csi" : "bai"
    """
    touch ${input}.${extension}
    """
}
