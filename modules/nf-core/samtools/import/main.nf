nextflow.enable.types = true

process SAMTOOLS_IMPORT {
    tag "${id}"
    label 'process_low'

    ext [
        args: { samtools_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    record(
        id: String,
        meta: Record?,
        single_end: Boolean,
        reads: List<Path>,
        samtools_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        sam: file("*.sam", optional: true),
        bam: file("*.bam", optional: true),
        cram: file("*.cram", optional: true)
    )

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args
    def prefix = task.ext.prefix
    def suffix = args.contains("--output-fmt sam")
        ? "sam"
        : args.contains("--output-fmt bam")
            ? "bam"
            : args.contains("--output-fmt cram")
                ? "cram"
                : "bam"
    def input = single_end
        ? "-0 ${reads.join(' -0 ')}"
        : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    samtools \\
        import \\
        ${input} \\
        ${args} \\
        -@ ${task.cpus} \\
        -o ${prefix}.${suffix}
    """

    stub:
    def suffix = task.ext.args.contains("--output-fmt sam")
        ? "sam"
        : task.ext.args.contains("--output-fmt cram")
            ? "cram"
            : "bam"
    """
    touch ${task.ext.prefix}.${suffix}
    """
}
