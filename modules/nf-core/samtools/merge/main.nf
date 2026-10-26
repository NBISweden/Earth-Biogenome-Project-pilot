nextflow.enable.types = true

process SAMTOOLS_MERGE {
    tag "${id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    ext [
        args: { samtools_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        input_files: List<Path>,
        index_files: List<Path>,
        fasta: Path?,
        fai: Path?,
        gzi: Path?,
        samtools_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        bam: file("${task.ext.prefix}.bam", optional: true),
        cram: file("${task.ext.prefix}.cram", optional: true),
        index: file("*.{bai,crai,csi}", optional: true)
    )

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix
    def file_type = input_files[0].extension
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        merge \\
        --threads ${task.cpus - 1} \\
        ${args} \\
        ${reference} \\
        ${prefix}.${file_type} \\
        ${input_files.toSorted { a, b -> a.name <=> b.name }.join(' ')}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix
    def file_type = input_files[0].extension
    def index_type = file_type == "bam" ? "csi" : "crai"
    def index = args.contains("--write-index") ? "touch ${prefix}.${index_type}" : ""
    """
    touch ${prefix}.${file_type}
    ${index}
    """
}
