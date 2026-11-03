nextflow.enable.types = true

process SAMTOOLS_VIEW {
    tag "${id}"
    label 'process_low'

    ext [
        args: { samtools_args ?: '' },
        args2: { samtools_args2 ?: '' },
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
        input: Path,
        index: Path?,
        fasta: Path?,
        fai: Path?,
        qname: Path?,
        bed: Path?,
        index_format: String?,
        samtools_args: String?,
        samtools_args2: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        bam: file("${task.ext.prefix}.bam", optional: true),
        cram: file("${task.ext.prefix}.cram", optional: true),
        sam: file("${task.ext.prefix}.sam", optional: true),
        bai: file("${task.ext.prefix}.*.bai", optional: true),
        csi: file("${task.ext.prefix}.*.csi", optional: true),
        crai: file("${task.ext.prefix}.*.crai", optional: true),
        unselected: file("${task.ext.prefix}.unselected.*", optional: true)
    )

    topic:
    tuple(task.process, 'samtools', eval('samtools version | sed "1!d;s/.* //"')) >> 'versions'

    script:
    def args = task.ext.args
    def prefix = task.ext.prefix
    def reference = fasta ? "--reference ${fasta}" : ""
    def file_type = args.contains("--output-fmt sam")
        ? "sam"
        : args.contains("--output-fmt bam")
            ? "bam"
            : args.contains("--output-fmt cram")
                ? "cram"
                : input.extension

    def output_file = index_format ? "${prefix}.${file_type}##idx##${prefix}.${file_type}.${index_format} --write-index" : "${prefix}.${file_type}"
    // Can't choose index type of unselected file
    def readnames = qname ? "--qname-file ${qname} --output-unselected ${prefix}.unselected.${file_type}" : ""
    def bedfile = bed ? "-L ${bed}" : ""

    if ("${input}" == "${prefix}.${file_type}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    if (index_format) {
        if (!['bai', 'csi', 'crai'].contains(index_format)) {
            error("Index format not one of bai, csi, crai.")
        }
        else if (file_type == "sam") {
            error("Indexing not compatible with SAM output")
        }
    }
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        view \\
        --threads ${task.cpus - 1} \\
        ${reference} \\
        ${readnames} \\
        ${bedfile} \\
        ${args} \\
        -o ${output_file} \\
        ${input} \\
        ${task.ext.args2}
    """

    stub:
    def args = task.ext.args
    def prefix = task.ext.prefix
    def file_type = args.contains("--output-fmt sam")
        ? "sam"
        : args.contains("--output-fmt bam")
            ? "bam"
            : args.contains("--output-fmt cram")
                ? "cram"
                : input.extension
    def default_index_format = file_type == "bam"
        ? "csi"
        : file_type == "cram" ? "crai" : ""
    def index_cmd = index_format ? "touch ${prefix}.${file_type}.${index_format}" : args.contains("--write-index") ? "touch ${prefix}.${file_type}.${default_index_format}" : ""
    def unselected = qname ? "touch ${prefix}.unselected.${file_type}" : ""

    if ("${input}" == "${prefix}.${file_type}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.${file_type}
    ${index_cmd}
    ${unselected}
    """
}
