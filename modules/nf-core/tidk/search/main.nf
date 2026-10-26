nextflow.enable.types = true

process TIDK_SEARCH {
    tag "${id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tidk:0.2.41--hdbdd923_0':
        'biocontainers/tidk:0.2.41--hdbdd923_0' }"

    ext [
        args: { tidk_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        fasta: Path,
        string: String,
        tidk_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        tsv: file("*.tsv", optional: true),
        bedgraph: file("*.bedgraph", optional: true)
    )

    topic:
    tuple(task.process, 'tidk', eval("tidk --version | sed 's/tidk //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    tidk \\
        search \\
        --string $string \\
        --output $prefix \\
        --dir tidk \\
        $args \\
        $fasta

    mv \\
        tidk/${prefix}_telomeric_repeat_windows.tsv \\
        ${prefix}.tsv \\
        || echo "TSV file was not produced"

    mv \\
        tidk/${prefix}_telomeric_repeat_windows.bedgraph \\
        ${prefix}.bedgraph \\
        || echo "BEDGRAPH file was not produced"
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    def extension = args.contains("--extension bedgraph") ? 'bedgraph' : 'tsv'
    """
    touch ${prefix}.${extension}
    """
}
