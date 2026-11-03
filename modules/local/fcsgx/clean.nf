nextflow.enable.types = true

process FCSGX_CLEAN {
    tag "${id}"
    label 'process_low'

    ext [
        prefix: { output_prefix ?: "${id}" }
    ]

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::ncbi-fcs-gx=0.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.4--h4ac6f70_1':
        'biocontainers/ncbi-fcs-gx:0.5.4--h4ac6f70_1' }"

    input:
    record(
        id: String,
        meta: Record?,
        assembly: Path,
        action_report: Path,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        clean_fasta: file("*.clean.fasta"),
        contaminants_fasta: file("*.contaminants.fasta")
    )

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'fcs_gx', '0.5.4') >> 'versions'

    script:
    def prefix = task.ext.prefix
    def infile = assembly.name.endsWith('.gz') ? "<( gzip -dc $assembly )" : "${assembly}"
    """
    gx \\
        clean-genome \\
        --input $infile \\
        --action-report $action_report \\
        --output ${prefix}.clean.fasta \\
        --contam-fasta-out ${prefix}.contaminants.fasta
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}.clean.fasta ${prefix}.contaminants.fasta
    """
}
