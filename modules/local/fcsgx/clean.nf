process FCSGX_CLEAN {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::ncbi-fcs-gx=0.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.4--h4ac6f70_1':
        'biocontainers/ncbi-fcs-gx:0.5.4--h4ac6f70_1' }"

    input:
    tuple val(meta), path(assembly), path(action_report)

    output:
    tuple val(meta), path("*.clean.fasta")       , emit: clean_fasta
    tuple val(meta), path("*.contaminants.fasta"), emit: contaminants_fasta
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fcs_gx'), val('0.5.4'), emit: versions_fcsgx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def infile = assembly.name.endsWith('.gz') ? "<( gzip -dc $assembly )" : assembly
    """
    gx \\
        clean-genome \\
        --input $infile \\
        --action-report $action_report \\
        --output ${prefix}.clean.fasta \\
        --contam-fasta-out ${prefix}.contaminants.fasta
    """
}
