process FCSGX_CLEAN {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.0--h4ac6f70_3':
        'biocontainers/ncbi-fcs-gx:0.5.0--h4ac6f70_3' }"

    input:
    tuple val(meta), path(assembly), path(action_report)

    output:
    tuple val(meta), path("*.clean.fasta")       , emit: clean_fasta
    tuple val(meta), path("*.contaminants.fasta"), emit: contaminants_fasta
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.5.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    gx \\
        clean-genome \\
        --input $assembly \\
        --action-report $action_report \\
        --output ${prefix}.clean.fasta \\
        --contam-fasta-out ${prefix}.contaminants.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs_gx: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.5.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.clean.fasta
    touch ${prefix}.contaminants.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs_gx: $VERSION
    END_VERSIONS
    """
}
