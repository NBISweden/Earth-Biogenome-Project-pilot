process FCSGX_RUNGX {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.0--h4ac6f70_3':
        'biocontainers/ncbi-fcs-gx:0.5.0--h4ac6f70_3' }"

    input:
    tuple val(meta), val(taxid), path(assembly)
    path gxdb

    output:
    tuple val(meta), path("*.fcs_gx_report.txt"), emit: fcs_gx_report
    tuple val(meta), path("*.taxonomy.rpt")     , emit: taxonomy_report
    tuple val(meta), path("*.hits.tsv.gz")      , emit: hits, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.5.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    export GX_NUM_CORES=$task.cpus
    run_gx.py \\
        --fasta $assembly \\
        --gx-db $gxdb \\
        --tax-id $taxid \\
        --out-basename $prefix \\
        --out-dir . \\
        $args

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
    touch ${prefix}.fcs_gx_report.txt
    touch ${prefix}.taxonomy.rpt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs_gx: $VERSION
    END_VERSIONS
    """
}
