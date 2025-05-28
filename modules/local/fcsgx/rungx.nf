process FCSGX_RUNGX {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::ncbi-fcs-gx=0.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.4--h4ac6f70_1':
        'biocontainers/ncbi-fcs-gx:0.5.4--h4ac6f70_1' }"

    input:
    tuple val(meta), val(taxid), path(assembly)
    path gxdb
    val ramdisk_path

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
    def database = ramdisk_path ?: gxdb // maxForks set to 1 to limit concurrent execution corruption
    def VERSION = '0.5.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    (ramdisk_path ?
    """
    if test -d "${database}"; then
        echo "ERROR: Database exists in memory, and may be in use by another process" >&2
        ls -l ${database}
        exit 1
    fi
    # Clean up shared memory on exit
    trap "rm -rf "${database}" EXIT
    # Copy DB to RAM-disk when supplied. Otherwise, rungx is very slow.
    rclone copy $gxdb ${database}

    """ : "")
    <<
    """
    export GX_NUM_CORES=$task.cpus
    run_gx.py \\
        --fasta $assembly \\
        --gx-db $database \\
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
    def VERSION = '0.5.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.fcs_gx_report.txt
    touch ${prefix}.taxonomy.rpt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs_gx: $VERSION
    END_VERSIONS
    """
}
