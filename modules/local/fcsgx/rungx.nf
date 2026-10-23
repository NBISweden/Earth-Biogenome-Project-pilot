nextflow.enable.types = true

process FCSGX_RUNGX {
    tag "${id}"
    label 'process_high'

    ext [
        args: { fcsgx_args ?: '' },
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
        taxid: String,
        assembly: Path,
        gxdb: Path,
        ramdisk_path: String?,
        fcsgx_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        fcs_gx_report: file("*.fcs_gx_report.txt"),
        taxonomy_report: file("*.taxonomy.rpt"),
        hits: file("*.hits.tsv.gz", optional: true)
    )

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'fcs_gx', '0.5.4') >> 'versions'

    script:
    def prefix = task.ext.prefix
    def database = ramdisk_path ?: "${gxdb}" // maxForks set to 1 to limit concurrent execution corruption
    def copy_database = ramdisk_path ?
    """
    if test -d "${database}"; then
        echo "ERROR: Database exists in memory, and may be in use by another process" >&2
        ls -l ${database}
        exit 1
    fi
    # Clean up shared memory on exit
    trap "rm -rf ${database}" EXIT
    # Copy DB to RAM-disk when supplied. Otherwise, rungx is very slow.
    rclone copy $gxdb ${database}

    """ : ""
    """
    ${copy_database}
    export GX_NUM_CORES=$task.cpus
    run_gx.py \\
        --fasta $assembly \\
        --gx-db $database \\
        --tax-id $taxid \\
        --out-basename $prefix \\
        --out-dir . \\
        ${task.ext.args}
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}.fcs_gx_report.txt ${prefix}.taxonomy.rpt
    echo "" | gzip > ${prefix}.hits.tsv.gz
    """
}
