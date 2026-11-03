nextflow.enable.types = true

process CREATE_CHROMOSOME_SIZES_FILE {
    tag "${id}"
    label 'process_low'

    container "community.wave.seqera.io/library/pip_awk:e0daab0638d06dfd"

    ext [
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        fai: Path,
        hic_map_sort_by: String,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        sizes: file("*.sizes")
    )

    topic:
    tuple(task.process, 'sort', eval("sort --version | sed '1!d; s/.* //'")) >> 'versions'
    tuple(task.process, 'awk', eval("awk --version | sed '1!d; s/mawk //; s/ .*//'")) >> 'versions'

    script:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    if [[ "$hic_map_sort_by" == "length" ]]
    then
        awk '{print \$1\"\t\"\$2}' ${fai} | \\
        sort -k2,2 -nr > ${prefix}.sizes
    elif [[ "$hic_map_sort_by" == "name" ]]
    then
        awk '{print \$1\"\t\"\$2}' ${fai} | \\
        sort -k1,1 -nr > ${prefix}.sizes
    else
        awk '{print \$1\"\t\"\$2}' ${fai} > ${prefix}.sizes
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}.sizes
    """
}
