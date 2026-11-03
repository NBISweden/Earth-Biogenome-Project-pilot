nextflow.enable.types = true

process CREATE_GAP_TRACKS {
    tag "${id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/pysam_pip_awk_clodius:79d5112d656076f9"

    ext [
        args: { gap_args ?: '' },
        args2: { gap_args2 ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        bed: Path,
        chrom_sizes: Path,
        gap_args: String?,
        gap_args2: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        bedgraph: file("*_gaps.bedgraph"),
        beddb: file("*_gaps.bedgraph.beddb")
    )

    topic:
    tuple(task.process, 'sort', eval("sort --version | sed '1!d; s/.* //'")) >> 'versions'
    tuple(task.process, 'awk', eval("awk --version | sed '1!d; s/mawk //; s/ .*//'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${id}"

    """
    # create gaps bedgraph file which can be loaded into Pretext
    awk '{print \$0"\t"sqrt((\$3-\$2)*(\$3-\$2))}' ${bed} | \\
    sort -k1,1V -k2,2n -k3,3n --parallel=${task.cpus} ${args} > ${prefix}_gaps.bedgraph
    # create sorted.bedgraph.beddb which can be ingested to HiGlass
    clodius aggregate bedfile ${args2} --chromsizes-filename ${chrom_sizes} ${prefix}_gaps.bedgraph
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}_gaps.bedgraph ${prefix}_gaps.bedgraph.beddb
    """
}
