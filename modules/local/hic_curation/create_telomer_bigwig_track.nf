nextflow.enable.types = true

process CREATE_TELOMER_BIGWIG_TRACK {
    tag "${id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/ucsc-bedgraphtobigwig:469--f66c00635e28f216"

    ext [
        args: { telomer_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        bedgraph: Path,
        chrom_sizes: Path,
        telomer_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        bw: file("*_telomer.bw")
    )

    topic:
    tuple(task.process, 'bedGraphToBigWig', eval("bedGraphToBigWig |& sed '1!d; s/.* v //; s/ .*//'")) >> 'versions'

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    sort -k1,1V -k2,2n -k3,3n ${args} ${bedgraph} > ${bedgraph}_sorted.bedgraph
    bedGraphToBigWig ${bedgraph}_sorted.bedgraph ${chrom_sizes} ${prefix}_telomer.bw;
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}_telomer.bw
    """
}
