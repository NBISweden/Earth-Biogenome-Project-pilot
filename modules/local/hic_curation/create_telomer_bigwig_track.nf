process CREATE_TELOMER_BIGWIG_TRACK {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/ucsc-bedgraphtobigwig:469--f66c00635e28f216"

    input:
    tuple val(meta), path(bedgraph)
    path(chrom_sizes)

    output:
    tuple val(meta), path("*_telomer.bw"), emit: bw
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sort -k1,1V -k2,2n -k3,3n ${args} ${bedgraph} > ${bedgraph}_sorted.bedgraph
    bedGraphToBigWig ${bedgraph}_sorted.bedgraph ${chrom_sizes} ${prefix}_telomer.bw;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedGraphToBigWig: \$(bedGraphToBigWig |& sed '1!d; s/.* v //; s/ .*//')
    END_VERSIONS
    """
}
