process CREATE_GAP_TRACKS {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/pysam_pip_awk_clodius:79d5112d656076f9"

    input:
    tuple val(meta), path(bed)
    path(chrom_sizes)

    output:
    tuple val(meta), path("*_gaps.bedgraph")           , emit: bedgraph
    tuple val(meta), path("*_gaps.bedgraph.beddb")     , emit: beddb

    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # create gaps bedgraph file which can be loaded into Pretext
    awk '{print \$0"\t"sqrt((\$3-\$2)*(\$3-\$2))}' ${bed} | \\
    sort -k1,1V -k2,2n -k3,3n --parallel=${task.cpus} ${args} > ${prefix}_gaps.bedgraph
    # create sorted.bedgraph.beddb which can be ingested to HiGlass
    clodius aggregate bedfile ${args2} --chromsizes-filename ${chrom_sizes} ${prefix}_gaps.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort:     \$(sort --version | sed '1!d; s/.* //')
        awk:      \$(awk --version | sed '1!d; s/mawk //; s/ .*//')
    END_VERSIONS
    """
}
