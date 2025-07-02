process BAM2COVERAGE_TRACKS {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/bedtools_samtools_ucsc-bedgraphtobigwig:ba8edacb0d7bff96"

    input:
    tuple val(meta), path(bam)
    path(chrom_sizes)
    val(hifi_coverage_cap)

    output:
    tuple val(meta), path("*_coverage.bed")            , emit: bed
    tuple val(meta), path("*_coverage.bw")             , emit: bigwig
    tuple val(meta), path("*_coverage_capped.bed")     , emit: capped_bed
    tuple val(meta), path("*_coverage_capped.bw")      , emit: capped_bigwig

    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # convert bam to sorted bed
    samtools sort -@ $task.cpus $args ${bam} | \\
    samtools view -@ $task.cpus $args2 - | \\
    genomeCoverageBed $args3 -ibam - | \\
    sort -k1,1 -k2,2n --parallel=${task.cpus} $args4  > ${prefix}_coverage.bed
    # convert bed to bigwig
    bedGraphToBigWig ${prefix}_coverage.bed ${chrom_sizes} ${prefix}_coverage.bw
    # coverage cap
    sort -k1,1V -k2,2n -k3,3n $args4 ${prefix}_coverage.bed | \\
    awk -v COVERAGE_CAP=$hifi_coverage_cap 'BEGIN {FS = " "; OFS = "\t";}{if(int(\$4)>COVERAGE_CAP){print \$1,\$2,\$3,COVERAGE_CAP} else {print \$1,\$2,\$3,\$4}}' > ${prefix}_coverage_capped.bed
    # convert bed to bigwig
    bedGraphToBigWig ${prefix}_coverage_capped.bed ${chrom_sizes} ${prefix}_coverage_capped.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/.* //')
        bedtools: \$(bedtools --version | sed 's/.*v//')
        sort:     \$(sort --version | sed '1!d; s/.* //')
        awk:      \$(awk --version | sed '1!d; s/mawk //; s/ .*//')
    END_VERSIONS
    """
}
