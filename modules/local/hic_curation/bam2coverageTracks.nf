process BAM2COVERAGE_TRACKS {
    tag "$meta.id"
    label 'process_medium'

    container = "community.wave.seqera.io/library/bedtools_pysam_samtools_ucsc-bedgraphtobigwig_pruned:cc5ee8bcdaac18a3"

    input:
    tuple val(meta), path(bam)
    path(chrom_sizes)

    output:
    tuple val(meta), path("*_coverage.bed")            , emit: bed
    tuple val(meta), path("*_coverage.bw")             , emit: bigwig
    tuple val(meta), path("*_coverage_capped.bed")     , emit: capped_bed
    tuple val(meta), path("*_coverage_capped.hitile")  , emit: capped_hitile    

    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def args5 = task.ext.args5 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # convert bam to sorted bed
    samtools sort -@ $task.cpus $args ${bam} | \\
    samtools view -@ $task.cpus $args2 - | \\
    genomeCoverageBed $args3 -ibam - | \\
    sort -k1,1 -k2,2n $args4  > ${prefix}_coverage.bed
    # convert bed to bigwig
    bedGraphToBigWig ${prefix}_coverage.bed ${chrom_sizes} ${prefix}_coverage.bw
    # coverage cap 
    sort -k1,1V -k2,2n -k3,3n $args4 ${prefix}_coverage.bed | \\
    awk -v COVERAGE_CAP=$params.hifi_coverage_cap 'BEGIN {FS = " "; OFS = "\t";}{if(int(\$4)>COVERAGE_CAP){print \$1,\$2,\$3,COVERAGE_CAP} else {print \$1,\$2,\$3,\$4}}' > ${prefix}_coverage_capped.bed
    # create coverage hitile 
    clodius aggregate bedgraph --chromsizes-filename ${chrom_sizes} ${args5} ${prefix}_coverage_capped.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //')        
        sort:     \$(echo \$(sort --version 2>&1) | head -n 1 | sed 's/^.*sort //; s/Copyright.*\$//')
        awk:      \$(echo \$(awk --version 2>&1) | head -n 1 | sed 's/Copyright.*\$//')
    END_VERSIONS
    """
}
