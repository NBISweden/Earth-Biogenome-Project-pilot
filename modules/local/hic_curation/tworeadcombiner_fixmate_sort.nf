process TWOREADCOMBINER_FIXMATE_SORT {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/samtools_perl:e2ce3f7265547d1f"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam")  , emit: bam
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | sed "1!d; s/samtools //"'), emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('perl'), eval("perl --version | sed '2!d; s/.*(v//; s/).*//'"), emit: versions_perl, topic: versions
    tuple val("${task.process}"), val('grep'), eval("grep --version | sed '1!d; s/.* //'"), emit: versions_grep, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    two_read_bam_combiner_sanger.pl ${bam} samtools ${args} | \\
    grep -v -e "^@HD" | \\
    samtools sort ${args2} -@${task.cpus} -T sort_tmp -o ${prefix}_comb.bam -
    """
}
