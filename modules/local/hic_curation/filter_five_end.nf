process FILTER_FIVE_END {
    tag "$meta.id"
    label 'process_medium'

    container = "community.wave.seqera.io/library/samtools_perl:e2ce3f7265547d1f"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam")  , emit: bam
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    samtools view -h -@ $task.cpus ${bam} | \\
    perl ${projectDir}/bin/filter_five_end.pl | \\
    samtools view -@ $task.cpus -Sb - > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
