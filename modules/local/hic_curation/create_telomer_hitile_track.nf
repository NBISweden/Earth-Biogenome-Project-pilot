process CREATE_TELOMER_HITILE_TRACK {
    tag "$meta.id"
    label 'process_medium'

    container = "community.wave.seqera.io/library/pysam_pip_awk_clodius:79d5112d656076f9"

    input:
    tuple val(meta), path(bedgraph)
    path(chrom_sizes)

    output:
    tuple val(meta), path("*_telomer.hitile")         , emit: hitile 
    
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    sort -k1,1V -k2,2n -k3,3n ${args} ${bedgraph} | \\
    clodius aggregate bedgraph ${args2} --chromsizes-filename ${chrom_sizes} -o ${prefix}_telomer.hitile -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort:     \$(echo \$(sort --version 2>&1) | sed 's/^.*sort //; s/Copyright.*\$//')
    END_VERSIONS
    """
}
