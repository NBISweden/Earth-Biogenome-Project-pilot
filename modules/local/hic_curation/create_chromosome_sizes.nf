process CREATE_CHROMOSOME_SIZES_FILE {
    tag "$meta.id"
    label 'process_low'

    container "community.wave.seqera.io/library/pip_awk:e0daab0638d06dfd"

    input:
    tuple val(meta), path(fai)

    output:
    tuple val(meta), path("*.sizes")  , emit: sizes
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    if [[ "$params.hic_map_sort_by" == "length" ]]
    then 
        awk '{print \$1\"\t\"\$2}' ${fai} | \\
        sort -k2,2 -nr > ${prefix}.sizes
    elif [[ "$params.hic_map_sort_by" == "name" ]]
    then 
        awk '{print \$1\"\t\"\$2}' ${fai} | \\
        sort -k1,1 -nr > ${prefix}.sizes
    else 
        awk '{print \$1\"\t\"\$2}' ${fai} > ${prefix}.sizes
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort:     \$(echo \$(sort --version 2>&1) | head -n 1 | sed 's/^.*sort //; s/Copyright.*\$//')
        awk:      \$(echo \$(awk --version 2>&1) | head -n 1 | sed 's/Copyright.*\$//')
    END_VERSIONS
    """
}
