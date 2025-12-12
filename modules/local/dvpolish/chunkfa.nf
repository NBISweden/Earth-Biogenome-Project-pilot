process DVPOLISH_CHUNKFA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"


    input:
    tuple val(meta), path(fai)

    output:
    tuple val(meta), path ("*.bed", arity: '1..*')        , emit: bed
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Split the string by spaces, then iterate to create a map
    // assumption all arguments are provided like key1 = value1 key2 = value2 ...
    def args_map = [:]
    args.replaceAll(' = ', '=').split().each {
        def (key, value) = it.split('=')
        args_map[key.trim()] = value.trim()
    }
    def chunk_size = MemoryUnit.of(args_map['chunk_size'] ?: '90MB').toBytes()
    
    """
    # convert chunk size into base pairs
    awk -v chunk_size_inBases=${chunk_size} -v prefix=$prefix 'BEGIN {
        block=1
        cum_basecount=0
    }{
        output_file = sprintf("%s_chunk_%d.bed", prefix, block)
        printf("%s\\t0\\t%s\\n", \$1, \$2) > output_file
        cum_basecount+=\$2

        if (cum_basecount >= chunk_size_inBases)
        {
            cum_basecount=0
            block+=1
        }
    }' $fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(awk -W version |& head -n 1 | sed 's/mawk //; s/ .*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_chunk_1.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(awk -W version |& head -n 1 | sed 's/mawk //; s/ .*//')
    END_VERSIONS
    """
}
