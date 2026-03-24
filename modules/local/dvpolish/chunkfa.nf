process DVPOLISH_CHUNKFA {
    tag "$meta.id"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"


    input:
    tuple val(meta), path(fai)
    val chunk_size

    output:
    tuple val(meta), path ("*.bed", arity: '1..*')        , emit: bed
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_size_bytes = MemoryUnit.of(chunk_size).toBytes()

    """
    # convert chunk size into base pairs
    awk -v chunk_size_bases=${chunk_size_bytes} -v prefix=${prefix} '
    BEGIN {
        block=1
        cumulative_bases=0
    }
    {
        output_file = sprintf("%s_chunk_%d.bed", prefix, block)
        printf("%s\\t0\\t%s\\n", \$1, \$2) > output_file
        cumulative_bases+=\$2

        if (cumulative_bases >= chunk_size_bases) {
            cumulative_bases=0
            block++
        }
    }
    ' ${fai}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(awk -W version |& sed '1!d; s/mawk //; s/ .*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_chunk_1.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(awk -W version |& sed '1!d; s/mawk //; s/ .*//')
    END_VERSIONS
    """
}
