nextflow.enable.types = true

process DVPOLISH_CHUNKFA {
    tag "${id}"
    label 'process_single'

    ext [
        prefix: { output_prefix ?: "${id}" }
    ]

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    record(
        id: String,
        meta: Record?,
        fai: Path,
        chunk_size_bytes: Long,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        bed: files("*.bed")
    )

    topic:
    tuple(task.process, 'mawk', eval("awk -W version |& sed '1!d; s/mawk //; s/ .*//'")) >> 'versions'

    script:
    def prefix = task.ext.prefix
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
    """

    stub:
    """
    touch ${task.ext.prefix}_chunk_1.bed
    """
}
