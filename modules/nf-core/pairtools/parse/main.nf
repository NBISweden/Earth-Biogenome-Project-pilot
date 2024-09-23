process PAIRTOOLS_PARSE {
    tag "$meta.id"
    label 'process_high'

    // Pinning numpy to 1.23 until https://github.com/open2c/pairtools/issues/170 is resolved
    // Not an issue with the biocontainers because they were built prior to numpy 1.24
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.1.0--py39hd5a99d8_2' :
        'biocontainers/pairtools:1.1.0--py39hd5a99d8_2' }"

    input:
    tuple val(meta), path(bam)
    path chromsizes

    output:
    tuple val(meta), path("*.pairsam.gz")  , emit: pairsam
    tuple val(meta), path("*.pairsam.stat"), emit: stat
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    int decompression_cpus = Math.max(1, Math.ceil(task.cpus/4.0) as int)
    int compression_cpus = Math.max(1, task.cpus-decompression_cpus-1) as int
    """
    pairtools \\
        parse \\
        -c $chromsizes \\
        --nproc-in $decompression_cpus \\
        --nproc-out $compression_cpus \\
        $args \\
        --output-stats ${prefix}.pairsam.stat \\
        -o ${prefix}.pairsam.gz \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version | tr '\\n' ',' | sed 's/.*pairtools.*version //' | sed 's/,\$/\\n/')
    END_VERSIONS
    """
}
