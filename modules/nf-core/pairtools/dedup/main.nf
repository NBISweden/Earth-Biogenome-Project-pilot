process PAIRTOOLS_DEDUP {
    tag "$meta.id"
    label 'process_high'

    // Pinning numpy to 1.23 until https://github.com/open2c/pairtools/issues/170 is resolved
    // Not an issue with the biocontainers because they were built prior to numpy 1.24
    conda "bioconda::pairtools=1.1.0 conda-forge::numpy=1.23"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.1.0--py39hd5a99d8_2' :
        'biocontainers/pairtools:1.1.0--py39hd5a99d8_2' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.pairs.gz")  , emit: pairs
    tuple val(meta), path("*.pairs.stat"), emit: stat
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pairtools dedup \\
        $args \\
        -o ${prefix}.pairs.gz \\
        --output-stats ${prefix}.pairs.stat \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version | tr '\\n' ',' | sed 's/.*pairtools.*version //' | sed 's/,\$/\\n/')
    END_VERSIONS
    """
}
