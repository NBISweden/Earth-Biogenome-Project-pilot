process ENTREZDIRECT_ESEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    tuple val(meta), val(species)
    val database

    output:
    tuple val(meta), val(species), env(TAXONOMY_COUNT), emit: count
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    # Get count of taxonomy hits for the query species name
    TAXONOMY_COUNT=\$(esearch \\
        -db $database \\
        -query "$species" \\
        $args | \\
    xtract \\
        -pattern ENTREZ_DIRECT \\
        -element Count \\
        $args2)

    export TAXONOMY_COUNT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esearch: \$(esearch -version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.xml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esearch: \$(esearch -version)
    END_VERSIONS
    """
}
