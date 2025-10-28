process ENTREZDIRECT_ESEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    tuple val(meta), val(term)
    val database

    output:
    val "success"      , emit: success
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    # Get count of taxonomy hits for species name
    TAXONOMY_HITS=\$(esearch \\
        -db $database \\
        -query "$term" \\
        $args | \\
    xtract \\
        -pattern ENTREZ_DIRECT \\
        -element Count \\
        $args2)

    # Test that only one taxonomy record was found
    # TAXONOMY_HITS = 0: no records found
    # TAXONOMY_HITS > 1: multiple records found
    if [ "\$TAXONOMY_HITS" == "0" ]; then
        echo "Error: No taxonomy records found for query: $term" >&2
        exit 1
    else
        if [ "\$TAXONOMY_HITS" -gt "1" ]; then
            echo "Error: More than one appropriate taxonomy entry for query: $term" >&2
            exit 1
        fi
    fi

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
