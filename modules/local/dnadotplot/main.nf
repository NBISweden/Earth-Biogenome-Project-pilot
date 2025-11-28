process DNADOTPLOT {
    tag "$meta.id"
    label 'process_low'
    container "oras://community.wave.seqera.io/library/dnadotplot:0.1.4--aa4389449e8456fc"

    input:
    tuple val(meta), path(ref), path(query), val(experiment)

    output:
    tuple val(meta), path("*.svg"), emit: svg
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dnadotplot \\
        $args \\
        -1 $ref \\
        -2 $query \\
        -o seq1-${ref.baseName}_seq2-${query.baseName}.${experiment}.svg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnadotplot: \$(dnadotplot --version | sed "s/dnadotplot //")
    END_VERSIONS
    """
}
