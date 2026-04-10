process DNADOTPLOT {
    tag "$meta.id"
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/dnadotplot:0.1.4--aa4389449e8456fc' :
        'community.wave.seqera.io/library/dnadotplot:0.1.4--165aba59778a399e' }"

    input:
    tuple val(meta), path(ref), path(query), val(experiment)

    output:
    tuple val(meta), path("*.svg"), emit: svg
    tuple val("${task.process}"), val('dnadotplot'), eval("dnadotplot --version | sed 's/dnadotplot //'"), topic: versions, emit: versions_dnadotplot

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
    """
}
