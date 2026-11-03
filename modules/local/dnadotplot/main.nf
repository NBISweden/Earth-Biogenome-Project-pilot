nextflow.enable.types = true

process DNADOTPLOT {
    tag "${id}"
    label 'process_low'

    ext [
        args: { dnadotplot_args ?: '' }
    ]

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/dnadotplot:0.1.4--aa4389449e8456fc' :
        'community.wave.seqera.io/library/dnadotplot:0.1.4--165aba59778a399e' }"

    input:
    record(
        id: String,
        meta: Record?,
        ref: Path,
        query: Path,
        experiment: String,
        dnadotplot_args: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        svg: file("*.svg")
    )

    topic:
    tuple(task.process, 'dnadotplot', eval("dnadotplot --version | sed 's/dnadotplot //'")) >> 'versions'

    script:
    """
    dnadotplot \\
        ${task.ext.args} \\
        -1 $ref \\
        -2 $query \\
        -o seq1-${ref.baseName}_seq2-${query.baseName}.${experiment}.svg
    """

    stub:
    """
    touch seq1-${ref.baseName}_seq2-${query.baseName}.${experiment}.svg
    """
}
