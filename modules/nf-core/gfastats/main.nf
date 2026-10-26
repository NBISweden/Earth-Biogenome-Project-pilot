nextflow.enable.types = true

process GFASTATS {
    tag "${id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7150110bea918bd523e41aab58f3be2eda949bbe7f715500681b01965e5f667/data':
        'community.wave.seqera.io/library/gfastats:1.3.11--692f8595aa6b06a8' }"

    ext [
        args: { gfastats_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        assembly: Path,
        out_fmt: String?,
        genome_size: String?,
        target: String,
        agpfile: Path?,
        include_bed: Path?,
        exclude_bed: Path?,
        instructions: Path?,
        gfastats_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        assembly_summary: file("*.assembly_summary"),
        assembly: file("*.${out_fmt}.gz", optional: true)
    )

    topic:
    tuple(task.process, 'gfastats', eval("gfastats -v | sed '1!d;s/.*v//'")) >> 'versions'

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    def agp    = agpfile ? "--agp-to-path $agpfile" : ""
    def ibed   = include_bed ? "--include-bed $include_bed" : ""
    def ebed   = exclude_bed ? "--exclude-bed $exclude_bed" : ""
    def sak    = instructions ? "--swiss-army-knife $instructions" : ""
    def output_sequences = out_fmt ? "--out-format ${prefix}.${out_fmt}.gz" : ""
    """
    gfastats \\
        $assembly \\
        ${genome_size ?: ''} \\
        $target \\
        $args \\
        --threads $task.cpus \\
        $agp \\
        $ibed \\
        $ebed \\
        $sak \\
        $output_sequences \\
        > ${prefix}.assembly_summary
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    echo | gzip > ${prefix}.${out_fmt ?: "fa"}.gz
    touch ${prefix}.assembly_summary
    """
}
