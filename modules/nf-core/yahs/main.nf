nextflow.enable.types = true

process YAHS {
    tag "${id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yahs:1.2.2--h577a1d6_1':
        'biocontainers/yahs:1.2.2--h577a1d6_1' }"

    ext [
        args: { yahs_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        fasta: Path,
        fai: Path,
        hic_map: Path,
        agp: Path?,
        yahs_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        scaffolds_fasta: file("${task.ext.prefix}_scaffolds_final.fa", optional: true),
        scaffolds_agp: file("${task.ext.prefix}_scaffolds_final.agp", optional: true),
        initial_break_agp: files("${task.ext.prefix}_{inital,no}_break*.agp", optional: true),
        round_agp: files("${task.ext.prefix}_r*_*.agp", optional: true),
        binary: file("${task.ext.prefix}.bin"),
        log: file("${task.ext.prefix}.log")
    )

    topic:
    tuple(task.process, 'yahs', eval("yahs --version 2>&1")) >> 'versions'

    script:
    def args      = task.ext.args
    def prefix    = task.ext.prefix
    def agp_input = agp ? "-a ${agp}" : ""
    """
    yahs \\
        -o ${prefix} \\
        ${agp_input} \\
        ${args} \\
        ${fasta} \\
        ${hic_map} \\
        2>| >( tee ${prefix}.log >&2 )
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}_scaffolds_final.fa
    touch ${prefix}_scaffolds_final.agp
    touch ${prefix}_inital_break_01.agp
    touch ${prefix}_no_break.agp
    touch ${prefix}_r01.agp
    touch ${prefix}_r01_break.agp
    touch ${prefix}.bin
    touch ${prefix}.log
    """
}
