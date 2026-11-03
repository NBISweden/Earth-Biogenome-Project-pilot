nextflow.enable.types = true

process SEQKIT_MERGE {
    tag "${id}"
    label 'process_low'
    // File IO can be a bottleneck. See: https://bioinf.shenwei.me/seqkit/usage/#parallelization-of-cpu-intensive-jobs

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    ext [
        args: { seqkit_args ?: '' },
        args2: { seqkit_args2 ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        fastx: List<Path>,
        out_ext: String?,
        seqkit_args: String?,
        seqkit_args2: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        fastx: file("${task.ext.prefix}.*")
    )

    topic:
    tuple(task.process, 'seqkit', eval('seqkit version | cut -d" " -f2')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix
    def extension = out_ext ?: "fasta"
    def call_gzip = extension.endsWith('.gz') ? "| gzip -c ${args2}" : ''
    """
    for FA in ${fastx.join(' ')}; do
        LABEL=\$(echo "\$FA" | grep -o 'hap[0-9]\\+')
        seqkit \\
            replace \\
            -p "^" -r "\${LABEL}_" \\
            --threads ${task.cpus} \\
            \$FA
    done | seqkit seq \\
        --threads ${task.cpus} \\
        ${args} \\
    ${call_gzip} > ${prefix}.${extension}
    """

    stub:
    def prefix = task.ext.prefix
    def extension = out_ext ?: "fasta"
    """
    touch ${prefix}.${extension}
    """
}
