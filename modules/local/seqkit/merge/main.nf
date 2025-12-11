process SEQKIT_MERGE {
    tag "${meta.id}"
    label 'process_low'
    // File IO can be a bottleneck. See: https://bioinf.shenwei.me/seqkit/usage/#parallelization-of-cpu-intensive-jobs

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(fastx)
    val out_ext

    output:
    tuple val(meta), path("${prefix}.*"), emit: fastx
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def extension = out_ext ?: "fasta"
    def call_gzip = extension.endsWith('.gz') ? "| gzip -c ${args2}" : ''
    """
    for FA in ${fastx}; do
        LABEL=\$(echo "\$FA" | grep -o 'hap[0-9]\\+')
        seqkit \\
            replace \\
            -p "^" -r "\${LABEL}_" \\
            --threads ${task.cpus} \\
            ${args} \\
            \$FA
    done ${call_gzip} > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}
