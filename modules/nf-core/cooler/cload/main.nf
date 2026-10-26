nextflow.enable.types = true

process COOLER_CLOAD {
    tag "${id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.10.4--pyhdfd78af_0' :
        'biocontainers/cooler:0.10.4--pyhdfd78af_0' }"

    ext [
        args: { cooler_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        contacts: Path,
        index: Path?,
        chromsizes: Path,
        mode: String,
        cool_bin: Integer,
        cooler_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        cool: file("*.cool")
    )

    topic:
    tuple(task.process, 'cooler', eval('cooler --version 2>&1 | sed "s/cooler, version //"')) >> 'versions'

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    def nproc  = mode in ["pairix", "tabix"] ? "--nproc ${task.cpus}" : ""
    """
    cooler cload ${mode} \\
        ${args} \\
        ${nproc} \\
        ${chromsizes}:${cool_bin} \\
        ${contacts} \\
        ${prefix}.cool
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}.cool
    """
}
