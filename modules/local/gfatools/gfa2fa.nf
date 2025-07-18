process GFATOOLS_GFA2FA {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::gfatools=0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfatools:0.5--he4a0461_4':
        'biocontainers/gfatools:0.5--he4a0461_4' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gfatools \\
        gfa2fa \\
        $args \\
        $gfa \\
        | gzip -c > ${prefix}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfatools: \$( gfatools version | sed '1!d; s/.* //' )
    END_VERSIONS
    """
}
