process WGET {
    tag "$url"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/fabc54c698bff97dd3cc4fdd16c6f5d9235f37675c4c08bc3b55ed98a5037e4f/data':
        'community.wave.seqera.io/library/wget:1.21.4--ec3657d9d2ab4987' }"

    input:
    val url

    output:
    path("$prefix"),     emit: download
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: url.tokenize('/').last()
    """
    wget \\
        ${args} \\
        -O ${prefix} \\
        ${url}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | sed '1!d; s/.*Wget //; s/ .*//')
    END_VERSIONS
    """

    stub:
    // def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: url.tokenize('/').last()
    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | sed '1!d; s/.*Wget //; s/ .*//')
    END_VERSIONS
    """
}
