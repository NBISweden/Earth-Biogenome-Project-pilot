process KMC_DUMP {

    conda "${params.enable_conda ? 'bioconda::kmc==3.1.2rc1' : '' }"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/kmc:3.1.2rc1--h2d02072_0' :
              'quay.io/biocontainers/kmc:3.1.2rc1--h2d02072_0' }"

    input:
    tuple val(meta), path(count_db), val(interval)

    output:
    tuple val(meta), path("*.hist"), emit: histogram
    path "versions.yml"            , emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    kmc_dump -ci${interval[0]} -cx${interval[1]} \\
        $args \\
        $count_db \\
        $prefix.${interval[0]}-${interval[1]}.hist

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        kmc_dump: \$( kmc_dump --version | sed '1 !d;s/[^0-9]*\\(\\([0-9]\\.\\)\\{0,4\\}[0-9][^.]\\).*/\\1/' )
    END_VERSIONS
    """
}
