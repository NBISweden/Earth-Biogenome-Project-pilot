process KMC_DUMP {

    conda "${params.enable_conda ? 'bioconda::kmc==3.1.2rc1' : '' }"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/kmc:3.1.2rc1--h2d02072_0' :
              'quay.io/biocontainers/kmc:3.1.2rc1--h2d02072_0' }"

    input:
    tuple val(meta), path(count_db), val(lower_bound), val(upper_bound)

    output:
    tuple val(meta), path("*.hist"), emit: histogram
    path "versions.yml"            , emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    kmc_dump -ci${lower_bound} -cx${upper_bound} \\
        $args \\
        $prefix \\
        $prefix.${lower_bound}-${upper_bound}.hist

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        kmc_dump: \$( kmc_dump --version | sed '1 !d;s/[^0-9]*\\(\\([0-9]\\.\\)\\{0,4\\}[0-9][^.]\\).*/\\1/' )
    END_VERSIONS
    """
}
