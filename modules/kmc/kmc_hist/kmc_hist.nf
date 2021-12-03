process KMC_HIST {

    conda "${params.enable_conda ? 'bioconda::kmc==3.1.2rc1' : '' }"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/kmc:3.1.2rc1--h2d02072_0' :
              'quay.io/biocontainers/kmc:3.1.2rc1--h2d02072_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("$prefix"), emit: count_db
    tuple val(meta), path("*.hist") , emit: histogram
    path "versions.yml"             , emit: versions

    script:
    prefix   = task.ext.prefix ?: meta.id
    def args = task.ext.args   ?: ''
    """
    mkdir kmc_workdir
    printf '%s\n' $reads > kmc_read_files.txt
    kmc $args \\
        -t$task.cpus \\
        -j$prefix.summary.json \\
        @kmc_read_files.txt \\
        $prefix \\
        kmc_workdir
    kmc_tools transform \\
        $prefix \\
        histogram ${prefix}.hist

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        kmc: \$( kmc --version | sed '1 !d;s/[^0-9]*\\(\\([0-9]\\.\\)\\{0,4\\}[0-9][^.]\\).*/\\1/' )
    END_VERSIONS
    """
}
