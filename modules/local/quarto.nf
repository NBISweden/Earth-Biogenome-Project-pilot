process QUARTO {
    label 'process_single'

    conda "conda-forge::quarto=1.3.450 conda-forge::r-base=4.3.2 conda-forge::r-tidyverse=2.0.0"
    container "docker.io/rocker/verse:4.3.2"

    input:
    path notebook
    path data, stageAs: 'data/*', arity: '1..*'

    output:
    path "*.html"       , emit: html
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // TODO add prefix to script
    def prefix = task.ext.args ?: notebook.baseName
    """
    quarto \\
        render \\
        $args \\
        $notebook

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$( quarto --version )
    END_VERSIONS
    """
}
