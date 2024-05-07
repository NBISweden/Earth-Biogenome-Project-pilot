process QUARTO {
    tag "$notebook.baseName"
    label 'process_single'

    conda "conda-forge::quarto=1.3.450 conda-forge::r-base=4.3.2 conda-forge::r-tidyverse=2.0.0"
    container "docker.io/rocker/verse:4.3.2"

    input:
    tuple (
        val (meta),
        path (notebook),
        path (data, stageAs: 'data/*', arity: '1..*')
    )

    output:
    path "*.html"      , emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.args ?: notebook.baseName
    """
    USERID=\$UID
    XDG_CACHE_HOME=tmp/quarto_cache_home
    XDG_DATA_HOME=tmp/quarto_data_home

    quarto \\
        render \\
        $notebook \\
        $args \\
        --output ${prefix}.html

    rm -rf tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$( quarto --version )
    END_VERSIONS
    """
}
