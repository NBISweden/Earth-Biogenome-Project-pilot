process QUARTO_NOTEBOOK {
    tag "$notebook.baseName"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // container "docker.io/rocker/verse:4.3.2"

    input:
    tuple val(meta), path(notebook, arity: '1')
    path (log_files, stageAs: 'log_files/*', arity: '1..*')
    path 'params.yml'

    output:
    path "*.html"      , arity: '1', emit: html
    path "versions.yml", arity: '1', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = notebook.baseName
    """
    export USERID=\$UID
    export XDG_CACHE_HOME=tmp/quarto_cache_home
    export XDG_DATA_HOME=tmp/quarto_data_home

    # Link params to meta data
    ln -s params.yml _quarto.yml

    quarto check
    quarto \\
        render \\
        $notebook \\
        --execute-params params.yml \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        multiqc: \$(multiqc --version | sed '1!d; s/.*version //')
    END_VERSIONS
    """
}
