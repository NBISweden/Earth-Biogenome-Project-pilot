process GENOMESCOPE {

    conda "${params.enable_conda ? 'bioconda::genomescope2==2.0' : '' }"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/genomescope2:2.0--py39r40hdfd78af_4' :
              'quay.io/biocontainers/genomescope2:2.0--py39r40hdfd78af_4' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("*linear_plot.png")            , emit: linear_plot_png
    tuple val(meta), path("*transformed_linear_plot.png"), emit: transformed_linear_plot_png
    tuple val(meta), path("*log_plot.png")               , emit: log_plot_png
    tuple val(meta), path("*transformed_log_plot.png")   , emit: transformed_log_plot_png
    tuple val(meta), path("*model.txt")                  , emit: model
    tuple val(meta), path("*summary.txt")                , emit: summary
    path "versions.yml"                                  , emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    genomescope2 \\
        -i $histogram \\
        $args \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        tool: \$( genomescope2 -v | sed '1 !d;s/[^0-9]*\\(\\([0-9]\\.\\)\\{0,4\\}[0-9]\\).*/\\1/' )
    END_VERSIONS
    """
}
