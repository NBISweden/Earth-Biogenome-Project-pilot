process SMUDGEPLOT_CUTOFF {

    conda "${task.ext.enable_conda ? 'bioconda::smudgeplot==0.2.4' : '' }"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/smudgeplot:0.2.4--py39r40h779adbc_1' :
              'quay.io/biocontainers/smudgeplot:0.2.4--py39r40h779adbc_1' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), env([LOWER_BOUND, UPPER_BOUND]), emit: bounds
    path "versions.yml"                             , emit: versions

    script:
    // def prefix = task.ext.prefix ?: meta.id
    // def args   = task.ext.args   ?: ''
    """
    LOWER_BOUND=\$( smudgeplot.py cutoff $histogram L )
    UPPER_BOUND=\$( smudgeplot.py cutoff $histogram U )

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        tool: \$( smudgeplot.py -v |& sed '1 !d;s/[^0-9]*\(\([0-9]\.\)\{0,4\}[0-9]\).*/\1/' )
    END_VERSIONS
    """
}
