process SMUDGEPLOT_HETKMERS {

    conda "${task.ext.enable_conda ? 'bioconda::smudgeplot==0.2.4' : '' }"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/smudgeplot:0.2.4--py39r40h779adbc_1' :
              'quay.io/biocontainers/smudgeplot:0.2.4--py39r40h779adbc_1' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("_coverages.tsv"), emit: coverage_tsv
    path "versions.yml"                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    // def args   = task.ext.args   ?: ''
    """
    smudgeplot.py hetkmers \\
        -o $prefix < $histogram

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        tool: \$(smudgeplot.py -v |& sed '1 !d;s/[^0-9]*\(\([0-9]\.\)\{0,4\}[0-9]\).*/\1/' )
    END_VERSIONS
    """
}
