process KAT_COMP {

    conda "${params.enable_conda ? 'bioconda::kat==2.4.2' : '' }"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/kat:2.4.2--py38hfc5f9d8_2' :
              'quay.io/biocontainers/kat:2.4.2--py38hfc5f9d8_2' }"

    input:
    // Mandatory
    tuple val(meta), path(reads), path(assembly) // [ [id: 'name'], [file(read1), file(read2)], file(assembly)]

    output:
    tuple val(meta), path('*-main.mx')               , emit: mx // Matrix of k-mer counts
    tuple val(meta), path("*-main.mx.spectra-cn.png"), emit: png // Copy number spectra plot
    tuple val(meta), path("*.dist_analysis.json")    , emit: json
    tuple val(meta), path("*.stats")                 , emit: stats
    path "versions.yml"                              , emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    kat comp \\
        -t $task.cpus \\
        $args \\
        -o ${prefix}.katcomp \\
        '$reads' \\
        $assembly

    cat <<-END_VERSIONS > versions.yml
    '${task.process}':
        kat: \$( kat --version | sed 's/kat //' )
    END_VERSIONS
    """
}
