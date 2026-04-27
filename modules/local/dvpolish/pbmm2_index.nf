process DVPOLISH_PBMM2_INDEX {
    label 'process_medium'

    // Note: the versions here need to match the versions used in pbmm2/align
    conda 'bioconda::pbmm2=1.13.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmm2:1.13.1--h9ee0642_0' :
        'biocontainers/pbmm2:1.13.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi"), emit: index
    tuple val("${task.process}"), val('pbmm2'), eval("pbmm2 --version 2>&1 | sed '1!d; s/pbmm2 //'"), emit: versions_pbmm2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    pbmm2 index\\
        -j $task.cpus \\
        $args \\
        $fasta \\
        ${fasta.baseName}.mmi
    """
}
