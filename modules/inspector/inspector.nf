process INSPECTOR {

    // TODO:: Update conda path to correct conda depedencies
    // conda "${task.ext.enable_conda ? 'bioconda::tool=0.0.0' : '' }"
    // TODO:: Update Singularity and Docker paths to correct container paths
    container 'ghcr.io/nbisweden/earth-biogenome-project-pilot/inspector:1.0'
    // container "${workflow.containerEngine == 'singularity' &&
    //               !task.ext.singularity_pull_docker_container ?
    //           'https://depot.galaxyproject.org/singularity/tool:0.0.0--0' :
    //           'quay.io/biocontainers/tool:0.0.0--0' }"

    // TODO:: Update inputs to allow all possible file inputs to be staged.
    //        When a file is optional, use [] to call the process.
    input:
    // Mandatory
    tuple val(meta), path(reads), path(assembly) // [ [id: 'name'], [file(read1), file(read2)], file(assembly)]
    // Optional
    path(reference)                              // file: reference genome for comparison

    // TODO:: Update outputs for each desired output.
    //        meta should be passed with each output to allow downstream processing
    //        Each output channel should be named using the emit directive.
    output:
    tuple val(meta), path("*"), emit: all

    // TODO:: Update script with command.
    //        Allow tool specific arguments to be passed using ext.args
    //        Allow tool output to be prefixed using $prefix
    script:
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def args   = task.ext.args ?: ''
    """
    inspector.py \\
        $args \\
        -o $prefix \\
        -t $task.cpus \\
        -c $assembly \\
        -r $reads
    """
}
