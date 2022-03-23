process INSPECTOR {

    tag "$assembly.simpleName"
    label 'process_medium'
    // TODO:: Update conda path to correct conda depedencies
    // conda "${params.enable_conda ? 'bioconda::tool=0.0.0' : '' }"
    // TODO:: Update Singularity and Docker paths to correct container paths
    container 'ghcr.io/nbisweden/earth-biogenome-project-pilot/inspector:1.0'
    // container "${workflow.containerEngine == 'singularity' &&
    //               !task.ext.singularity_pull_docker_container ?
    //           'https://depot.galaxyproject.org/singularity/tool:0.0.0--0' :
    //           'quay.io/biocontainers/tool:0.0.0--0' }"

    input:
    // Mandatory
    tuple val(meta), path(reads), path(assembly) // [ [id: 'name'], [file(read1), file(read2)], file(assembly)]
    // Optional
    path(reference)                              // file: reference genome for comparison

    output:
    tuple val(meta), path("$prefix/Inspector.log")          , emit: log
    tuple val(meta), path("$prefix/summary_statistics")     , emit: stats
    tuple val(meta), path("$prefix/valid_contig.{fa,fai}")  , emit: valid_contigs
    tuple val(meta), path("$prefix/small_scale_error.bed")  , emit: small_error_bed
    tuple val(meta), path("$prefix/structural_error.bed")   , emit: struct_error_bed
    // tuple val(meta), path("$prefix/read_to_contig.{bam,bai}"), emit: read_to_contigs // Omitting for time being
    path "versions.yml"                                     , emit: versions

    script:
    prefix   = task.ext.prefix ?: meta.id
    def args = task.ext.args   ?: ''
    def ref  = reference ? "--ref $reference" : ''
    """
    inspector.py \\
        $args \\
        -o $prefix \\
        -t $task.cpus \\
        -c $assembly \\
        -r $reads \\
        $ref

    cat <<-END_VERSIONS > versions.yml
    '${task.process}':
        inspector: \$( inspector.py --version |& sed 's/Inspector_v//' )
    END_VERSIONS
    """
}
