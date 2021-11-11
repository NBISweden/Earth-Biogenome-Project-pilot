process GETHIFI {

    conda "${task.ext.enable_conda ? 'bioconda::samtools==1.14-0' : '' }"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/samtools%3A1.14--hb421002_0' :
              'quay.io/biocontainers/samtools:1.14--hb421002_0' }"

    input:
    // Mandatory
    tuple val(meta), path(bamfile) // [ [id: 'name'], file(bamfile)]
    
    output:
    tuple val(meta), path('*.fasta') , emit: fasta // reads into a fasta file
    
    script:
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def args   = task.ext.args ?: ''

    """
    samtools fasta '$bamfile' > "${prefix}.fasta"
    """
}
