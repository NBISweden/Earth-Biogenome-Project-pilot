process SAMTOOLS_FASTA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'biocontainers/samtools:1.14--hb421002_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    // Set endedness to -o to output data from testprofile
    def endedness = meta.single_end ? "-0 ${prefix}.fasta.gz" : "-1 ${prefix}_1.fasta.gz -2 ${prefix}_2.fasta.gz"
    """
    samtools fasta \\
        $args \\
        --threads $task.cpus \\
        $endedness \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version |& sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
