nextflow.enable.types = true

process SAMTOOLS_FASTA {
    tag "${id}"
    label 'process_low'

    ext [
        args: { samtools_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    conda "bioconda::samtools=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'biocontainers/samtools:1.14--hb421002_0' }"

    input:
    record(
        id: String,
        meta: Record?,
        single_end: Boolean,
        bam: Path,
        samtools_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        fasta: files("*.fasta.gz")
    )

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def prefix    = task.ext.prefix
    // Set endedness to -o to output data from testprofile
    def endedness = single_end ? "-0 ${prefix}.fasta.gz" : "-1 ${prefix}_1.fasta.gz -2 ${prefix}_2.fasta.gz"
    """
    samtools fasta \\
        ${task.ext.args} \\
        --threads $task.cpus \\
        $endedness \\
        $bam
    """

    stub:
    def prefix    = task.ext.prefix
    def endedness = single_end ? "${prefix}.fasta.gz" : "${prefix}_1.fasta.gz ${prefix}_2.fasta.gz"
    """
    touch $endedness
    """
}
