nextflow.enable.types = true

process DVPOLISH_PBMM2_ALIGN {
    tag "${id}"
    label 'process_medium'

    ext [
        args: { pbmm2_args ?: '' }
    ]

    // Note: the versions here need to match the versions used in pbmm2/index
    conda 'bioconda::pbmm2=1.13.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmm2:1.13.1--h9ee0642_0' :
        'biocontainers/pbmm2:1.13.1--h9ee0642_0' }"

    input:
    record(
        id: String,
        meta: Record?,
        reads: Path,
        index: Path,
        pbmm2_args: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        bam: file("*.bam"),
        bai: file("*.bai")
    )

    topic:
    tuple(task.process, 'pbmm2', eval("pbmm2 --version 2>&1 | sed '1!d; s/pbmm2 //'")) >> 'versions'

    script:
    def out_name_part1 = index.baseName
    def out_name_part2 = reads.name.replaceAll(/\.(fastq|fq|fasta|fa)(\.gz)?$/, '')
    """
    pbmm2 align --sort \\
        ${task.ext.args} \\
        -j $task.cpus \\
        "$index" \\
        "$reads" \\
        ${out_name_part1}_${out_name_part2}.bam
    """

    stub:
    def out_name = "${index.baseName}_${reads.name.replaceAll(/\.(fastq|fq|fasta|fa)(\.gz)?$/, '')}"
    """
    touch ${out_name}.bam
    touch ${out_name}.bam.bai
    """
}
