nextflow.enable.types = true

process TWOREADCOMBINER_FIXMATE_SORT {
    tag "${id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/samtools_perl:e2ce3f7265547d1f"

    ext [
        args: { combiner_args ?: '' },
        args2: { combiner_args2 ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        bam: List<Path>,
        combiner_args: String?,
        combiner_args2: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        bam: file("*.bam")
    )

    topic:
    tuple(task.process, 'samtools', eval('samtools --version | sed "1!d; s/samtools //"')) >> 'versions'
    tuple(task.process, 'perl', eval("perl --version | sed '2!d; s/.*(v//; s/).*//'")) >> 'versions'
    tuple(task.process, 'grep', eval("grep --version | sed '1!d; s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    two_read_bam_combiner_sanger.pl ${bam.join(' ')} samtools ${args} | \\
    grep -v -e "^@HD" | \\
    samtools sort ${args2} -@${task.cpus} -T sort_tmp -o ${prefix}_comb.bam -
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}_comb.bam
    """
}
