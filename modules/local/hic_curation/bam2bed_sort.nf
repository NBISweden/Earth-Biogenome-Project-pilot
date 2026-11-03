nextflow.enable.types = true

process BAM2BED_SORT {
    tag "${id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/bedtools_samtools_pip_awk:61be4cfccef27592"

    ext [
        args: { bam2bed_args ?: '' },
        args2: { bam2bed_args2 ?: '' },
        args3: { bam2bed_args3 ?: '1' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        bam: Path,
        bam2bed_args: String?,
        bam2bed_args2: String?,
        bam2bed_args3: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        pairs: file("*.pairs")
    )

    topic:
    tuple(task.process, 'samtools', eval('samtools --version | sed "1!d; s/samtools //"')) >> 'versions'
    tuple(task.process, 'bedtools', eval("bedtools --version | sed -e 's/bedtools v//g'")) >> 'versions'
    tuple(task.process, 'sort', eval("sort --version | sed '1!d; s/.* //'")) >> 'versions'
    tuple(task.process, 'awk', eval("awk --version | sed '1!d; s/mawk //; s/ .*//'")) >> 'versions'
    tuple(task.process, 'paste', eval("paste --version | sed '1!d; s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: '1'
    def prefix = task.ext.prefix ?: "${id}"

    """
    samtools view -@$task.cpus $args ${bam} | \\
    bamToBed | \\
    sort -k4 --parallel=$task.cpus $args2 | \\
    paste -d '\t' - - | \\
    awk -v q=$args3 'BEGIN {FS=\"\t\"; OFS=\"\t\"} { if(int(\$5) >= int(q) && int(\$11) >= int(q)) { if (\$1 > \$7) { print substr(\$4,1,length(\$4)-2),\$7,\$8,\$1,\$2,\$12,\$6,\"UU\"} else { print substr(\$4,1,length(\$4)-2),\$1,\$2,\$7,\$8,\$6,\$12,\"UU\"} } }' | \\
    sort -k2,2V -k4,4V -k3,3n -k5,5n --parallel=$task.cpus $args2 > ${prefix}.pairs
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}.pairs
    """
}
