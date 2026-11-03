nextflow.enable.types = true

process BAM2COVERAGE_TRACKS {
    tag "${id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/bedtools_samtools_ucsc-bedgraphtobigwig:ba8edacb0d7bff96"

    ext [
        args: { coverage_args ?: '' },
        args2: { coverage_args2 ?: '' },
        args3: { coverage_args3 ?: '' },
        args4: { coverage_args4 ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        bam: Path,
        chrom_sizes: Path,
        hifi_coverage_cap: Integer,
        coverage_args: String?,
        coverage_args2: String?,
        coverage_args3: String?,
        coverage_args4: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        bed: file("*_coverage.bed"),
        bigwig: file("*_coverage.bw"),
        capped_bed: file("*_coverage_capped.bed"),
        capped_bigwig: file("*_coverage_capped.bw")
    )

    topic:
    tuple(task.process, 'samtools', eval('samtools --version | sed "1!d; s/samtools //"')) >> 'versions'
    tuple(task.process, 'bedtools', eval("bedtools --version | sed -e 's/bedtools v//g'")) >> 'versions'
    tuple(task.process, 'sort', eval("sort --version | sed '1!d; s/.* //'")) >> 'versions'
    tuple(task.process, 'awk', eval("awk --version | sed '1!d; s/mawk //; s/ .*//'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    """
    # convert bam to sorted bed
    samtools sort -@ $task.cpus $args ${bam} | \\
    samtools view -@ $task.cpus $args2 - | \\
    genomeCoverageBed $args3 -ibam - | \\
    sort -k1,1 -k2,2n --parallel=${task.cpus} $args4  > ${prefix}_coverage.bed
    # convert bed to bigwig
    bedGraphToBigWig ${prefix}_coverage.bed ${chrom_sizes} ${prefix}_coverage.bw
    # coverage cap
    sort -k1,1V -k2,2n -k3,3n $args4 ${prefix}_coverage.bed | \\
    awk -v COVERAGE_CAP=$hifi_coverage_cap 'BEGIN {FS = " "; OFS = "\t";}{if(int(\$4)>COVERAGE_CAP){print \$1,\$2,\$3,COVERAGE_CAP} else {print \$1,\$2,\$3,\$4}}' > ${prefix}_coverage_capped.bed
    # convert bed to bigwig
    bedGraphToBigWig ${prefix}_coverage_capped.bed ${chrom_sizes} ${prefix}_coverage_capped.bw
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}_coverage.bed ${prefix}_coverage.bw
    touch ${prefix}_coverage_capped.bed ${prefix}_coverage_capped.bw
    """
}
