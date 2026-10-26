nextflow.enable.types = true

process EXTRACT_CHROMSIZES {
    tag "${id}"
    executor 'local'

    input:
    record(
        id: String,
        meta: Record?,
        faidx: Path
    )

    output:
    record(
        id: id,
        meta: meta,
        chr_sizes: file("chrom_sizes.tsv")
    )

    exec:
    def chr_sizes = faidx.splitCsv( sep: '\t', header: false )
        .collect{ row -> row[ 0..1 ].join('\t') }
        .join('\n')
    file("${task.workDir}/chrom_sizes.tsv").text = chr_sizes
}
