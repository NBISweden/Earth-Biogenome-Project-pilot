process EXTRACT_CHROMSIZES {
    executor 'local'
    
    input:
    tuple val(meta), val( faidx )

    exec:
    def chr_sizes = faidx.splitCsv( sep: '\t', header: false )
        .collect{ row -> row[ 0..1 ].join('\t') }
        .join('\n')
    file("${task.workDir}/chrom_sizes.tsv").text = chr_sizes

    output:
    tuple val(meta), path("chrom_sizes.tsv"), emit: chr_sizes
}