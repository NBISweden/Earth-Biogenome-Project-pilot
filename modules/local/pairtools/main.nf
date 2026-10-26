nextflow.enable.types = true

process PAIRTOOLS {
    tag "${id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.1.0--py39hd5a99d8_2' :
        'biocontainers/pairtools:1.1.0--py39hd5a99d8_2' }"

    ext [
        args: { pairtools_args ?: '' },
        args2: { pairtools_args2 ?: '' },
        args3: { pairtools_args3 ?: '' },
        args4: { pairtools_args4 ?: '' },
        args5: { pairtools_args5 ?: '' },
        args6: { pairtools_args6 ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        bam: List<Path>,
        chromsizes: Path,
        fasta: Path,
        sort_bam: Boolean,
        pairtools_args: String?,
        pairtools_args2: String?,
        pairtools_args3: String?,
        pairtools_args4: String?,
        pairtools_args5: String?,
        pairtools_args6: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        pairs: file("*.split.pairs.gz"),
        bam: file("*.split.pairs.bam", optional: true),
        cram: file("*.split.pairs.cram", optional: true),
        crai: file("*.crai", optional: true),
        csi: file("*.csi", optional: true),
        stat: files("*.stat")
    )

    topic:
    tuple(task.process, 'samtools', eval('samtools --version | sed "1!d; s/samtools //"')) >> 'versions'
    tuple(task.process, 'pairtools', eval("pairtools --version 2>&1 | sed 's/pairtools, version //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def args5 = task.ext.args5 ?: ''
    def args6 = task.ext.args6 ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    def halfMemory = task.memory.toGiga().intdiv(2)
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args6 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    """
    export MPLCONFIGDIR=tmp

    BAM_FILES=(${bam.join(' ')})
    TEMP_FILES=()

    for BAM in "\${BAM_FILES[@]}"; do
        BAM_PREFIX=\$(basename "\$BAM" .bam)
        TEMP_FILE="\$BAM_PREFIX.pairs_temp.gz"
        TEMP_FILES+=("\$TEMP_FILE")
        pairtools parse \\
            -c ${chromsizes} \\
            ${args} \\
            --output-stats \$BAM_PREFIX.pairsam.stat \\
            \$BAM | \\
        pairtools sort \\
            ${args2} \\
            --nproc ${task.cpus} \\
            --memory ${halfMemory}G \\
            --output \$TEMP_FILE
    done

    pairtools merge \\
        ${args3} \\
        --nproc ${task.cpus} \\
        --memory ${halfMemory}G \\
        "\${TEMP_FILES[@]}" | \\
    pairtools dedup \\
        ${args4} \\
        --output-stats ${prefix}_dedup.pairs.stat | \\
    pairtools split \\
        --nproc-in ${task.cpus} \\
        --nproc-out ${task.cpus} \\
        --output-pairs ${prefix}.split.pairs.gz \\
        --output-sam - \\
        ${args5} | \\
    samtools ${samtools_command} \\
        ${args6} \\
        -@ ${task.cpus} \\
        ${reference} \\
        -o ${prefix}.split.pairs.${extension} \\
        -

    rm -f "\${TEMP_FILES[@]}"
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}.split.pairs.gz
    touch ${prefix}.split.pairs.bam
    touch ${prefix}.pairsam.stat
    touch ${prefix}_dedup.pairs.stat
    """
}
