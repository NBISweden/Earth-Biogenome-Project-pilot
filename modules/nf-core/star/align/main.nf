nextflow.enable.types = true

process STAR_ALIGN {
    tag "${id}"
    label 'process_high'

    conda "bioconda::star=2.7.10a bioconda::samtools=1.16.1 conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    ext [
        args: { star_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        single_end: Boolean,
        reads: List<Path>,
        index: Path,
        gtf: Path?,
        star_ignore_sjdbgtf: Boolean,
        seq_platform: String?,
        seq_center: String?,
        star_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        log_final: file("*Log.final.out"),
        log_out: file("*Log.out"),
        log_progress: file("*Log.progress.out"),
        bam: file('*d.out.bam', optional: true),
        bam_sorted: file('*sortedByCoord.out.bam', optional: true),
        bam_transcript: file('*toTranscriptome.out.bam', optional: true),
        bam_unsorted: file('*Aligned.unsort.out.bam', optional: true),
        fastq: files('*fastq.gz', optional: true),
        tab: files('*.tab', optional: true),
        spl_junc_tab: file('*.SJ.out.tab', optional: true),
        read_per_gene_tab: file('*.ReadsPerGene.out.tab', optional: true),
        junction: file('*.out.junction', optional: true),
        sam: file('*.out.sam', optional: true),
        wig: files('*.wig', optional: true),
        bedgraph: files('*.bg', optional: true)
    )

    topic:
    tuple(task.process, 'star', eval('STAR --version | sed "s/STAR_//"')) >> 'versions'
    tuple(task.process, 'samtools', eval("samtools --version | sed -n '1s/samtools //p'")) >> 'versions'
    tuple(task.process, 'gawk', eval("gawk --version | sed -n '1s/GNU Awk \\([0-9.]*\\).*/\\1/p'")) >> 'versions'

    script:
    def args = task.ext.args
    def prefix = task.ext.prefix
    def n = reads.size()
    def reads1 = single_end ? reads : (0..<n).findAll { i -> i % 2 == 0 }.collect { i -> reads[i] }.toList()
    def reads2 = single_end ? [] : (0..<n).findAll { i -> i % 2 == 1 }.collect { i -> reads[i] }.toList()
    def ignore_gtf      = star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"
    def seq_platform_   = seq_platform ? "'PL:$seq_platform'" : ""
    def seq_center_     = seq_center ? "'CN:$seq_center'" : ""
    def attrRG          = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:$prefix' $seq_center_ 'SM:$prefix' $seq_platform_"
    def out_sam_type    = (args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'
    def mv_unsorted_bam = (args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $out_sam_type \\
        $ignore_gtf \\
        $attrRG \\
        $args

    $mv_unsorted_bam

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}Xd.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.sortedByCoord.out.bam
    touch ${prefix}.toTranscriptome.out.bam
    touch ${prefix}.Aligned.unsort.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.unmapped_1.fastq.gz
    touch ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}.tab
    touch ${prefix}.SJ.out.tab
    touch ${prefix}.ReadsPerGene.out.tab
    touch ${prefix}.Chimeric.out.junction
    touch ${prefix}.out.sam
    touch ${prefix}.Signal.UniqueMultiple.str1.out.wig
    touch ${prefix}.Signal.UniqueMultiple.str1.out.bg
    """
}
