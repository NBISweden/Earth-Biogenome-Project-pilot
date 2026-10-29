nextflow.enable.types = true

process MERQURYFK_MERQURYFK {
    tag "${id}"
    label 'process_medium'

    //conda "${moduleDir}/environment.yml"
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    ext [
        args: { merquryfk_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        fastk_hist: Path,
        fastk_ktab: List<Path>,
        assembly: Path,
        haplotigs: Path?,
        mathaptab: List<Path>,
        pathaptab: List<Path>,
        merquryfk_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        stats: file("${task.ext.prefix}.completeness.stats"),
        bed: files("${task.ext.prefix}.*_only.bed"),
        assembly_qv: files("${task.ext.prefix}.*.qv"),
        part_spectra_cn_fl: files("${task.ext.prefix}.*.spectra-cn.fl.{png,pdf}"),
        part_spectra_cn_ln: files("${task.ext.prefix}.*.spectra-cn.ln.{png,pdf}"),
        part_spectra_cn_st: files("${task.ext.prefix}.*.spectra-cn.st.{png,pdf}"),
        full_spectra_cn_fl: file("${task.ext.prefix}.spectra-cn.fl.{png,pdf}", optional: true),
        full_spectra_cn_ln: file("${task.ext.prefix}.spectra-cn.ln.{png,pdf}", optional: true),
        full_spectra_cn_st: file("${task.ext.prefix}.spectra-cn.st.{png,pdf}", optional: true),
        qv: file("${task.ext.prefix}.qv"),
        spectra_asm_fl: file("${task.ext.prefix}.spectra-asm.fl.{png,pdf}", optional: true),
        spectra_asm_ln: file("${task.ext.prefix}.spectra-asm.ln.{png,pdf}", optional: true),
        spectra_asm_st: file("${task.ext.prefix}.spectra-asm.st.{png,pdf}", optional: true),
        phased_block_bed: file("${task.ext.prefix}.phased_block.bed", optional: true),
        phased_block_stats: file("${task.ext.prefix}.phased_block.stats", optional: true),
        continuity_N: file("${task.ext.prefix}.continuity.N.{pdf,png}", optional: true),
        block_N: file("${task.ext.prefix}.block.N.{pdf,png}", optional: true),
        block_blob: file("${task.ext.prefix}.block.blob.{pdf,png}", optional: true),
        hapmers_blob: file("${task.ext.prefix}.hapmers.blob.{pdf,png}", optional: true),
        false_duplications: file("${task.ext.prefix}.false_duplications.tsv"),
        cn_histogram: file("${task.ext.prefix}.spectra-cn.cni.gz")
    )

    topic:
    tuple(task.process, 'merquryfk', '41451f8fb146158c5b747ae7915e69975c61ddd9') >> 'versions'
    tuple(task.process, 'fastk', '1.2') >> 'versions'
    tuple(task.process, 'R', eval('R --version | sed "1!d; s/.*version //; s/ .*//"')) >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MERQURYFK_MERQURYFK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    if ((mathaptab || pathaptab) && !(mathaptab && pathaptab)) {
        log.error("Error: Only one of the maternal and paternal hap tabs have been provided!")
    }

    def args        = task.ext.args
    def prefix      = task.ext.prefix
    def fk_ktab     = fastk_ktab ? "${fastk_ktab.find { path -> path.name.endsWith(".ktab") }}" : ''
    def mat_hapktab = mathaptab  ? "${mathaptab.find { path -> path.name.endsWith(".ktab") }}"  : ''
    def pat_hapktab = pathaptab  ? "${pathaptab.find { path -> path.name.endsWith(".ktab") }}"  : ''
    """
    MerquryFK \\
        $args \\
        -T$task.cpus \\
        ${fk_ktab} \\
        ${mat_hapktab} \\
        ${pat_hapktab} \\
        $assembly \\
        ${haplotigs ? "${haplotigs}" : ''} \\
        $prefix

    # Providing 1 assembly outputs 1 cni that requires a rename
    # Providing assembly + haplotigs outputs 3 cnis, of which one is correctly named already
    CNI_COUNT=\$(find . -maxdepth 1 -name "*.spectra-cn.cni" -type f | wc -l)
    if [ "\$CNI_COUNT" -eq 1 ]; then
      mv ./*.spectra-cn.cni ${prefix}.spectra-cn.cni
    fi

    awk -v asm_ploidy=1 \\
        -f $projectDir/bin/false_duplications.awk ${prefix}.spectra-cn.cni \\
        > ${prefix}.false_duplications.tsv

    gzip ${prefix}.spectra-cn.cni
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}.completeness.stats
    touch ${prefix}.qv
    touch ${prefix}._.qv
    touch ${prefix}._only.bed
    touch ${prefix}.false_duplications.tsv
    echo "" | gzip > ${prefix}.spectra-cn.cni.gz
    """
}
