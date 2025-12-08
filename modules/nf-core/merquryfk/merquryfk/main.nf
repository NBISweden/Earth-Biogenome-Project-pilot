process MERQURYFK_MERQURYFK {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    input:
    tuple val(meta), path(fastk_hist),path(fastk_ktab),path(assembly),path(haplotigs)
    path matktab                                                                        //optional
    path patktab                                                                        //optional

    output:
    tuple val(meta), path("${prefix}.completeness.stats")         , emit: stats
    tuple val(meta), path("${prefix}.*_only.bed")                 , emit: bed
    tuple val(meta), path("${prefix}.*.qv")                       , emit: assembly_qv
    tuple val(meta), path("${prefix}.*.spectra-cn.fl.{png,pdf}")  , emit: part_spectra_cn_fl,  optional: true
    tuple val(meta), path("${prefix}.*.spectra-cn.ln.{png,pdf}")  , emit: part_spectra_cn_ln,  optional: true
    tuple val(meta), path("${prefix}.*.spectra-cn.st.{png,pdf}")  , emit: part_spectra_cn_st,  optional: true
    tuple val(meta), path("${prefix}.spectra-cn.fl.{png,pdf}")    , emit: full_spectra_cn_fl,  optional: true
    tuple val(meta), path("${prefix}.spectra-cn.ln.{png,pdf}")    , emit: full_spectra_cn_ln,  optional: true
    tuple val(meta), path("${prefix}.spectra-cn.st.{png,pdf}")    , emit: full_spectra_cn_st,  optional: true
    tuple val(meta), path("${prefix}.qv")                         , emit: qv
    tuple val(meta), path("${prefix}.spectra-asm.fl.{png,pdf}")   , emit: spectra_asm_fl,      optional: true
    tuple val(meta), path("${prefix}.spectra-asm.ln.{png,pdf}")   , emit: spectra_asm_ln,      optional: true
    tuple val(meta), path("${prefix}.spectra-asm.st.{png,pdf}")   , emit: spectra_asm_st,      optional: true
    tuple val(meta), path("${prefix}.phased_block.bed")           , emit: phased_block_bed,    optional: true
    tuple val(meta), path("${prefix}.phased_block.stats")         , emit: phased_block_stats,  optional: true
    tuple val(meta), path("${prefix}.continuity.N.{pdf,png}")     , emit: continuity_N,        optional: true
    tuple val(meta), path("${prefix}.block.N.{pdf,png}")          , emit: block_N,             optional: true
    tuple val(meta), path("${prefix}.block.blob.{pdf,png}")       , emit: block_blob,          optional: true
    tuple val(meta), path("${prefix}.hapmers.blob.{pdf,png}")     , emit: hapmers_blob,        optional: true
    tuple val(meta), path("${prefix}.false_duplications.tsv")     , emit: false_duplications
    tuple val(meta), path("${prefix}.spectra-cn.cni.gz")          , emit: cn_histogram
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MERQURYFK_MERQURYFK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def mat_ktab = matktab ? "${matktab.find{ it.toString().endsWith(".ktab") }}" : ''
    def pat_ktab = patktab ? "${patktab.find{ it.toString().endsWith(".ktab") }}" : ''
    def FASTK_VERSION = '0e24fb45b71c4e14382ae1e1bc063bf66ea4e112' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '41451f8fb146158c5b747ae7915e69975c61ddd9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    MerquryFK \\
        $args \\
        -T$task.cpus \\
        ${fastk_ktab.find{ it.toString().endsWith(".ktab") }} \\
        ${mat_ktab} \\
        ${pat_ktab} \\
        $assembly \\
        $haplotigs \\
        $prefix

    # Providing 1 assembly outputs 1 cni that requires a rename
    # Providing assembly + haplotigs outputs 3 cnis, of which one is correctly named already
    CNI_COUNT=\$(find . -maxdepth 1 -name "*.spectra-cn.cni" -type f | wc -l)
    if [ "\$CNI_COUNT" -eq 1 ]; then
      mv ./*.spectra-cn.cni ${prefix}.spectra-cn.cni
    fi

    awk -v asm_ploidy=${assembly instanceof List ? assembly.size() : 1} \\
        -f $projectDir/bin/false_duplications.awk ${prefix}.spectra-cn.cni \\
        > ${prefix}.false_duplications.tsv

    gzip ${prefix}.spectra-cn.cni

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '0e24fb45b71c4e14382ae1e1bc063bf66ea4e112' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '41451f8fb146158c5b747ae7915e69975c61ddd9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.completeness.stats
    touch ${prefix}.qv
    touch ${prefix}._.qv
    touch ${prefix}._only.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """
}
