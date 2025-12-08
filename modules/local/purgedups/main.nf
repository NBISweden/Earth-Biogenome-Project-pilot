process ITERATIVE_PURGEDUPS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/14ace48b0cddd87c685607131b30552d25486769fc31edb7cc59748969cf5b1a/data':
        'community.wave.seqera.io/library/purge_dups_seqkit:613e35b8286c1a62' }"

    input:
    tuple val(meta), path(reads), path(assemblies)

    output:
    tuple val(meta), path("purged/*.purged.fa")   , emit: purged
    tuple val(meta), path("purged/*.haplotigs.fa"), emit: haplotigs
    tuple val(meta), path("logs/*")               , emit: logs
    path "versions.yml"                           , emit: versions

    script:
    def prefix              = task.ext.prefix ?: "${meta.id}"
    def minimap2_reads_args = task.ext.minimap2_reads_args ?: ''
    def minimap2_self_args  = task.ext.minimap2_self_args  ?: ''
    def pbcstat_args        = task.ext.pbcstat_args        ?: ''
    def calcuts_args        = task.ext.calcuts_args        ?: ''
    def hist_plot_args      = task.ext.hist_plot_args      ?: ''
    def split_fa_args       = task.ext.split_fa_args       ?: ''
    def purge_dups_args     = task.ext.purge_dups_args     ?: ''
    def get_seqs_args       = task.ext.get_seqs_args       ?: ''
    def seqkit_args         = task.ext.seqkit_args         ?: ''
    def reads_input         = reads instanceof List ? reads.join(' ') : reads
    """
    # Create output directories
    mkdir -p purged logs work

    # Get list of assemblies (sorted for consistent ordering)
    mapfile -t ASSEMBLY_FILES < <(printf '%s\\n' ${assemblies.join(' ')} | sort)
    N_HAPS=\${#ASSEMBLY_FILES[@]}
    echo "Processing \${N_HAPS} haplotype(s)"

    PREVIOUS_HAPLOTIGS=""

    # ITERATIVE PURGING: Process each haplotype sequentially
    for HAP_IDX in \$(seq 1 \${N_HAPS}); do
        echo "Processing haplotype \${HAP_IDX} / \${N_HAPS}"
        ARRAY_IDX=\$((HAP_IDX - 1))
        CURRENT_ASM="\${ASSEMBLY_FILES[\${ARRAY_IDX}]}"
        CURRENT_PREFIX="${prefix}_hap\${HAP_IDX}"

        # Merge assembly with previous haplotigs
        if [ -n "\${PREVIOUS_HAPLOTIGS}" ] && [ -f "\${PREVIOUS_HAPLOTIGS}" ]; then
            echo "Merging with previous haplotigs: \${PREVIOUS_HAPLOTIGS}"
            seqkit seq "\${CURRENT_ASM}" "\${PREVIOUS_HAPLOTIGS}" > "work/\${CURRENT_PREFIX}_merged.fa"
            TARGET_ASM="work/\${CURRENT_PREFIX}_merged.fa"
        else
            echo "First haplotype - no previous haplotigs to merge"
            TARGET_ASM="\${CURRENT_ASM}"
        fi

        echo "Aligning reads to assembly"
        minimap2 \\
            -t ${task.cpus} \\
            ${minimap2_reads_args} \\
            "\${TARGET_ASM}" \\
            ${reads_input} \\
            > "work/\${CURRENT_PREFIX}_reads.paf"

        echo "Calculating coverage statistics"
        pbcstat \\
            ${pbcstat_args} \\
            "work/\${CURRENT_PREFIX}_reads.paf"
        for PBFILE in PB.*; do
            mv \$PBFILE "work/\${CURRENT_PREFIX}.\${PBFILE}"
        done

        calcuts \\
            ${calcuts_args} \\
            "work/\${CURRENT_PREFIX}.PB.stat" \\
            > "work/\${CURRENT_PREFIX}.cutoffs" \\
            2>| >(tee work/\${CURRENT_PREFIX}.calcuts.log >&2)

        hist_plot.py \\
            ${hist_plot_args} \\
            -c "work/\${CURRENT_PREFIX}.cutoffs" \\
            "work/\${CURRENT_PREFIX}.PB.stat" \\
            "logs/\${CURRENT_PREFIX}.hist_plot.png"

        echo "Splitting assembly"
        split_fa \\
            ${split_fa_args} \\
            "\${TARGET_ASM}" \\
            > "work/\${CURRENT_PREFIX}.split.fa"

        echo "Self-aligning assembly"
        minimap2 \\
            -t ${task.cpus} \\
            ${minimap2_self_args} \\
            "work/\${CURRENT_PREFIX}.split.fa" \\
            "work/\${CURRENT_PREFIX}.split.fa" \\
            > "work/\${CURRENT_PREFIX}_self.paf"

        echo "Identifying duplicates"
        purge_dups \\
            ${purge_dups_args} \\
            -T "work/\${CURRENT_PREFIX}.cutoffs" \\
            -c "work/\${CURRENT_PREFIX}.PB.base.cov" \\
            "work/\${CURRENT_PREFIX}_self.paf" \\
            > "logs/\${CURRENT_PREFIX}.dups.bed" \\
            2>| >(tee "work/\${CURRENT_PREFIX}_purge.log" >&2)

        echo "Extracting sequences"
        get_seqs \\
            ${get_seqs_args} \\
            -p "\${CURRENT_PREFIX}" \\
            "logs/\${CURRENT_PREFIX}.dups.bed" \\
            "\${TARGET_ASM}"

        seqkit seq \\
            ${seqkit_args} \\
            "\${CURRENT_PREFIX}.purged.fa" \\
            > "purged/\${CURRENT_PREFIX}.purged.fa"
        seqkit seq \\
            ${seqkit_args} \\
            "\${CURRENT_PREFIX}.hap.fa" \\
            > "purged/\${CURRENT_PREFIX}.haplotigs.fa"
        PREVIOUS_HAPLOTIGS="purged/\${CURRENT_PREFIX}.haplotigs.fa"
        echo "Completed haplotype \${HAP_IDX}"
    done

    if [ "\${N_HAPS}" -eq 1 ]; then
        # Copy haplotigs as hap2 for Merqury
        cp purged/${prefix}_hap1.haplotigs.fa purged/${prefix}_hap2.purged.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        purge_dups: \$(purge_dups -h |& sed '3!d; s/.*: //')
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}