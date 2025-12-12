process DVPOLISH_CREATE_FINALASM {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::samtools=1.22.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(unpol_fasta), path(unpol_merqury_csv)  // meta map, unpolished assembly, corresponding merqury qv file
    tuple val(meta2), path(pol_fasta), path(pol_merqury_csv)     // meta map, polished assembly, corresponding merqury qv file

    output:
    tuple val(meta), path('*.fasta.gz')     , emit: fasta_gz
    tuple val(meta), path('*_selection.tsv'), emit: dvpolish_selection_tsv
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    # Validate merqury files have same number of lines & are not empty
    nl_unpol_ASM=\$(wc -l < ${unpol_merqury_csv})
    nl_pol_ASM=\$(wc -l < ${pol_merqury_csv})

    if [[ \${nl_unpol_ASM} -ne \${nl_pol_ASM} ]]; then
        >&2 echo "[ERROR] DVPOLISH_CREATE_FINALASM: merqury files have different line counts: ${unpol_merqury_csv}: \${nl_unpol_ASM} != ${pol_merqury_csv}: \${nl_pol_ASM}"
        exit 1
    fi
    if [[ \${nl_unpol_ASM} -eq 0 ]]; then
        >&2 echo "[ERROR] DVPOLISH_CREATE_FINALASM: merqury files are empty"
        exit 1
    fi

    # Create selection report header
    echo -e "contig\\tsource_asm\\tunpolished_asm_errors\\tpolished_asm_errors\\tunpolished_asm_qv\\tpolished_asm_qv" > ${prefix}_selection.tsv

    # Index both assemblies
    samtools faidx ${unpol_fasta}
    samtools faidx ${pol_fasta}

    # Process contigs in single pass (merqury files = file descriptor 0 and 3)
    line_count=0
    while IFS=\$'\\t' read -r -a l_uasm && IFS=\$'\\t' read -r -a l_pasm <&3; do
        line_count=\$((line_count + 1))
        contig="\${l_uasm[0]}"

        # Validate contig name match between unpolished & polished assemblies
        if [[ "\$contig" != "\${l_pasm[0]}" ]]; then
            >&2 echo "[ERROR] DVPOLISH_CREATE_FINALASM: merqury files are not in same order!"
            >&2 echo "[ERROR]        ${unpol_merqury_csv} line \$line_count: \${l_uasm[0]}"
            >&2 echo "[ERROR]        ${pol_merqury_csv} line \$line_count: \${l_pasm[0]}"
            exit 2
        fi

        # Compare error counts (column 2) and extract best contig
        if [[ \${l_uasm[1]} -le \${l_pasm[1]} ]]; then
            samtools faidx ${unpol_fasta} "\$contig"
            source="unpolished"
        else
            samtools faidx ${pol_fasta} "\$contig"
            source="polished"
        fi

        # Append selection info to report
        echo -e "\$contig\\t\$source\\t\${l_uasm[1]}\\t\${l_pasm[1]}\\t\${l_uasm[3]}\\t\${l_pasm[3]}" >> ${prefix}_selection.tsv

    done < ${unpol_merqury_csv} 3< ${pol_merqury_csv} | bgzip -@ ${task.cpus} -c > ${prefix}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """
}
