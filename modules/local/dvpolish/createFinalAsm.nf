process DVPOLISH_CREATE_FINALASM {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::seqkit=2.8.2 bioconda::samtools=1.22.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools_seqkit:d4c7f2ace72c42af' :
        'community.wave.seqera.io/library/samtools_seqkit:27168e18ad389b27' }"

    input:
    tuple val(meta), path(unpol_fasta), path(unpol_merqury_csv)  // meta map, unpolished assembly, corresponding merqury qv file
    tuple val(meta2), path(pol_fasta), path(pol_merqury_csv)     // meta map, polished assembly, corresponding merqury qv file

    output:
    tuple val(meta), path('*.fasta.gz')  , emit: fasta_gz
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def unpol_name = unpol_fasta.getBaseName()
    def unpol_ext  = unpol_fasta.getExtension()
    def pol_name   = pol_fasta.getBaseName()
    def pol_ext    = pol_fasta.getExtension()

    """
    nl_unpol_ASM=\$(wc -l < ${unpol_merqury_csv})
    nl_pol_ASM=\$(wc -l < ${pol_merqury_csv})

    if [[ \${nl_unpol_ASM} -ne \${nl_pol_ASM} ]]
    then 
        >&2 echo "[ERROR] DVPOLISH_CREATE_FINALASM: merqury files have different line counts: ${unpol_merqury_csv}: \${nl_unpol_ASM} != ${pol_merqury_csv}: \${nl_pol_ASM}"
        exit 1
    fi

    if [[ \${nl_unpol_ASM} -eq 0 ]]; then
        >&2 echo "[ERROR] DVPOLISH_CREATE_FINALASM: merqury files have no lines"
        exit 1
    fi

    #split unpolished assembly by sequence ID
    seqkit split -i -O unpolished_asm ${unpol_fasta}

    #split polished assembly by sequence ID
    seqkit split -i -O polished_asm ${pol_fasta}

    l=1
    while [[ \$l -le \${nl_pol_ASM} ]]
    do 
        # read qv string in line p from unpolished and polished merqury files into bash array l_uasm and l_pasm respectively
        IFS='\t' read -r -a l_uasm <<< "\$(sed -n \${l}p ${unpol_merqury_csv})"
        IFS='\t' read -r -a l_pasm <<< "\$(sed -n \${l}p ${pol_merqury_csv})"

        # check if the contig names (column 1) are the same
        if [[ "\${l_uasm[0]}" != "\${l_pasm[0]}" ]]
        then 
            >&2 echo "[ERROR] DVPOLISH_CREATE_FINALASM: merqury files are not in the same order!"
            >&2 echo "[ERROR]        file:  ${unpol_fasta} line \$l: \${l_uasm[*]}"
            >&2 echo "[ERROR]        file:  ${pol_fasta}   line \$l: \${l_pasm[*]}"
            exit 2
        fi 

        # compare number of erroneous kmers (column 2)
        if [[ \${l_uasm[1]} -le \${l_pasm[1]} ]]    # unpolished assembly has fewer errors, or no difference -> go with the unpolished assembly
        then
            cat unpolished_asm/${unpol_name}.part_\${l_uasm[0]}.${unpol_ext}
            >&2 echo "[WARNING] DVPOLISH_CREATE_FINALASM: unpolished contig \${l_uasm[0]} has better or equal QV: \${l_uasm[3]} vs \${l_pasm[3]}"
        else                                      # polished assembly has fewer errors
            cat polished_asm/${pol_name}.part_\${l_pasm[0]}.${pol_ext}
        fi

        l=\$((l+1))
    done | bgzip -c > ${prefix}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/^seqkit v//')
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """
}
