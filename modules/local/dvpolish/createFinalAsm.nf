process DVPOLISH_CREATE_FINALASM {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::seqkit=2.8.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.2--h9ee0642_0' :
        'nf-core/seqkit:2.8.2--h9ee0642_0' }"

    input:
    tuple val(meta), path(unpol_fasta), path(unpol_merqury_cxv)  // meta map, unpolished assembly, corresponding merqury qv file
    tuple val(meta2), path(pol_fasta), path(pol_merqury_csv)  // meta map, polished assembly, corresponding merqury qv file

    output:
    tuple val(meta), path('*.fa.gz'), emit: fa_gz
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def unpol_name = unpol_fasta.getSimpleName()
    def unpol_ext  = unpol_fasta.getExtension()
    def pol_name   = pol_fasta.getSimpleName()
    def pol_ext    = pol_fasta.getExtension()

    """
    nl_unpol_ASM=$(wc -l < ${unpol_merqury_cxv})
    nl_pol_ASM=$(wc -l < ${pol_merqury_csv})

    if [[ ${nl_unpol_ASM} -ne ${nl_pol_ASM} ]]
    then 
        >&2 echo "[ERROR] DVPOLISH_CREATE_FINALASM: merqury files have different lines: ${unpol_merqury_cxv}: ${${nl_unpol_ASM}} != ${pol_merqury_csv}: ${${nl_pol_ASM}}"
        exit 1
    fi
    
    #split unpolished assembly by sequence ID
    seqkit split -i -O unpolished_asm ${unpol_fasta}

    #split polished assembly by sequence ID
    seqkit split -i -O polished_asm ${pol_fasta}

    l=1 
    while [[ $l -le ${nl_pol_ASM} ]]
    do 
        IFS='\t' read -r -a l_uasm <<< "$(sed -n ${l}p ${pol_merqury_cxv})"
        IFS='\t' read -r -a l_pasm <<< "$(sed -n ${l}p ${unpol_merqury_cxv})"

        // check if the contig names (column 1) are the same
        if [[ "${l_uasm[0]}" != "${l_pasm[0]}" ]]
        then 
            >&2 echo "[ERROR] DVPOLISH_CREATE_FINALASM: merqury files are not in the same order!"
            >&2 echo "[ERROR]        file:  ${unpol_fasta} line $l: ${l_uasm[*]}"
            >&2 echo "[ERROR]        file:  ${pol_fasta}   line $l: ${l_pasm[*]}"
            exit 2
        fi 

        // compare number of errorneous kmers (column 2)
        if [[ ${l_uasm[1]} -le ${l_pasm[1]} ]]    // unpolished assembly has fewer errors, or no difference -> go with the unpolished assembly 
        then 
            cat unpolished_asm/${unpol_name}_part_${l_uasm[0]}.${unpol_ext}
            >&2 echo "[WARNING] DVPOLISH_CREATE_FINALASM: unpolished contig ${l_uasm[0]} has better or equal QV: ${l_uasm[3]} vs ${l_pasm[3]}"
        else                                      // polished assembly has fewer errors
            cat polished_asm/${pol_name}_part_${l_pasm[0]}.${pol_ext}
        fi | gzip -c > ${prefix}.fa.gz

        l=$((l+1))
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dvpolish: \$(seqkit version)
    END_VERSIONS
    """
}
