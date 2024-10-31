process PRETEXT_TRACKS_INGESTION {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/pretextgraph:0.0.6--561a906eedf096a7"

    input:
    tuple val(meta), path(pretext_in), path(cov_bedgraph), path(telomer_bedgraph), path(gap_bedgraph)

    output:
    tuple val(meta), path("*_wTracks.pretext") , emit: hitile 
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    cat ${cov_bedgraph} | PretextGraph -i ${pretext_in} -n coverage -o ${prefix}_cov.pretext
    if [[ ${telomer_bedgraph.size()} -gt 0 && ${gap_bedgraph.size()} -gt 0  ]]
    then 
        cat ${telomer_bedgraph} | PretextGraph -i ${prefix}_cov.pretext -n telomer -o ${prefix}_wTracks.pretext
        cat ${gap_bedgraph} | PretextGraph -i ${prefix}_wTracks.pretext -n gap
    elif [[ ${telomer_bedgraph.size()} -eq 0 && ${gap_bedgraph.size()} -gt 0  ]]
    then 
        cat ${gap_bedgraph} | PretextGraph -i ${prefix}_cov.pretext -n gap -o ${prefix}_wTracks.pretext
    elif [[ ${telomer_bedgraph.size()} -gt 0 && ${gap_bedgraph.size()} -eq 0  ]]
    then 
        cat ${telomer_bedgraph} | PretextGraph -i ${prefix}_cov.pretext -n telomer -o ${prefix}_wTracks.pretext
    else 
        mv ${prefix}_cov.pretext ${prefix}_wTracks.pretext
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextGraph:     \$(echo \$(PretextGraph 2>&1 | grep Version) | awk '{print \$NF}')
    END_VERSIONS
    """
}
