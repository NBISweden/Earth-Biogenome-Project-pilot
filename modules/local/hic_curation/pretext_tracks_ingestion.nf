nextflow.enable.types = true

process PRETEXT_TRACKS_INGESTION {
    tag "${id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/pretextgraph:0.0.6--561a906eedf096a7"

    ext [
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        pretext_in: Path,
        cov_bedgraph: Path,
        telomer_bedgraph: Path?,
        gap_bedgraph: Path?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        hitile: file("*_wTracks.pretext")
    )

    topic:
    tuple(task.process, 'PretextGraph', eval("PretextGraph | sed '2!d; s/.* //'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${id}"
    """
    cat ${cov_bedgraph} | PretextGraph -i ${pretext_in} -n coverage -o ${prefix}_cov.pretext
    if [[ ${telomer_bedgraph ? telomer_bedgraph.size() : 0} -gt 0 && ${gap_bedgraph ? gap_bedgraph.size() : 0} -gt 0  ]]
    then
        cat ${telomer_bedgraph} | PretextGraph -i ${prefix}_cov.pretext -n telomer -o ${prefix}_wTracks.pretext
        cat ${gap_bedgraph} | PretextGraph -i ${prefix}_wTracks.pretext -n gap
    elif [[ ${telomer_bedgraph ? telomer_bedgraph.size() : 0} -eq 0 && ${gap_bedgraph ? gap_bedgraph.size() : 0} -gt 0  ]]
    then
        cat ${gap_bedgraph} | PretextGraph -i ${prefix}_cov.pretext -n gap -o ${prefix}_wTracks.pretext
    elif [[ ${telomer_bedgraph ? telomer_bedgraph.size() : 0} -gt 0 && ${gap_bedgraph ? gap_bedgraph.size() : 0} -eq 0  ]]
    then
        cat ${telomer_bedgraph} | PretextGraph -i ${prefix}_cov.pretext -n telomer -o ${prefix}_wTracks.pretext
    else
        mv ${prefix}_cov.pretext ${prefix}_wTracks.pretext
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}_wTracks.pretext
    """
}
