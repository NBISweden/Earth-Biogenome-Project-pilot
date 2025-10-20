process MITOHIFI_FINDMITOREFERENCE {
    tag "$species"
    label 'process_low'

    // Docker image available at the project github repository
    container 'ghcr.io/marcelauliano/mitohifi:master'

    input:
    tuple val(meta), val(species)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta, optional: true
    tuple val(meta), path("*.gb")   , emit: gb, optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set +e
    findMitoReference.py \\
        --species "$species" \\
        --outfolder . \\
        $args
    set -e

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py -v | sed 's/.* //' )
    END_VERSIONS

    # Test: task success if both, or no outputs are found, and fail upon partial output (one of two filetypes found)

    FASTA_COUNT=\$(ls *.fasta 2>/dev/null | wc -l)
    GB_COUNT=\$(ls *.gb 2>/dev/null | wc -l)

    if [[ \$FASTA_COUNT -eq \$GB_COUNT ]]; then
        exit 0
    else
        exit 1
    fi
    """

    stub:
    """
    touch accession.fasta
    touch accession.gb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py -v | sed 's/.* //' )
    END_VERSIONS
    """
}
