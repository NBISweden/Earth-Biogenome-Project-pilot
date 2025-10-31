process MITOHIFI_FINDMITOREFERENCE {
    tag "$species"
    label 'process_low'

    // Docker image available at the project github repository
    container 'ghcr.io/marcelauliano/mitohifi:master'

    input:
    tuple val(meta), val(species)
    val _taxonomy_check_success

    output:
    tuple val(meta), path("*.fasta"), emit: fasta, optional: true
    tuple val(meta), path("*.gb")   , emit: gb, optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # override workflow exit on findMitoReference.py error (issues: #277, #171, #220)
    set +e
    findMitoReference.py \\
        --species "$species" \\
        --outfolder . \\
        $args
    set -e

    # Test for mitohifi reference files:
    # *.fasta && *.gb: Mitohifi complete output, proceed with workflow
    # ! *.fasta && ! *.gb: Mitohifi no output: species not found or other errors, proceed with workflow
    # *.fasta && ! *.gb: Mitohifi partial output: program crashed or core dumped, exit workflow

    FASTA=\$(find . -maxdepth 1 -name "*.fasta" -type f | wc -l)
    GB=\$(find . -maxdepth 1 -name "*.gb" -type f | wc -l)

    if [[ \$FASTA -ne \$GB ]]; then
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py -v | sed 's/.* //' )
    END_VERSIONS
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
