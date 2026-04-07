process MITOHIFI_FINDMITOREFERENCE {
    tag "$species"
    label 'process_single'
    secret secrets.NCBI_API_KEY ? "NCBI_API_KEY" : ""

    // NOTE: An optional NCBI API key can be supplied to MITOHIFI_FINDMITOREFERENCE.
    // This should be set using Nextflow's secrets functionality:
    // `nextflow secrets set NCBI_API_KEY <key>`
    //
    // See https://www.nextflow.io/docs/latest/secrets.html for more information.

    // Docker image available at the project github repository
    container 'ghcr.io/marcelauliano/mitohifi:3.2.3'

    input:
    tuple val(meta), val(species)

    output:
    tuple val(meta), path("*.fasta"), path("*.gb"), emit: reference, optional: true
    // WARN: Incorrect version information is provided by tool on CLI. Please update this string when bumping container versions.
    // old version command: \$(mitohifi.py -v | sed 's/.* //')
    tuple val("${task.process}"), val('mitohifi'), eval('echo 3.2.3'), emit: versions_mitohifi, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MitoHiFi module does not support Conda. Please use Docker / Singularity instead."
    }

    def args         = task.ext.args ?: ''
    def ncbi_api_key = secrets.NCBI_API_KEY ? "--ncbi-api-key \$NCBI_API_KEY" : ""
    """
    # override workflow exit on findMitoReference.py error (issues: #277, #171, #220)
    set +e
    findMitoReference.py \\
        ${ncbi_api_key} \\
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
    """

    stub:
    """
    touch accession.fasta
    touch accession.gb
    """
}
