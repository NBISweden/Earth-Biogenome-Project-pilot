// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
process FINDMITOREF {
    tag "$meta.id"
    label 'process_low'

    // There is no conda package available for MitoHiFi
    // Adding the Docker container from biocontainers found on DockerHub
    container 'docker://biocontainers/mitohifi:2.2_cv1'

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/bwa/index/main.nf
    tuple val(meta)

    output:
    tuple val(meta), path("*.fasta"), emit: mitoref
    tuple val(meta), path("*.gb")   , emit: gb

    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    """
    findMitoReference.py \\
        $args \\
        --outfolder /data/ \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        findmitoref: \$( mitohifi.py --version )
    END_VERSIONS
    """
}
