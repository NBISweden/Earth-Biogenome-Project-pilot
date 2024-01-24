process TAXONKIT_NAME2LINEAGE {
    tag "$meta.sample.name"
    label 'process_low'

    conda "bioconda::taxonkit=0.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/taxonkit:0.15.1--h9ee0642_0':
        'biocontainers/taxonkit:0.15.1--h9ee0642_0' }"

    input:
    val meta
    path taxdb

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: '-f "{k}"'
    def prefix = task.ext.prefix ?: "${meta.sample.name.replace(" ","_")}"
    """
    taxonkit \\
        name2taxid \\
        $args \\
        --data-dir $taxdb \\
        --threads $task.cpus \\
        <<< '${meta.sample.name}' \\
    | taxonkit \\
        reformat \\
        --taxid-field 2 \\
        $args2 \\
        --data-dir $taxdb
        --threads $task.cpus \\
        --out-file $prefix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$( taxonkit version | sed 's/.* v//' )
    END_VERSIONS
    """
}
