nextflow.enable.types = true

process GOAT_TAXONSEARCH {
    tag "${id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/goat:0.2.5--h9d3141d_2':
        'biocontainers/goat:0.2.5--h9d3141d_2' }"

    ext [
        args: { goat_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        taxon: String,
        taxa_file: Path?,
        goat_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        taxonsearch: file("*.tsv")
    )

    topic:
    tuple(task.process, 'goat', eval("goat-cli --version | cut -d' ' -f2")) >> 'versions'

    script:
        def args = task.ext.args
        def prefix = task.ext.prefix
        def input = taxa_file ? "-f ${taxa_file}" : "-t \"${taxon}\""
        if (!taxon && !taxa_file) error "No input. Valid input: single taxon identifier or a .txt file with identifiers"
        if (taxon && taxa_file ) error "Only one input is required: a single taxon identifier or a .txt file with identifiers"
        """
        goat-cli taxon search \\
            $args \\
            $input > ${prefix}.tsv
        """

    stub:
        def prefix = task.ext.prefix
        if (!taxon && !taxa_file) error "No input. Valid input: single taxon identifier or a .txt file with identifiers"
        if (taxon && taxa_file ) error "Only one input is required: a single taxon identifier or a .txt file with identifiers"
        """
        touch ${prefix}.tsv
        """
}
