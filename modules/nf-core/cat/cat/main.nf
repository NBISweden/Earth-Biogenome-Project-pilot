nextflow.enable.types = true

process CAT_CAT {
    tag "${id}"
    label 'process_low'

    ext [
        args: { cat_args ?: '' },
        args2: { cat_args2 ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.8' :
        'biocontainers/pigz:2.8' }"

    input:
    record(
        id: String,
        meta: Record?,
        files_in: List<Path>,
        cat_args: String?,
        cat_args2: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        file_out: file("${task.ext.prefix}")
    )

    topic:
    tuple(task.process, "pigz", eval("pigz --version 2>&1 | sed 's/pigz //g'")) >> 'versions'

    script:
    def file_list = files_in.collect { f -> f.name }.toList()

    // choose appropriate concatenation tool depending on input and output format

    // | input     | output     | command1 | command2 |
    // |-----------|------------|----------|----------|
    // | gzipped   | gzipped    | cat      |          |
    // | ungzipped | ungzipped  | cat      |          |
    // | gzipped   | ungzipped  | zcat     |          |
    // | ungzipped | gzipped    | cat      | pigz     |

    def prefix   = task.ext.prefix
    def out_zip  = prefix.endsWith('.gz')
    def in_zip   = file_list[0].endsWith('.gz')
    def command1 = (in_zip && !out_zip) ? 'zcat' : 'cat'
    def command2 = (!in_zip && out_zip) ? "| pigz -c -p $task.cpus ${task.ext.args2}" : ''
    if(file_list.contains(prefix.trim())) {
        error "The name of the input file can't be the same as for the output prefix in the " +
        "module CAT_CAT (currently `$prefix`). Please choose a different one."
    }
    """
    $command1 \\
        ${task.ext.args} \\
        ${file_list.join(' ')} \\
        $command2 \\
        > ${prefix}
    """

    stub:
    def file_list = files_in.collect { f -> f.name }.toList()
    def prefix    = task.ext.prefix
    if(file_list.contains(prefix.trim())) {
        error "The name of the input file can't be the same as for the output prefix in the " +
        "module CAT_CAT (currently `$prefix`). Please choose a different one."
    }
    """
    touch $prefix
    """
}
