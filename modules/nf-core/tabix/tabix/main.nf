nextflow.enable.types = true

process TABIX_TABIX {
    tag "${id}"
    label 'process_single'

    ext [
        args: { tabix_args ?: '' }
    ]

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data' :
        'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa' }"

    input:
    record(
        id: String,
        meta: Record?,
        tab: Path,
        tabix_args: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        index: file("*.{tbi,csi}")
    )

    topic:
    tuple(task.process, 'tabix', eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'")) >> 'versions'

    script:
    """
    tabix \\
        --threads $task.cpus \\
        ${task.ext.args} \\
        $tab

    """

    stub:
    def index = task.ext.args.contains("-C ") || task.ext.args.contains("--csi") ? "csi" : "tbi"
    """
    touch ${tab}.${index}
    """
}
