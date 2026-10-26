nextflow.enable.types = true

process PRETEXTMAP {
    tag "${id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6b88d972aaf27c7e748e2c3b5ee86065dc5ed7824d7d13937c65844242211e2/data':
        'community.wave.seqera.io/library/htslib_pretextmap_samtools:8a29f6d0f55f98f9' }"

    ext [
        args: { pretextmap_args ?: '' },
        args2: { pretextmap_args2 ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        input: Path,
        fasta: Path?,
        fai: Path?,
        pretextmap_args: String?,
        pretextmap_args2: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        pretext: file("*.pretext")
    )

    topic:
    tuple(task.process, 'PretextMap', eval('PretextMap | sed "/Version/!d; s/.*Version //"')) >> 'versions'
    tuple(task.process, 'samtools', eval('samtools --version | sed "1!d; s/samtools //"')) >> 'versions'

    script:
    def args          = task.ext.args     ?: '' // PretextMap args
    def args2         = task.ext.args2    ?: '' // Samtools view args
    def prefix        = task.ext.prefix   ?: "${id}"
    def reference     = fasta             ? "--reference ${fasta}" : ""
    def pairs_input   = input.name.endsWith(".pairs.gz")
    def input_command = pairs_input ? "zcat ${input}" : "samtools view $args2 $reference -h ${input}"
    """
    ${input_command} | PretextMap \\
        $args \\
        -o ${prefix}.pretext
    """

    stub:
    def prefix      = task.ext.prefix ?: "${id}"
    """
    touch ${prefix}.pretext
    """
}
