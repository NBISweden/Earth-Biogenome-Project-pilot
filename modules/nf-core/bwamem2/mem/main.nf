nextflow.enable.types = true

process BWAMEM2_MEM {
    tag "${id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05ce34b46ad42810eb29f74e4e304c0cb592b2ca15572929ed8bbaee58faf01/data'
        : 'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113'}"

    ext [
        args: { bwamem2_args ?: '' },
        args2: { bwamem2_args2 ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        reads: List<Path>,
        index: Path,
        fasta: Path?,
        sort_bam: Boolean,
        bwamem2_args: String?,
        bwamem2_args2: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        sam: file("*.sam", optional: true),
        bam: file("*.bam", optional: true),
        cram: file("*.cram", optional: true),
        crai: file("*.crai", optional: true),
        csi: file("*.csi", optional: true)
    )

    topic:
    tuple(task.process, 'bwamem2', eval('bwa-mem2 version | grep -o -E "[0-9]+(\\.[0-9]+)+"')) >> 'versions'
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    def samtools_command = sort_bam ? 'sort' : 'view'

    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher = (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension == "cram" ? "--reference ${fasta}" : ""
    if (!fasta && extension == "cram") {
        error("Fasta reference is required for CRAM output")
    }
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa-mem2 \\
        mem \\
        ${args} \\
        -t ${task.cpus} \\
        \$INDEX \\
        ${reads.join(' ')} \\
        | samtools ${samtools_command} ${args2} -@ ${task.cpus} ${reference} -o ${prefix}.${extension} -
    """

    stub:
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher = (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    if (!fasta && extension == "cram") {
        error("Fasta reference is required for CRAM output")
    }

    def create_index = ""
    if (extension == "cram") {
        create_index = "touch ${prefix}.crai"
    }
    else if (extension == "bam") {
        create_index = "touch ${prefix}.csi"
    }
    """
    touch ${prefix}.${extension}
    ${create_index}
    """
}
