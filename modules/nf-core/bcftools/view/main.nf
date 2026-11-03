nextflow.enable.types = true

process BCFTOOLS_VIEW {
    tag "${id}"
    label 'process_medium'

    ext [
        args: { bcftools_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    record(
        id: String,
        meta: Record?,
        vcf: Path,
        index: Path,
        regions: Path?,
        targets: Path?,
        samples: Path?,
        bcftools_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        vcf: file("*.{vcf,vcf.gz,bcf,bcf.gz}"),
        tbi: file("*.tbi", optional: true),
        csi: file("*.csi", optional: true)
    )

    topic:
    tuple(task.process, 'bcftools', eval("bcftools --version | sed '1!d; s/^.*bcftools //'")) >> 'versions'

    script:
    def args = task.ext.args
    def prefix = task.ext.prefix
    def regions_file = regions ? "--regions-file ${regions}" : ""
    def targets_file = targets ? "--targets-file ${targets}" : ""
    def samples_file = samples ? "--samples-file ${samples}" : ""
    def extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    """
    bcftools view \\
        --output ${prefix}.${extension} \\
        ${regions_file} \\
        ${targets_file} \\
        ${samples_file} \\
        ${args} \\
        --threads ${task.cpus} \\
        ${vcf}
    """

    stub:
    def args = task.ext.args
    def prefix = task.ext.prefix
    def extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def stub_index = args.contains("--write-index=tbi") || args.contains("-W=tbi")
        ? "tbi"
        : args.contains("--write-index=csi") || args.contains("-W=csi")
            ? "csi"
            : args.contains("--write-index") || args.contains("-W")
                ? "csi"
                : ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && ["csi", "tbi"].contains(stub_index) ? "touch ${prefix}.${extension}.${stub_index}" : ""
    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}
    """
}
