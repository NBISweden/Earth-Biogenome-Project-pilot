nextflow.enable.types = true

process DEEPVARIANT_RUNDEEPVARIANT {
    tag "${id}"
    label 'process_high'

    ext [
        args: { deepvariant_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    // FIXME Conda is not supported at the moment
    // https://github.com/bioconda/bioconda-recipes/pull/45214#issuecomment-1890937836
    // BUG https://github.com/nf-core/modules/issues/1754
    // BUG https://github.com/bioconda/bioconda-recipes/issues/30310
    container "docker.io/google/deepvariant:1.10.0"

    input:
    record(
        id: String,
        meta: Record?,
        input: Path,
        index: Path,
        intervals: Path?,
        fasta: Path,
        fai: Path,
        gzi: Path?,
        par_bed: Path?,
        deepvariant_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        vcf: file("${task.ext.prefix}.vcf.gz"),
        vcf_tbi: file("${task.ext.prefix}.vcf.gz.tbi"),
        gvcf: file("${task.ext.prefix}.g.vcf.gz"),
        gvcf_tbi: file("${task.ext.prefix}.g.vcf.gz.tbi"),
        report: file("${task.ext.prefix}.visual_report.html", optional: true)
    )

    topic:
    tuple(task.process, 'deepvariant', eval("/opt/deepvariant/bin/run_deepvariant --version | sed 's/^.*version //'")) >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix
    def regions = intervals ? "--regions=${intervals}" : ""
    def par_regions = par_bed ? "--par_regions_bed=${par_bed}" : ""

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --ref=${fasta} \\
        --reads=${input} \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        ${task.ext.args} \\
        ${regions} \\
        ${par_regions} \\
        --intermediate_results_dir=tmp \\
        --num_shards=${task.cpus}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix
    """
    echo "stub" | gzip > ${prefix}.vcf.gz
    echo "stub" | gzip > ${prefix}.vcf.gz.tbi
    echo "stub" | gzip > ${prefix}.g.vcf.gz
    echo "stub" | gzip > ${prefix}.g.vcf.gz.tbi
    """
}
