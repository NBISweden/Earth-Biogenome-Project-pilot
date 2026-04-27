process MERQURYFK_KATGC {
    tag "$meta.id"
    label 'process_medium'

    //conda "${moduleDir}/environment.yml"
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    input:
    tuple val(meta), path(fastk_hist), path(fastk_ktab)

    output:
    tuple val(meta), path("*.fi.png"), emit: filled_gc_plot_png , optional: true
    tuple val(meta), path("*.fi.pdf"), emit: filled_gc_plot_pdf , optional: true
    tuple val(meta), path("*.ln.png"), emit: line_gc_plot_png   , optional: true
    tuple val(meta), path("*.ln.pdf"), emit: line_gc_plot_pdf   , optional: true
    tuple val(meta), path("*.st.png"), emit: stacked_gc_plot_png, optional: true
    tuple val(meta), path("*.st.pdf"), emit: stacked_gc_plot_pdf, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('merquryfk'), val('41451f8fb146158c5b747ae7915e69975c61ddd9'), emit: versions_merquryfk, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastk'), val('1.2'), emit: versions_fastk, topic: versions
    tuple val("${task.process}"), val('R'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MERQURYFK_KATGC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    """
    KatGC \\
        $args \\
        -T$task.cpus \\
        ${fastk_ktab.find { path -> path.toString().endsWith(".ktab") }} \\
        $prefix
    """

    stub:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def outfmt          = args.contains('-pdf') ? "pdf" : "png"
    """
    touch ${prefix}.test.fi.${outfmt}
    touch ${prefix}.test.ln.${outfmt}
    touch ${prefix}.test.st.${outfmt}
    """
}
