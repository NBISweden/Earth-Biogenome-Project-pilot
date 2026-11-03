nextflow.enable.types = true

process MERQURYFK_KATGC {
    tag "${id}"
    label 'process_medium'

    ext [
        args: { katgc_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    //conda "${moduleDir}/environment.yml"
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    input:
    record(
        id: String,
        meta: Record?,
        fastk_hist: Path,
        fastk_ktab: List<Path>,
        katgc_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        filled_gc_plot_png: file("*.fi.png", optional: true),
        filled_gc_plot_pdf: file("*.fi.pdf", optional: true),
        line_gc_plot_png: file("*.ln.png", optional: true),
        line_gc_plot_pdf: file("*.ln.pdf", optional: true),
        stacked_gc_plot_png: file("*.st.png", optional: true),
        stacked_gc_plot_pdf: file("*.st.pdf", optional: true)
    )

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'merquryfk', '41451f8fb146158c5b747ae7915e69975c61ddd9') >> 'versions'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'fastk', '1.2') >> 'versions'
    tuple(task.process, 'R', eval('R --version | sed "1!d; s/.*version //; s/ .*//"')) >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MERQURYFK_KATGC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    KatGC \\
        ${task.ext.args} \\
        -T$task.cpus \\
        ${fastk_ktab.find { path -> path.name.endsWith(".ktab") }} \\
        ${task.ext.prefix}
    """

    stub:
    def prefix = task.ext.prefix
    def outfmt = task.ext.args.contains('-pdf') ? "pdf" : "png"
    """
    touch ${prefix}.test.fi.${outfmt}
    touch ${prefix}.test.ln.${outfmt}
    touch ${prefix}.test.st.${outfmt}
    """
}
