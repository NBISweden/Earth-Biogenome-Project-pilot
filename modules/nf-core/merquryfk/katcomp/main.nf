nextflow.enable.types = true

process MERQURYFK_KATCOMP {
    tag "${id}"
    label 'process_medium'

    ext [
        args: { katcomp_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    //conda "${moduleDir}/environment.yml"
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    input:
    record(
        id: String,
        meta: Record?,
        fastk1_hist: Path,
        fastk1_ktab: List<Path>,
        fastk2_hist: Path,
        fastk2_ktab: List<Path>,
        katcomp_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        filled_png: file("*.fi.png", optional: true),
        line_png: file("*.ln.png", optional: true),
        stacked_png: file("*.st.png", optional: true),
        filled_pdf: file("*.fi.pdf", optional: true),
        line_pdf: file("*.ln.pdf", optional: true),
        stacked_pdf: file("*.st.pdf", optional: true)
    )

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'merquryfk', '41451f8fb146158c5b747ae7915e69975c61ddd9') >> 'versions'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'fastk', '1.2') >> 'versions'
    tuple(task.process, 'R', eval('R --version | sed "1!d; s/.*version //; s/ .*//"')) >> 'versions'

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MERQURYFK_KATCOMP module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def input_fk1 = fastk1_ktab.find { path -> path.name.endsWith(".ktab") }.baseName
    def input_fk2 = fastk2_ktab.find { path -> path.name.endsWith(".ktab") }.baseName
    """
    KatComp \\
        ${task.ext.args} \\
        -T$task.cpus \\
        ${input_fk1} \\
        ${input_fk2} \\
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
