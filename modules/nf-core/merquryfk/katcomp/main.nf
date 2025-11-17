process MERQURYFK_KATCOMP {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.3'

    input:
    tuple val(meta), path(fastk1_hist), path(fastk1_ktab), path(fastk2_hist), path(fastk2_ktab)

    output:
    tuple val(meta), path("*.fi.png"), emit: filled_png , optional: true
    tuple val(meta), path("*.ln.png"), emit: line_png   , optional: true
    tuple val(meta), path("*.st.png"), emit: stacked_png, optional: true
    tuple val(meta), path("*.fi.pdf"), emit: filled_pdf , optional: true
    tuple val(meta), path("*.ln.pdf"), emit: line_pdf   , optional: true
    tuple val(meta), path("*.st.pdf"), emit: stacked_pdf, optional: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MERQURYFK_KATCOMP module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '0e24fb45b71c4e14382ae1e1bc063bf66ea4e112' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '41451f8fb146158c5b747ae7915e69975c61ddd9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    KatComp \\
        $args \\
        -T$task.cpus \\
        ${fastk1_ktab.find{ it.toString().endsWith(".ktab") }} \\
        ${fastk2_ktab.find{ it.toString().endsWith(".ktab") }} \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """
}
