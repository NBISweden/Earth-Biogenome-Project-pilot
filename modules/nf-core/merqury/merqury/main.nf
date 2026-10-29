nextflow.enable.types = true

process MERQURY_MERQURY {
    tag "${id}"
    label 'process_low'

    ext [
        prefix: { output_prefix ?: "${id}" }
    ]

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/merqury:1.3--hdfd78af_4':
        'biocontainers/merqury:1.3--hdfd78af_4' }"

    input:
    record(
        id: String,
        meta: Record?,
        meryl_db: Path,
        assembly: List<Path>,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        assembly_only_kmers_bed: files("*_only.bed"),
        assembly_only_kmers_wig: files("*_only.wig"),
        stats: file("*.completeness.stats"),
        dist_hist: file("*.dist_only.hist"),
        spectra_cn_fl_png: file("*.spectra-cn.fl.png"),
        spectra_cn_hist: file("*.spectra-cn.hist"),
        spectra_cn_ln_png: file("*.spectra-cn.ln.png"),
        spectra_cn_st_png: file("*.spectra-cn.st.png"),
        spectra_asm_fl_png: file("*.spectra-asm.fl.png"),
        spectra_asm_hist: file("*.spectra-asm.hist"),
        spectra_asm_ln_png: file("*.spectra-asm.ln.png"),
        spectra_asm_st_png: file("*.spectra-asm.st.png"),
        assembly_qv: file("${task.ext.prefix}.qv"),
        scaffold_qv: files("${task.ext.prefix}.*.qv"),
        read_ploidy: file("*.hist.ploidy"),
        hapmers_blob_png: file("*.hapmers.blob.png", optional: true)
    )

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'merqury', '1.3') >> 'versions'

    script:
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi
    # limit meryl to use the assigned number of cores.
    export OMP_NUM_THREADS=$task.cpus

    merqury.sh \\
        $meryl_db \\
        ${assembly.join(' ')} \\
        ${task.ext.prefix}
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}_only.bed
    touch ${prefix}_only.wig
    touch ${prefix}.completeness.stats
    touch ${prefix}.dist_only.hist
    touch ${prefix}.spectra-cn.fl.png
    touch ${prefix}.spectra-cn.hist
    touch ${prefix}.spectra-cn.ln.png
    touch ${prefix}.spectra-cn.st.png
    touch ${prefix}.spectra-asm.fl.png
    touch ${prefix}.spectra-asm.hist
    touch ${prefix}.spectra-asm.ln.png
    touch ${prefix}.spectra-asm.st.png
    touch ${prefix}.qv
    ${assembly.collect { asm -> "touch ${prefix}.${asm.baseName}.qv" }.join('\n    ')}
    touch ${prefix}.hist.ploidy
    """
}
