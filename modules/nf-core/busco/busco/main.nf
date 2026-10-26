nextflow.enable.types = true

process BUSCO_BUSCO {
    tag "${id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/6456c1880785adefb4fc9b480bb7662479d5662c17f70d5e8715b7f2a63ee28b/data'
        : 'community.wave.seqera.io/library/busco_numpy:b66937518a305dd7'}"
    // Note: one test had to be disabled when switching to Busco 6.0.0, cf https://github.com/nf-core/modules/pull/8781/files
    // Try to restore it when upgrading Busco to a later version

    ext [
        args: { busco_args ?: '' },
        prefix: { output_prefix ?: "${id}-${lineage}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        fasta: Path,
        mode: String,
        lineage: String,
        busco_lineages_path: Path?,
        config_file: Path?,
        clean_intermediates: Boolean,
        busco_args: String?,
        output_prefix: String?
    )

    stage:
    stageAs fasta, 'tmp_input/*'

    output:
    record(
        id: id,
        meta: meta,
        batch_summary: file("*-busco.batch_summary.txt"),
        batch_summary_failed: file("*-busco.batch_summary.failed.txt", optional: true),
        short_summaries_txt: file("short_summary.*.txt", optional: true),
        short_summaries_json: file("short_summary.*.json", optional: true),
        log: file("*-busco.log", optional: true),
        full_table: file("*-busco/*/run_*/full_table.tsv", optional: true),
        missing_busco_list: file("*-busco/*/run_*/missing_busco_list.tsv", optional: true),
        single_copy_proteins: file("*-busco/*/run_*/single_copy_proteins.faa", optional: true),
        seq_dir: file("*-busco/*/run_*/busco_sequences", optional: true),
        translated_dir: file("*-busco/*/translated_proteins", optional: true),
        busco_dir: file("*-busco"),
        downloaded_lineages: files("busco_downloads/lineages/*", optional: true),
        single_copy_faa: files("*-busco/*/run_*/busco_sequences/single_copy_busco_sequences/*.faa", optional: true),
        single_copy_fna: files("*-busco/*/run_*/busco_sequences/single_copy_busco_sequences/*.fna", optional: true)
    )

    topic:
    tuple(task.process, 'busco', eval("busco --version 2> /dev/null | sed 's/BUSCO //g'")) >> 'versions'

    script:
    if (mode !in ['genome', 'proteins', 'transcriptome']) {
        error("Mode must be one of 'genome', 'proteins', or 'transcriptome'.")
    }
    def args = task.ext.args
    def prefix = task.ext.prefix
    def busco_config = config_file ? "--config ${config_file}" : ''
    def busco_lineage = lineage in ['auto', 'auto_prok', 'auto_euk']
        ? lineage.replaceFirst('auto', '--auto-lineage').replaceAll('_', '-')
        : "--lineage_dataset ${lineage}"
    def busco_lineage_dir = busco_lineages_path ? "--download_path ${busco_lineages_path}" : ''
    def intermediate_files = [
        './*-busco/*/auto_lineage',
        './*-busco/*/**/{miniprot,hmmer,.bbtools}_output',
        './*-busco/*/prodigal_output/predicted_genes/tmp/',
    ]
    def clean_cmd = clean_intermediates ? "rm -fr ${intermediate_files.join(' ')}" : ''
    """
    # Fix Augustus for Apptainer
    ENV_AUGUSTUS=/opt/conda/etc/conda/activate.d/augustus.sh
    set +u
    if [ -z "\${AUGUSTUS_CONFIG_PATH}" ] && [ -f "\${ENV_AUGUSTUS}" ]; then
        source "\${ENV_AUGUSTUS}"
    fi
    set -u

    # If the augustus config directory is not writable, then copy to writeable area
    if [ ! -w "\${AUGUSTUS_CONFIG_PATH}" ]; then
        # Create writable tmp directory for augustus
        AUG_CONF_DIR=\$( mktemp -d -p \$PWD )
        cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
        export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
        echo "New AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_PATH}"
    fi

    # Ensure the input is uncompressed
    INPUT_SEQS=input_seqs
    mkdir "\$INPUT_SEQS"
    cd "\$INPUT_SEQS"
    for FASTA in ../tmp_input/*; do
        if [ "\${FASTA##*.}" == 'gz' ]; then
            gzip -cdf "\$FASTA" > \$( basename "\$FASTA" .gz )
        else
            ln -s "\$FASTA" .
        fi
    done
    cd ..

    busco \\
        --cpu ${task.cpus} \\
        --in "\$INPUT_SEQS" \\
        --out ${prefix}-busco \\
        --mode ${mode} \\
        ${busco_lineage} \\
        ${busco_lineage_dir} \\
        ${busco_config} \\
        ${args}

    # clean up
    rm -rf "\$INPUT_SEQS"
    ${clean_cmd}
    # find and remove broken symlinks from the cleanup
    find . -xtype l -delete

    # Move files to avoid staging/publishing issues
    grep -v 'Run failed; check logs' ${prefix}-busco/batch_summary.txt > ${prefix}-busco.batch_summary.txt
    mv ${prefix}-busco/*/short_summary.*.{json,txt} . || echo "Short summaries were not available: No genes were found."
    mv ${prefix}-busco/logs/busco.log ${prefix}-busco.log

    if grep -q 'Run failed; check logs' ${prefix}-busco/batch_summary.txt; then
        grep 'Run failed; check logs' ${prefix}-busco/batch_summary.txt > ${prefix}-busco.batch_summary.failed.txt
    fi
    """

    stub:
    def prefix = task.ext.prefix
    def fasta_name = fasta.name.replaceAll(/\.gz$/, '')
    """
    touch ${prefix}-busco.batch_summary.txt
    mkdir -p ${prefix}-busco/${fasta_name}/run_${lineage}/busco_sequences
    """
}
