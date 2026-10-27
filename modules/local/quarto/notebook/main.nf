nextflow.enable.types = true

process QUARTO_NOTEBOOK {
    tag "${notebook.baseName}"
    label 'process_single'

    ext [
        args: { quarto_args ?: '' },
        prefix: { output_prefix ?: notebook.baseName }
    ]

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ca/ca05bf35df4f90442c46f9a54f8b99aebe136608d7c78dd56fde39520b3bd620/data':
        'community.wave.seqera.io/library/multiqc_jupyter_pandas_papermill_pruned:c3a45031dd77805d' }"

    input:
    record(
        id: String,
        meta: Record?,
        notebook: Path,
        aux_files: List<Path>,
        log_files: List<Path>,
        report_params: Path,
        quarto_args: String?,
        output_prefix: String?
    )

    stage:
    stageAs log_files, 'log_files/*'
    stageAs report_params, 'params.yml'

    output:
    record(
        id: id,
        meta: meta,
        html: file("${task.ext.prefix}.html"),
        github_markdown: file("${task.ext.prefix}.md"),
        multiqc_summary: file("multiqc*.html"),
        versions: file("versions.yml")
    )

    script:
    def args = task.ext.args
    """
    export USERID=\$UID
    export XDG_CACHE_HOME=tmp/quarto_cache_home
    export XDG_DATA_HOME=tmp/quarto_data_home
    # Fix Quarto for apptainer
    ENV_QUARTO="\${ENV_QUARTO:-/opt/conda/etc/conda/activate.d/quarto.sh}"
    set +u
    if [ -z "\${QUARTO_DENO}" ] && [ -f "\${ENV_QUARTO}" ]; then
        source "\${ENV_QUARTO}"
    fi
    set -u

    # Link params to meta data
    ln -s params.yml _quarto.yml

    quarto check
    quarto \\
        render \\
        $notebook \\
        --execute-params params.yml \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
      quarto: \$(quarto --version)
      multiqc: \$(multiqc --version | sed '1!d; s/.*version //')
    \$(cat log_files/assembly_mqc_versions.yml)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}.html ${prefix}.md multiqc_stub.html versions.yml
    """
}
