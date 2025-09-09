process QUARTO_NOTEBOOK {
    tag "$notebook.baseName"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ca/ca05bf35df4f90442c46f9a54f8b99aebe136608d7c78dd56fde39520b3bd620/data':
        'community.wave.seqera.io/library/multiqc_jupyter_pandas_papermill_pruned:c3a45031dd77805d' }"

    input:
    tuple val(meta), path(notebook, arity: '1'), path(_aux_files)
    path (log_files, stageAs: 'log_files/*', arity: '1..*')
    path 'params.yml'

    output:
    path "${prefix}.html", arity: '1', emit: html
    path "${prefix}.md"  , arity: '1', emit: github_markdown
    path "multiqc*.html" , arity: '1', emit: multiqc_summary
    path "versions.yml"  , arity: '1', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: notebook.baseName
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
    "${task.process}":
        quarto: \$(quarto --version)
        multiqc: \$(multiqc --version | sed '1!d; s/.*version //')
    \$(cat log_files/assembly_mqc_versions.yml)
    END_VERSIONS
    """
}