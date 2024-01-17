process FCSGX_FETCHDB {
    tag "$manifest"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.4.0/fcs-gx.sif':
        'docker.io/ncbi/fcs-gx:0.4.0' }"

    input:
    path manifest

    output:
    path "$prefix"      , emit: database
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FCS_FCSGX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: 'gxdb'
    """
    sync_files \\
        get \\
        --mft "$manifest" \\
        --dir "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx -h | sed '/^build.*/ !d; s/.*git:v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: 'gxdb'
    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx -h | sed '/^build.*/ !d; s/.*git:v//' )
    END_VERSIONS
    """
}
