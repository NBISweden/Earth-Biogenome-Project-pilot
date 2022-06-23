process BLOBTOOLKIT_VIEW {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "python=3.6 docopt psutil pyyaml ujson tqdm nodejs=10 yq pysam seqtk geckodriver selenium pyvirtualdisplay fastjsonschema" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'genomehubs/blobtoolkit:3.1.3':
        'genomehubs/blobtoolkit:3.1.3' }"

    input:
    tuple val(meta), path(blobdir)
    each plot_type // Plot type (blob|cumulative|snail)

    output:
    tuple val(meta), path("*.png"), optional: true, emit: png
    tuple val(meta), path("*.svg"), optional: true, emit: svg
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blobtools view \\
        --view $plot_type \\
        $args \\
        $blobdir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( blobtools --version | sed 's/.*v//' )
    END_VERSIONS
    """
}
