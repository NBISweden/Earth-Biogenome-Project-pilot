process BLOBTOOLKIT_CREATE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "python=3.6 docopt psutil pyyaml ujson tqdm nodejs=10 yq pysam seqtk geckodriver selenium pyvirtualdisplay fastjsonschema" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'genomehubs/blobtoolkit:3.1.3':
        'genomehubs/blobtoolkit:3.1.3' }"

    input:
    tuple val(meta), path(assembly)
    path( blobtools_meta )

    output:
    tuple val(meta), path("$prefix"), emit: blobdir
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bmeta = blobtools_meta ? "--meta $blobtools_meta" : ''
    """
    blobtools create \\
        --threads $task.cpus \\
        --fasta $assembly \\
        $bmeta \\
        $args \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( blobtools --version | sed 's/.*v//' )
    END_VERSIONS
    """
}
