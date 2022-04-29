// nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//          list (`[]`) instead of a file can be used to work around this issue.

process BLOBTOOLKIT_ADD {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "python=3.6 docopt psutil pyyaml ujson tqdm nodejs=10 yq pysam seqtk geckodriver selenium pyvirtualdisplay fastjsonschema" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'genomehubs/blobtoolkit:3.1.3':
        'genomehubs/blobtoolkit:3.1.3' }"

    input:
    tuple val(meta), path(blobdir)
    path (hits)
    path (coverage)
    path (busco)
    path (bed)
    path (beddir)
    path (bedtsv)
    path (taxdump)

    output:
    tuple val(meta), path("$blobdir"), emit: blobdb
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    tax_dump = taxdump ? "--taxdump ${taxdump}" : ''
    bed_dir = beddir ? "--beddir ${beddir}" : ''  
    """
    blobtools add \\
        --threads $task.cpus \\
        ${hits.collect{ "--hits $it" }.join("\\\\\n") }
        ${coverage.collect{ "--cov $it "}.join("\\\\\n") }
        ${busco.collect{ "--busco $it" }.join("\\\\\n") }
        ${bed.collect{ "--busco $it" }.join("\\\\\n") }
        ${bedtsv.collect{ "--bedtsv $it" }.join("\\\\\n") }
        $tax_dump \\
        $bed_dir \\
        $args \\
        $blobdir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( blobtools --version | sed 's/.*v//' )
    END_VERSIONS
    """
}
