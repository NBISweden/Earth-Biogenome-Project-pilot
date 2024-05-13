// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process DVPOLISH_CHUNKFA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fai)

    output:
    tuple val(meta), path ("*.bed", arity: '1..*')        , emit: bed
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    oneKilobyte = MemoryUnit.of(1000)
    val chunk_size = getBytes(args.find { it.key == 'chunk_size' }?.value ?: '100.MB')
    """
    # convert chunk size into base pairs
    awk -v chunk_size_inBases=${chunk_size} -v prefix=$prefix 'BEGIN {
        block=1
        cum_basecount=0
    }{
        output_file = sprintf("%s_chunk_%d.bed", prefix, block)
        printf("%s\\t0\\t%s\\n", \$1, \$2) > output_file
        cum_basecount+=\$2

        if (cum_basecount >= chunk_size_inBases)
        {
            cum_basecount=0
            block+=1
        }
    }' $fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dvpolish: \$(awk --version |& sed '1!d')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_chunk_1.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dvpolish: \$(awk --version |& sed '1!d')
    END_VERSIONS
    """
}
