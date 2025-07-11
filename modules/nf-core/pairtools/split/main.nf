process PAIRTOOLS_SPLIT {
    tag "${meta.id}"
    label 'process_medium'

    // Pinning numpy to 1.23 until https://github.com/open2c/pairtools/issues/170 is resolved
    // Not an issue with the biocontainers because they were built prior to numpy 1.24
    conda "bioconda::pairtools=1.1.0 conda-forge::numpy=1.23"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.1.0--py39hd5a99d8_2' :
        'biocontainers/pairtools:1.1.0--py39hd5a99d8_2' }"

    input:
    tuple val(meta), path(pairs)
    tuple val(meta2), path(fasta)
    val   sort_bam

    output:
    tuple val(meta), path("*.split.pairs.gz"), emit:pairs
    tuple val(meta), path("*.bam")           , emit: bam , optional:true
    tuple val(meta), path("*.cram")          , emit: cram, optional:true
    tuple val(meta), path("*.crai")          , emit: crai, optional:true
    tuple val(meta), path("*.csi")           , emit: csi , optional:true
    path("versions.yml")                     , emit:versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'

    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"

    """
    export MPLCONFIGDIR=tmp
    pairtools split \
        --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
        --output-pairs ${prefix}.split.pairs.gz \
        --output-sam - \
        ${args} \
        ${pairs} \
        | samtools $samtools_command $args2 -@ $task.cpus ${reference} -o ${prefix}.split.pairs.${extension} -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools, version //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=tmp
    echo "" | gzip > ${prefix}.split.pairs.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version | sed 's/pairtools, version //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}
