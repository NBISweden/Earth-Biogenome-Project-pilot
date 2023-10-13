process PURGEDUPS_SPLITFA {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::purge_dups=1.2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--h7132678_0':
        'biocontainers/purge_dups:1.2.6--h7132678_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.split.fasta.gz") , emit: split_fasta
    tuple val(meta), path("*.merged.fasta.gz"), emit: merged_fasta
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def useGzip = !( assembly instanceof List ? assembly.every{ it.name.endsWith(".gz") } : assembly.name.endsWith(".gz") )
    """
    ${useGzip ? 'gzip -c' : 'cat'} $assembly > ${prefix}.merged.fasta.gz
    split_fa $args ${prefix}.merged.fasta.gz | gzip -c > ${prefix}.split.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: \$( purge_dups -h |& sed '3!d; s/.*: //' )
    END_VERSIONS
    """
}
