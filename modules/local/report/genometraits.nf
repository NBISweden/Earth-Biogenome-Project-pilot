process REPORT_GENOMETRAITS {
    tag "${meta.sample.sampleName()}"
    label 'process_single'

    input:
    val meta

    output:
    path "Genome_traits.tsv", emit: tsv

    when:
    task.ext.when == null || task.ext.when

    exec:
    def trait_table = """\
        Genome traits\tExpected\tObserved
        Haploid Size\t${meta.sample.genomeSize()}\tunknown
        Haploid Number\t${meta.sample.haploidNumber()}\tunknown
        Ploidy\t${meta.sample.ploidy()}\tunknown
        """.stripIndent()
    file("$task.workDir/Genome_traits.tsv").text = trait_table
}
