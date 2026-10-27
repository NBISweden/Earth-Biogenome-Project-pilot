nextflow.enable.types = true

process REPORT_GENOMETRAITS {
    tag "${sample_name}"
    label 'process_single'

    input:
    record(
        id: String,
        sample_name: String,
        genome_size: String?,
        haploid_number: String?,
        ploidy: String?
    )

    output:
    record(
        id: id,
        tsv: file("Genome_traits.tsv")
    )

    exec:
    def trait_table = """\
        Genome traits\tExpected\tObserved
        Haploid Size\t${genome_size}\tunknown
        Haploid Number\t${haploid_number}\tunknown
        Ploidy\t${ploidy}\tunknown
        """.stripIndent()
    file("${task.workDir}/Genome_traits.tsv").text = trait_table
}
