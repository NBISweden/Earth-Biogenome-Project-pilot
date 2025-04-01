import groovy.transform.RecordOptions

@RecordOptions(components = true, copyWith = true)
record SampleInfo (
    String sampleName,
    Integer taxId,
    Integer geneticCode,
    Integer mitoCode,
    String domain,
    Integer genomeSize,
    Integer ploidy,
    Integer haploidNumber
) {
    SampleInfo { assert !sampleName.blank }

    boolean hasComponents(List<String> components) {
        return components.every { this[it] }
    }

    boolean isEukaryota() {
        return domain == 'Eukaryota'
    }
}
