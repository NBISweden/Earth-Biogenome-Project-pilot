include { QUAST } from "$projectDir/modules/nf-core/quast/main"

workflow COMPARE_ASSEMBLIES {

    take:
    assembly_ch        // input type: [ [ id: 'sample_name' ], [ id:'assemblerX_build1', primary_asm_path: '/path/to/primary_asm', alternate_asm_path: '/path/to/alternate_asm' ] ]
    reference_ch       // optional: file( reference_genome ) for comparison

    main:
    QUAST (
        assembly_ch.map { sample, assembly -> assembly.pri_fasta }
            .collect(),
        reference_ch,
        []              // No GFF
    )
    versions_ch = QUAST.out.versions

    emit:
    versions = versions_ch

}