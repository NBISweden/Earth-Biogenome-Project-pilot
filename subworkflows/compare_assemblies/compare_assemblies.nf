include { QUAST } from "$projectDir/modules/nf-core/modules/quast/main"

workflow COMPARE_ASSEMBLIES {

    take:
    assembly_ch        // input type: [ [ id: 'sample_name' ], [ id:'assemblerX_build1', primary_asm_path: '/path/to/primary_asm', alternate_asm_path: '/path/to/alternate_asm' ] ]
    reference_ch       // optional: file( reference_genome ) for comparison

    main:
    QUAST (
        assembly_ch.map { sample, assembly -> assembly.primary_asm_path }
            .collect(),
        reference_ch,
        [], // gff
        reference_ch, // true / false to use reference_ch
        []
    )
    versions_ch = QUAST.out.versions

    emit:
    versions = versions_ch

}