include { STAR_ALIGN          } from "$projectDir/modules/nf-core/star/align/main"
include { STAR_GENOMEGENERATE } from "$projectDir/modules/nf-core/star/genomegenerate/main"

workflow ALIGN_RNASEQ {

    take:
    reads_and_assembly_ch // [ meta, primary assembly, [ reads ] ] 

    main:
    reads_and_assembly_ch
        .multiMap { meta, assembly, reads -> 
            assembly_ch: [ meta, assembly ]
            reads_ch: [ meta, reads ]
        }.set { input }
    STAR_GENOMEGENERATE (
        input.assembly_ch,
        []                   // No valid annotations for draft genome
    )
    input.reads_ch.join( STAR_GENOMEGENERATE.out.index )
        .multiMap { meta, reads, index ->
            index_ch: index
            reads_ch: [ meta, reads ]
        }.set { index }
    STAR_ALIGN (
        index.reads_ch,
        index.index_ch,
        [],           // No gtf
        true,         // Ignore gtf
        'Illumina',   // Sequencing Platform
        ''            // Sequencing Center
    )

}