include { STAR_ALIGN          } from "$projectDir/modules/nf-core/star/align/main"
include { STAR_GENOMEGENERATE } from "$projectDir/modules/nf-core/star/genomegenerate/main"

workflow ALIGN_RNASEQ {

    take:
    reads_ch
    assembly_ch

    main:
    STAR_GENOMEGENERATE (
        assembly_ch,
        [] // No valid annotations for draft genome
    )
    STAR_GENOMEGENERATE.out.index
        .map { meta, index -> def keylist = ['build']; [ meta.findAll { !(it.key in keylist) }, meta.subMap(keylist), index ] }
        .combine( reads_ch.map { meta, reads -> [ meta.findAll { !(it.key in ['single_end'])}, reads ] }, by: 0 )
        .multiMap { sample_meta, build_meta, index, reads ->
            reads_ch: [ sample_meta + build_meta , reads ]
            index_ch: index
        }.set { indexed }
    STAR_ALIGN (
        indexed.reads_ch,
        indexed.index_ch,
        [],           // No gtf
        true,         // Ignore gtf
        'Illumina',   // Sequencing Platform
        ''            // Sequencing Center
    )

}