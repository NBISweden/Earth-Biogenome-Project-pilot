#! /usr/bin/env nextflow

// include { BLOBTOOLKIT     } from "$projectDir/subworkflows/modules/blobtoolkit/blobtoolkit"

include { MERQURYFK_MERQURYFK } from "$projectDir/modules/local/merquryfk/merquryfk"
include { INSPECTOR           } from "$projectDir/modules/local/inspector/inspector"


workflow ASSEMBLY_VALIDATION {

    take:
    assembly_ch        // input type: [ [ id: 'sample_name' ], [ id:'assemblerX_build1', path:'/path/to/assembly' ] ]
    reads_ch           // input type: [ [ id: 'sample_name' ], [ file('path/to/reads') ] ]
    fastk_db           // input type: [ [ id: 'sample_name' ], [ file('path/to/reads.hist') ], [ file('/path/to/reads.ktab') ] ]
    meryl_db           // input type: [ [ id: 'sample_name' ], [ file('path/to/meryl_db') ] ]
    reference_ch       // optional: file( reference_genome ) for comparison
    busco_lineages     // Busco lineages to check against
    busco_lineage_path // Path to Busco lineage files
    diamond_db         // Paths to Diamond databases
    blast_db           // Paths to Blast databases
    ncbi_taxonomy      // Path to ncbi taxonomy database


    /* Assembly validation workflow:
        - Contamination check ( BLOBTOOLKIT )
        - K-mer spectra check
        - Coverage check ( BLOBTOOLKIT )
        - Gene space check ( BUSCO )
        - Mis assembly signal check
    */
    main:
    // BUSCO( assembly )
    QUAST (
        assembly_ch.map { sample, assembly -> assembly.path }
            .collect(),
        reference_ch,
        [], // gff
        reference_ch, // true / false to use reference_ch
        []
    )
    versions_ch = QUAST.out.versions

    // Construct input channel = [ [id: 'name'], [ file(read1), file(read2) ], file(assembly) ]
    id_reads_asm_ch = reads_ch.combine( assembly_ch.map { sample, assembly -> [ sample, assembly.path ] }, by: 0 )
    
    INSPECTOR (
        id_reads_asm_ch,
        reference_ch
    )
    // BLOBTOOLKIT( 
    //     id_reads_asm_ch,
    //     busco_lineages,
    //     busco_lineage_path,
    //     diamond_db,
    //     blast_db,
    //     ncbi_taxonomy 
    // )
    MERQURYFK_MERQURYFK (
        fastk_db.combine( assembly_ch.map { sample, assembly -> [ sample, assembly.path ] }, by: 0 )
    )
    
    versions_ch = versions_ch.mix ( 
        INSPECTOR.out.versions.first(),
        BLOBTOOLKIT.out.versions.first(),
        MERQURY.out.versions.first(),
        MERQURYFK_MERQURYFK.out.versions.first()
    )

    emit:
    versions = versions_ch

}
