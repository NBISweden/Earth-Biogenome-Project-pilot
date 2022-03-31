#! /usr/bin/env nextflow

// include { BUSCO       } from "$projectDir/modules/busco"
include { BLOBTOOLKIT    } from "$projectDir/modules/local/blobtoolkit/blobtoolkit"
include { QUAST          } from "$projectDir/modules/nf-core/modules/quast/main"
include { INSPECTOR      } from "$projectDir/modules/local/inspector/inspector"

workflow ASSEMBLY_VALIDATION {

    take:
    assembly_ch        // input type: [ [ id: 'sample_name' ], [ id:'assemblerX_build1', path:'/path/to/assembly' ] ]
    reads_ch           // input type: [ [ id: 'sample_name' ], [ file('path/to/reads') ] ]
    reference_ch       // optional: file( reference_genome ) for comparison
    busco_lineages     // Busco lineages to check against
    // Is this better as a database map? What are implications - these must be channels!
    busco_lineage_path // Path to Busco lineage files
    uniprot_db         // Path to Uniprot database
    ncbi_nt_db         // Path the ncbi nt database
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

    // Construct input channel = [ [id: 'name'], [ file(read1), file(read2) ], file(assembly) ]
    id_reads_asm_ch = reads_ch.combine( assembly_ch.map { sample, assembly -> [ sample, assembly.path ] }, by: 0 )
    
    INSPECTOR (
        id_reads_asm_ch,
        reference_ch
    )
    BLOBTOOLKIT( 
        id_reads_asm_ch,
        busco_lineages,
        busco_lineage_path,
        uniprot_db,
        ncbi_nt_db,
        ncbi_taxonomy 
    )

}
