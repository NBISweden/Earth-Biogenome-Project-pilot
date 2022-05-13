#! /usr/bin/env nextflow

include { BLOBTOOLKIT     } from "$projectDir/subworkflows/modules/blobtoolkit/blobtoolkit"
include { QUAST           } from "$projectDir/modules/nf-core/modules/quast/main"
include { MERYL_COUNT     } from "$projectDir/modules/nf-core/modules/meryl/count/main"
include { MERYL_UNIONSUM  } from "$projectDir/modules/nf-core/modules/meryl/unionsum/main"
include { MERYL_HISTOGRAM } from "$projectDir/modules/nf-core/modules/meryl/histogram/main"
include { GENOMESCOPE2    } from "$projectDir/modules/nf-core/modules/genomescope2/main"
include { MERQURY         } from "$projectDir/modules/local/merqury"
include { INSPECTOR       } from "$projectDir/modules/local/inspector/inspector"


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
    MERYL_COUNT ( reads_ch )
    MERYL_UNIONSUM ( MERYL_COUNT.out.meryl_db )
    MERYL_HISTOGRAM ( MERYL_UNIONSUM.out.meryl_db )
    GENOMESCOPE2 ( MERYL_HISTOGRAM.out.hist )
    MERQURY (
        MERYL_UNIONSUM.out.meryl_db
            .combine( assembly_ch.map { sample, assembly -> [ sample, assembly.path ] }, by: 0 )
    )

}
