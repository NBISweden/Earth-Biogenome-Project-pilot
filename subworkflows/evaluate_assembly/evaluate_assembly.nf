include { BUSCO               } from "$projectDir/modules/nf-core/modules/busco/main"
include { MERQURYFK_MERQURYFK } from "$projectDir/modules/nf-core/modules/merquryfk/merquryfk/main"
include { INSPECTOR           } from "$projectDir/modules/local/inspector/inspector"

workflow EVALUATE_ASSEMBLY {

    take:
    assembly_ch        // input type: [ [ id: 'sample_name' ], [ id:'assemblerX_build1', path:'/path/to/assembly' ] ]
    reads_ch           // input type: [ [ id: 'sample_name' ], [ file('path/to/reads') ] ]
    fastk_db           // input type: [ [ id: 'sample_name' ], [ file('path/to/reads.hist') ], [ file('/path/to/reads.ktab') ] ]
    reference_ch       // optional: file( reference_genome ) for comparison
    busco_lineages     // Busco lineages to check against
    busco_lineage_path // Path to Busco lineage files

    main:

    // Kmer consistency check
    MERQURYFK_MERQURYFK (
        fastk_db.combine( assembly_ch.map { sample, assembly -> [ sample, assembly.path ] }, by: 0 )
    )

    // Read consistency check
    INSPECTOR (
        reads_ch.combine( assembly_ch.map { sample, assembly -> [ sample, assembly.path ] }, by: 0 ),
        reference_ch
    )

    // Evaluate core gene space coverage
    BUSCO (
        assembly_ch.map { sample, assembly -> [ sample, assembly.path ] },
        busco_lineages,
        busco_lineage_path,
        []
    )
}