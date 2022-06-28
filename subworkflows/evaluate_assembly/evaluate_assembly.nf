include { BUSCO               } from "$projectDir/modules/nf-core/modules/busco/main"
include { MERQURYFK_MERQURYFK } from "$projectDir/modules/nf-core/modules/merquryfk/merquryfk/main"
include { INSPECTOR           } from "$projectDir/modules/local/inspector/inspector"

workflow EVALUATE_ASSEMBLY {

    take:
    assembly_ch        // input type: [ [ id: 'sample_name' ], [ id:'assemblerX_build1', primary_asm_path: '/path/to/primary_asm', alternate_asm_path: '/path/to/alternate_asm' ] ]
    reads_ch           // input type: [ [ id: 'sample_name' ], [ file('path/to/reads') ] ]
    fastk_db           // input type: [ [ id: 'sample_name' ], [ file('path/to/reads.hist') ], [ file('/path/to/reads.ktab') ] ]
    reference_ch       // optional: file( reference_genome ) for comparison
    busco_lineages     // Busco lineages to check against
    busco_lineage_path // Path to Busco lineage files

    main:

    // Kmer consistency check
    MERQURYFK_MERQURYFK (
        fastk_db.combine( assembly_ch.map { sample, assembly ->
            [
                sample, 
                ( assembly.alternate_asm_path ? [ assembly.primary_asm_path, assembly.alternate_asm_path ] : assembly.primary_asm_path ),
                assembly.id
            ] 
        }, by: 0 ).map {
            sample, fastk_hist, fastk_ktab, asm_files, build_name -> 
                [ 
                    [ id: sample.id , build: build_name ],
                    fastk_hist,
                    fastk_ktab,
                    asm_files
                ]
        }
    )

    // Read consistency check
    INSPECTOR (
        reads_ch.combine( assembly_ch.map { sample, assembly -> [ sample, assembly.primary_asm_path, assembly.id ] }, by: 0 )
            .map { sample, reads, asm_files, build_name -> 
                [
                    [ id: sample.id , build: build_name ],
                    reads,
                    asm_files
                ]
            },
        reference_ch
    )

    // Evaluate core gene space coverage
    BUSCO (
        assembly_ch.map { sample, assembly -> [ [ id: sample.id, build: assembly.id ] , assembly.primary_asm_path ] },
        busco_lineages,
        busco_lineage_path,
        []
    )
}