include { BUSCO               } from "$projectDir/modules/nf-core/modules/busco/main"
include { MERQURYFK_MERQURYFK } from "$projectDir/modules/nf-core/modules/merquryfk/merquryfk/main"
include { INSPECTOR           } from "$projectDir/modules/local/inspector/inspector"

workflow EVALUATE_ASSEMBLY {

    take:
    assembly_ch        // input type: [ meta, [ id:'assemblerX_build1', pri_fasta: '/path/to/primary_asm', alt_fasta: '/path/to/alternate_asm' ] ]
    reads_ch           // input type: [ meta, [ 'path/to/reads' ] ]
    fastk_db           // input type: [ meta, [ 'path/to/reads.hist' ], [ '/path/to/reads.ktab' ] ]
    reference_ch       // optional: file( reference_genome ) for comparison
    busco_lineages     // Busco lineages to check against
    busco_lineage_path // Path to Busco lineage files

    main:

    // Kmer consistency check
    MERQURYFK_MERQURYFK (
        fastk_db.combine( assembly_ch.map { sample, assembly ->
            [
                sample, 
                ( assembly.alt_fasta ? [ assembly.pri_fasta, assembly.alt_fasta ] : assembly.pri_fasta ),
                assembly.id
            ] 
        }, by: 0 ).map {
            sample, fastk_hist, fastk_ktab, asm_files, build_name -> 
                [ 
                    sample + [ build: build_name ],
                    fastk_hist,
                    fastk_ktab,
                    asm_files
                ]
        }
    )

    // Read consistency check
    INSPECTOR (
        reads_ch.map { meta, files -> [ meta.findAll { it.key != 'single_end' }, files ] } // strip read pairing info from meta map
            .combine( assembly_ch.map { sample, assembly -> [ sample, assembly.pri_fasta, assembly.id ] }, by: 0 )
            .map { sample, reads, asm_files, build_name -> 
                [
                    sample + [ build: build_name ],
                    reads,
                    asm_files
                ]
            },
        reference_ch
    )

    // Evaluate core gene space coverage
    BUSCO (
        assembly_ch.map { sample, assembly -> [ sample + [ build: assembly.id ] , assembly.pri_fasta ] },
        busco_lineages,
        busco_lineage_path,
        []
    )
}