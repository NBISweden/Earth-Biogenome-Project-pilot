include { BUSCO               } from "$projectDir/modules/nf-core/busco/main"
include { MERQURYFK_MERQURYFK } from "$projectDir/modules/local/merquryfk/merquryfk"
// include { INSPECTOR           } from "$projectDir/modules/local/inspector/inspector"

workflow EVALUATE_ASSEMBLY {

    take:
    assembly_ch        // input type: [ meta, [ id:'assemblerX_build1', pri_fasta: '/path/to/primary_asm', alt_fasta: '/path/to/alternate_asm' ] ]
    reads_ch           // input type: [ meta, [ 'path/to/reads' ] ]
    fastk_db           // input type: [ meta, [ 'path/to/reads.hist' ], [ '/path/to/reads.ktab' ] ]
    reference_ch       // optional: file( reference_genome ) for comparison
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

    // Evaluate core gene space coverage
    busco_input = assembly_ch.map { sample, assembly -> [ sample + [ build: assembly.id ] , assembly.pri_fasta ] }
        .flatMap { meta, asm -> meta.settings?.busco?.lineages ? meta.settings.busco.lineages.tokenize(',').collect{ [ meta, asm, it ] } : [ [ meta, asm, 'auto' ] ] }
        .dump(tag: 'BUSCO')
        .multiMap { meta, asm, line ->
            assembly_ch: [ meta, asm ]
            lineage_ch: line
        }
    BUSCO (
        busco_input.assembly_ch,
        'genome',
        busco_input.lineage_ch,
        busco_lineage_path,
        []
    )
}