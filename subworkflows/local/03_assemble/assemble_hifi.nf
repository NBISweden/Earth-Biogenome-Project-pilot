include { combineByMetaKeys } from "../../../modules/local/functions"
include { HIFIASM           } from "../../../modules/nf-core/hifiasm/main"
include { GFATOOLS_GFA2FA   } from "../../../modules/local/gfatools/gfa2fa"
include { deepMergeMaps     } from "../../../modules/local/functions"

workflow ASSEMBLE_HIFI {
    take:
    hifi_reads     // [ meta, fastx ]
    hic_reads      // [ meta, hic]

    main:
    // Add build ID.
    reads_ch = hifi_reads
        .flatMap { meta, reads ->
            if (params.hifiasm) {
                params.hifiasm.collect { key, value -> [ deepMergeMaps(meta,
                    [
                        settings: [ hifiasm: [ id: key, args: value ] ],
                        assembly: [ assembler: 'hifiasm', stage: 'raw', id: key, build: "hifiasm-raw-$key" ]
                    ]
                    ), reads, [] ]
                }
            } else {
                def key = "default"
                [ [ deepMergeMaps(meta,
                    [
                        settings: [ hifiasm: [ id: key, args: "" ] ],
                        assembly: [ assembler: 'hifiasm', stage: 'raw', id: key, build: "hifiasm-raw-$key" ]
                    ]
                    ), reads, [] ]
                ]
            }
        }
    if( params.hifiasm_with_hic ){
        ch_hifiasm_input = combineByMetaKeys(
            reads_ch,
            hic_reads.map{ meta, hic -> tuple(meta, hic.first(), hic.last()) }.groupTuple(),
            keySet: ['id','sample'],
            meta: 'lhs'
        )
        .multiMap { meta, hifi, ul_ont, hic1, hic2 ->
            long_reads : tuple(meta, hifi, ul_ont)
            short_reads: tuple(meta, hic1, hic2)
        }
        HIFIASM(
            ch_hifiasm_input.long_reads,
            [[],[],[]], // meta, paternal k-mers, maternal k-mers
            ch_hifiasm_input.short_reads, // meta, Hi-C r1, Hi-C r2
            [[], []],   // meta, bin files
        )
    } else {
        HIFIASM(
            reads_ch,
            [[],[],[]], // meta, paternal k-mers, maternal k-mers
            [[],[],[]], // meta, Hi-C r1, Hi-C r2
            [[], []],   // meta, bin files
        )
    }
    raw_assembly_ch = params.use_phased ? HIFIASM.out.hap1_contigs.mix( HIFIASM.out.hap2_contigs ) : HIFIASM.out.primary_contigs
    GFATOOLS_GFA2FA( raw_assembly_ch )

    gfa_ch = params.use_phased ?
        HIFIASM.out.hap1_contigs
            .join( HIFIASM.out.hap2_contigs )
            .map { meta, pgfa, mgfa -> [ meta, [ pgfa, mgfa ] ] } :
        HIFIASM.out.primary_contigs
    assemblies_ch = GFATOOLS_GFA2FA.out.fasta.groupTuple( sort: { fa -> fa.name } )
        .join( gfa_ch )
        .map { meta, fasta, gfa ->
            [ meta, meta.assembly + (
                params.use_phased ?
                [
                    pri_fasta: fasta[0],
                    alt_fasta: fasta[1],
                    pri_gfa: gfa[0],
                    alt_gfa: gfa[1]
                ] :
                [
                    pri_fasta: fasta[0],
                    alt_fasta: null,
                    pri_gfa: gfa,
                    alt_gfa: null
                ]
            ) ]
        }
        .dump( tag: "Assemblies: Pre-purge", pretty: true )

    HIFIASM.out.log
        .map { _meta, log -> log }
        .set { logs }
    versions_ch = HIFIASM.out.versions.first()
        .mix( GFATOOLS_GFA2FA.out.versions.first() )

    emit:
    assemblies = assemblies_ch
    logs
    versions = versions_ch
}