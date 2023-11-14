include { HIFIASM         } from "$projectDir/modules/nf-core/hifiasm/main"
include { GFATOOLS_GFA2FA } from "$projectDir/modules/local/gfatools/gfa2fa"
include { GFASTATS        } from "$projectDir/modules/nf-core/gfastats/main"

workflow ASSEMBLE_HIFI {
    take:
    hifi_reads

    main:
        // Need to include a build ID here. 
        HIFIASM(
            hifi_reads.flatMap { meta, reads -> params.hifiasm ? params.hifiasm.collect { [ meta, reads ] } : [ [ meta, reads ] ] },
            [], // paternal k-mers
            [], // maternal k-mers
            [], // Hi-C r1
            []  // Hi-C r2
        )
        raw_assembly_ch = params.use_phased ? HIFIASM.out.paternal_contigs.mix( HIFIASM.out.maternal_contigs ) : HIFIASM.out.processed_contigs
        GFASTATS(
            raw_assembly_ch,
            "gfa",   // output format
            "",      // genome size
            "",      // target
            [],      // AGP file
            [],      // include bed
            [],      // exclude bed
            []       // SAK instructions
        )
        GFATOOLS_GFA2FA( raw_assembly_ch )

        gfa_ch = params.use_phased ? 
            HIFIASM.out.paternal_contigs
                .join( HIFIASM.out.maternal_contigs )
                .map { meta, pgfa, mgfa -> [ meta, [ pgfa, mgfa ] ] } : 
            HIFIASM.out.processed_contigs
        assemblies_ch = GFATOOLS_GFA2FA.out.fasta.groupTuple( sort: { it.name } )
            .join( gfa_ch )
            .map { meta, fasta, gfa -> 
                [ meta, 
                    params.use_phased ? 
                    [ 
                        id: 'hifiasm_phased', 
                        pri_fasta: fasta[0],
                        alt_fasta: fasta[1],
                        pri_gfa: gfa[0],
                        alt_gfa: gfa[1]
                    ] :
                    [
                        id: 'hifiasm_consensus', 
                        pri_fasta: fasta[0],
                        alt_fasta: null,
                        pri_gfa: gfa,
                        alt_gfa: null
                    ]
                ]
            }
            .dump( tag: "Assemblies: Pre-purge", pretty: true )

        // Find mitochondria
            // Need to check options to mitohifi modules.

    versions_ch = HIFIASM.out.versions.first()
        .mix( GFASTATS.out.versions.first() )
        .mix( GFATOOLS_GFA2FA.out.versions.first() )

    emit:
    assemblies = assemblies_ch
    versions = versions_ch

}