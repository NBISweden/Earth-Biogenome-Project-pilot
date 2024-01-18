include { HIFIASM         } from "$projectDir/modules/nf-core/hifiasm/main"
include { GFATOOLS_GFA2FA } from "$projectDir/modules/local/gfatools/gfa2fa"
include { GFASTATS        } from "$projectDir/modules/nf-core/gfastats/main"

workflow ASSEMBLE_HIFI {
    take:
    hifi_reads

    main:
        // Add build ID.
        reads_ch = hifi_reads
            .flatMap { meta, reads ->
                if (params.hifiasm) {
                    params.hifiasm.collect { key, value -> [ meta.deepMerge(
                        [
                            settings: [ hifiasm: [ id: key, args: value ] ],
                            assembly: [ assembler: 'hifiasm', stage: 'raw', id: key, build: "hifiasm-raw-$key" ]
                        ]
                        ), reads ]
                    }
                } else {
                    def key = "default"
                    [ [ meta.deepMerge(
                        [
                            settings: [ hifiasm: [ id: key, args: "" ] ],
                            assembly: [ assembler: 'hifiasm', stage: 'raw', id: key, build: "hifiasm-raw-$key" ]
                        ]
                        ), reads ]
                    ]
                }
            }
        HIFIASM(
            reads_ch,
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

    versions_ch = HIFIASM.out.versions.first()
        .mix( GFASTATS.out.versions.first() )
        .mix( GFATOOLS_GFA2FA.out.versions.first() )

    emit:
    assemblies = assemblies_ch
    versions = versions_ch

}