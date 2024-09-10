include { constructAssemblyRecord } from "$projectDir/modules/local/functions"
include { BWAMEM2_INDEX           } from "$projectDir/modules/nf-core/bwamem2/index/main"
include { BWAMEM2_MEM             } from "$projectDir/modules/nf-core/bwamem2/mem/main"
include { YAHS                    } from "$projectDir/modules/nf-core/yahs/main.nf"


workflow SCAFFOLD {
    take:
    ch_assemblies // [ meta, assembly ]
    ch_hic        // [ meta, hic-pairs ]

    main:
    ch_versions = Channel.empty()

    // TODO Fill in workflow
    // Scaffold Assembly using HiC

    // Map HiC to assembly
    // BWAMEM2_INDEX ( ch_assemblies )
    // Use join and multiMap operators to ensure correctly paired input.
    // BWAMEM2_MEM ( index, hic )

    // Yahs
    // YAHS( channels )
    //

    emit:
    assemblies = ch_scaffolded_assemblies
    versions   = ch_versions
}