include { constructAssemblyRecord } from "$projectDir/modules/local/functions"
include { getPrimaryAssembly      } from "$projectDir/modules/local/functions"
include { joinByMetaKeys          } from "$projectDir/modules/local/functions"
include { BWAMEM2_INDEX           } from "$projectDir/modules/nf-core/bwamem2/index/main"
include { BWAMEM2_MEM             } from "$projectDir/modules/nf-core/bwamem2/mem/main"
include { SAMTOOLS_FAIDX          } from "$projectDir/modules/nf-core/samtools/faidx/main"
include { PAIRTOOLS_PARSE         } from "$projectDir/modules/nf-core/pairtools/parse/main"
include { PAIRTOOLS_SORT          } from "$projectDir/modules/nf-core/pairtools/sort/main"
include { PAIRTOOLS_MERGE         } from "$projectDir/modules/nf-core/pairtools/merge/main"
include { PAIRTOOLS_DEDUP         } from "$projectDir/modules/nf-core/pairtools/dedup/main"
include { PAIRTOOLS_SPLIT         } from "$projectDir/modules/nf-core/pairtools/split/main"
include { YAHS                    } from "$projectDir/modules/nf-core/yahs/main.nf"


workflow SCAFFOLD {
    take:
    ch_assemblies // [ meta, assembly ]
    ch_hic        // [ meta, hic-pairs ]

    main:
    ch_versions = Channel.empty()
    ch_scaffolded_assemblies = Channel.empty()


    ch_assemblies.view()
    ch_hic.view()

    BWAMEM2_INDEX ( getPrimaryAssembly(ch_assemblies) )

    SAMTOOLS_FAIDX (
        getPrimaryAssembly(ch_assemblies),
        [ [ ], [ ] ]
    )

    joinByMetaKeys(joinByMetaKeys(ch_hic, BWAMEM2_INDEX.out.index, 
        keySet: ['id','sample'],
        meta: 'rhs'
        ),
        getPrimaryAssembly(ch_assemblies),
        keySet: ['id','sample'],
        meta: 'rhs'
    ).multiMap{ meta, hic_reads, index, fasta -> 
        reads: [ meta, hic_reads ]
        index: [ meta, index ]
        fasta: [ meta, fasta ]
    }.set{ bwamem2_input } 

    BWAMEM2_MEM (bwamem2_input.reads, bwamem2_input.index, bwamem2_input.fasta, false)

    SAMTOOLS_FAIDX.out.fai.map{ meta, fai ->
        fai.splitCsv( sep: '\t', header: false )
            .collect{ row -> 
                row[ 0..1 ]
                .join('\t')
            }.join('\n')        
    }.collectFile()
    .set {chrom_sizes}

    PAIRTOOLS_PARSE(BWAMEM2_MEM.out.bam, chrom_sizes)

    PAIRTOOLS_SORT(PAIRTOOLS_PARSE.out.pairsam)

    PAIRTOOLS_SORT.out.sorted.groupTuple()
        .branch { meta, pairsam ->
            single: pairsam.size() == 1
                return [ meta, *pairsam ]
            multi: true
                return [ meta, pairsam ]
        }
        .set { pairtools_sorted }

    PAIRTOOLS_MERGE(pairtools_sorted.multi)

    PAIRTOOLS_DEDUP(
        PAIRTOOLS_MERGE.out.pairs.mix(
            pairtools_sorted.single
        )
    )



    

/*
    reads_plus_assembly_ch = combineByMetaKeys (
            ch_hic,
            ch_assemblies,
            keySet: ['id','sample'],
            meta: 'rhs'
        )
    reads_plus_assembly_ch
        // Add single_end for minimap module
        .flatMap { meta, reads, assembly -> reads instanceof List ?
            reads.collect{ [ meta + [ single_end: true ], it, assembly.pri_fasta ] }
            : [ [ meta + [ single_end: true ], reads, assembly.pri_fasta ] ] }
        .multiMap { meta, reads, assembly ->
            reads_ch: [ meta, reads ]
            assembly_ch: assembly
        }
        .set { input }
*/
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