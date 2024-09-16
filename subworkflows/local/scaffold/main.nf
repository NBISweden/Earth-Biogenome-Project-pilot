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

    joinByMetaKeys(PAIRTOOLS_DEDUP.out.pairs, 
        getPrimaryAssembly(ch_assemblies),
        keySet: ['id','sample'],
        meta: 'rhs'
    ).multiMap{ meta, pairs, fasta -> 
        pairs: [ meta, pairs ]
        fasta: [ meta, fasta ]
    }.set { pairtools_split_input }

    PAIRTOOLS_SPLIT(pairtools_split_input.pairs,
        pairtools_split_input.fasta,
        true    // sort_bam
    )
    

    joinByMetaKeys(
        joinByMetaKeys(PAIRTOOLS_SPLIT.out.bam, 
            getPrimaryAssembly(ch_assemblies), 
            keySet: ['id','sample'],
            meta: 'rhs'
        ),
        SAMTOOLS_FAIDX.out.fai,
        keySet: ['id','sample'],
        meta: 'rhs'
    ).multiMap{ meta, bam, fasta, fai -> 
        bam:   [ meta, bam ]
        fasta: [ fasta ]
        fai:   [ fai ]
    }.set{ yahs_input } 

    YAHS(yahs_input.bam,
        yahs_input.fasta,
        yahs_input.fai
    )

    ch_scaffolded_assemblies = constructAssemblyRecord( YAHS.out.scaffolds_fasta )

    emit:
    assemblies = ch_scaffolded_assemblies
    versions   = ch_versions
}