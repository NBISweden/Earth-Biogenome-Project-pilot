include { constructAssemblyRecord               } from "../../../modules/local/functions"
include { getPrimaryAssembly                    } from "../../../modules/local/functions"
include { combineByMetaKeys                     } from "../../../modules/local/functions"
include { joinByMetaKeys                        } from "../../../modules/local/functions"
include { BWAMEM2_INDEX                         } from "../../../modules/nf-core/bwamem2/index/main"
include { BWAMEM2_MEM                           } from "../../../modules/nf-core/bwamem2/mem/main"
include { FILTER_FIVE_END                       } from "../../../modules/local/hic_curation/filter_five_end"
include { TWOREADCOMBINER_FIXMATE_SORT          } from "../../../modules/local/hic_curation/tworeadcombiner_fixmate_sort"
include { BIOBAMBAM_BAMMARKDUPLICATES2          } from "../../../modules/nf-core/biobambam/bammarkduplicates2/main"
include { SAMTOOLS_FAIDX                        } from "../../../modules/nf-core/samtools/faidx/main"
include { SAMTOOLS_MERGE                        } from "../../../modules/nf-core/samtools/merge/main"
include { SAMTOOLS_SORT                         } from "../../../modules/nf-core/samtools/sort/main"
include { YAHS                                  } from "../../../modules/nf-core/yahs/main"

workflow SCAFFOLD_ARIMA_BWA {
    take:
    ch_assemblies // [ meta, assembly ]
    ch_hic        // [ meta, hic-pairs ]

    main:

    ch_versions = Channel.empty()

    BWAMEM2_INDEX(getPrimaryAssembly(ch_assemblies))
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

    SAMTOOLS_FAIDX(
        getPrimaryAssembly(ch_assemblies),
        [[], []],
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    combineByMetaKeys(
        combineByMetaKeys(
            ch_hic,
            BWAMEM2_INDEX.out.index,
            keySet: ['id', 'sample'],
            meta: 'merge',
        ),
        getPrimaryAssembly(ch_assemblies),
        keySet: ['id', 'sample'],
        meta: 'merge',
    ).transpose(by: 1).multiMap { meta, hic_reads, index, fasta ->
        reads: [meta, hic_reads]
        index: [meta, index]
        fasta: [meta, fasta]
    }.set { bwamem2_input }

    BWAMEM2_MEM(
        bwamem2_input.reads,
        bwamem2_input.index,
        bwamem2_input.fasta,
        false,
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    // filter alignments
    FILTER_FIVE_END(BWAMEM2_MEM.out.bam)
    ch_versions = ch_versions.mix(FILTER_FIVE_END.out.versions.first())

    // combine reads
    FILTER_FIVE_END.out.bam
        .map { meta, bam -> [meta.subMap('id', 'sample', 'assembly', 'pair_id'), bam] }
        .groupTuple(sort: { a, b -> a.name <=> b.name })
        .set { combine_input }

    TWOREADCOMBINER_FIXMATE_SORT(combine_input)
    ch_versions = ch_versions.mix(TWOREADCOMBINER_FIXMATE_SORT.out.versions.first())

    // merge bam files in case multiple HIC paired-end libraries are present
    TWOREADCOMBINER_FIXMATE_SORT.out.bam
        .map { meta, bam_list -> [meta - meta.subMap('pair_id'), bam_list] }
        .groupTuple()
        .branch { meta, bam_list ->
            multiples: bam_list.size() > 1
            singleton: true
        }
        .set { merge_bam }

    SAMTOOLS_MERGE(
        merge_bam.multiples,
        [[], []],
        [[], []],
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    merge_bam.singleton
        .map { meta, bam -> [meta, bam.first()] }
        .mix(SAMTOOLS_MERGE.out.bam)
        .set { dedup_bam }

    // dedupliucate bam file
    BIOBAMBAM_BAMMARKDUPLICATES2(dedup_bam)
    ch_versions = ch_versions.mix(BIOBAMBAM_BAMMARKDUPLICATES2.out.versions.first())

    // convert bam to sorted bed file
    SAMTOOLS_SORT(
        BIOBAMBAM_BAMMARKDUPLICATES2.out.bam,  // [ meta,  bam ]
        [[],[]]                                // [ meta2, reference_fasta ] --> not used 
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    combineByMetaKeys(
        combineByMetaKeys(
            SAMTOOLS_SORT.out.bam,
            getPrimaryAssembly(ch_assemblies),
            keySet: ['id', 'sample'],
            meta: 'rhs',
        ),
        SAMTOOLS_FAIDX.out.fai,
        keySet: ['id', 'sample'],
        meta: 'rhs',
    ).multiMap { meta, bam, fasta, fai ->
        bam: [meta, bam]
        fasta: [fasta]
        fai: [fai]
    }.set { yahs_input }

    YAHS(
        yahs_input.bam,
        yahs_input.fasta,
        yahs_input.fai,
    )
    ch_versions = ch_versions.mix(YAHS.out.versions.first())

    // Consensus case:
    // Preserve haplotigs from purge dups
    ch_scaff_and_alt = ch_assemblies
        .map { meta, assembly -> [meta, assembly.alt_fasta] }
        .mix(YAHS.out.scaffolds_fasta)
    ch_scaffolded_assemblies = constructAssemblyRecord(ch_scaff_and_alt, false)

    emit:
    assemblies = ch_scaffolded_assemblies
    logs       = []
    versions   = ch_versions
}
