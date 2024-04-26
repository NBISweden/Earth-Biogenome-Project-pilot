/*
 * Workflow based around the DeepVariant tool to polish homozygous variants.
 * https://git.mpi-cbg.de/assembly/programs/polishing
 */

include { combineByMetaKeys                       } from "$projectDir/modules/local/functions"
include { DVPOLISH_CHUNKFA                        } from "$projectDir/modules/local/dvpolish/chunkfa"
include { DVPOLISH_PBMM2_INDEX                    } from "$projectDir/modules/local/dvpolish/pbmm2_index"
include { DVPOLISH_PBMM2_ALIGN                    } from "$projectDir/modules/local/dvpolish/pbmm2_align"
include { SAMTOOLS_FAIDX                          } from "$projectDir/modules/nf-core/samtools/faidx/main"
include { SAMTOOLS_VIEW                           } from "$projectDir/modules/nf-core/samtools/view/main"
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTER } from "$projectDir/modules/nf-core/samtools/index/main"
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGE  } from "$projectDir/modules/nf-core/samtools/index/main"
include { SAMTOOLS_MERGE                          } from "$projectDir/modules/nf-core/samtools/merge/main"
include { DEEPVARIANT                             } from "$projectDir/modules/nf-core/deepvariant/main"
include { BCFTOOLS_VIEW                           } from "$projectDir/modules/nf-core/bcftools/view/main"
include { TABIX_TABIX as TABIX_TABIX              } from "$projectDir/modules/nf-core/tabix/tabix/main"
include { TABIX_TABIX as TABIX_TABIX_MERGED       } from "$projectDir/modules/nf-core/tabix/tabix/main"
include { BCFTOOLS_MERGE                          } from "$projectDir/modules/nf-core/bcftools/merge/main"
include { BCFTOOLS_CONSENSUS                      } from "$projectDir/modules/nf-core/bcftools/consensus/main"

workflow DVPOLISH {

    take:
    ch_assemblies // [ meta, assembly ]
    ch_hifi       // [ meta, hifi ]

    main:
    //
    // MODULE: Run SAMTOOLS_FAIDX
    //
    SAMTOOLS_FAIDX (
        ch_assemblies,
        [[],[]]
    )

    //
    // MODULE: Run SPLIT_FA
    //
    DVPOLISH_CHUNKFA (
        SAMTOOLS_FAIDX.out.fai
    )

    //
    // MODULE: PBMM2_INDEX
    //
    DVPOLISH_PBMM2_INDEX (
        ch_assemblies
    )

    //
    // MODULE: PBMM2_ALIGN
    //
    DVPOLISH_PBMM2_ALIGN (
        ch_assemblies,
        ch_hifi
    )

    //
    // MODULE: Run SAMTOOLS_VIEW
    //

    def path_closure = {meta, files -> files.collect(){[meta, it ]}}

    DVPOLISH_PBMM2_ALIGN.out.bam
    .flatMap(path_closure)
    .combine(DVPOLISH_CHUNKFA.out.bed.flatMap(path_closure), by:0)
    .map { map, bam, bed -> [ map + [mergeID: bed.baseName], bam, bed]
    }
    .set { bam_bed_ch }

    DVPOLISH_PBMM2_ALIGN.out.bai
    .flatMap(path_closure)
    .combine(DVPOLISH_CHUNKFA.out.bed.flatMap(path_closure), by:0)
    .map { meta, bai, bed -> bai }
    .set { bai_ch }

    SAMTOOLS_VIEW (bam_bed_ch, [[],[]], bai_ch)

    //
    // MODULE: Run SAMTOOLS_INDEX
    //
    SAMTOOLS_INDEX_FILTER(SAMTOOLS_VIEW.out.bam)

    SAMTOOLS_VIEW.out.bam
    .groupTuple(by:0)
    .branch { meta, bam_list ->
        merge: bam_list.size() > 1
        link: true
    }
    .set { bam_merge_ch }

    //
    // MODULE: Run SAMTOOLS_MERGE
    //
    SAMTOOLS_MERGE(
        bam_merge_ch.merge,
        [[],[]],
        [[],[]]
    )
    SAMTOOLS_INDEX_MERGE(SAMTOOLS_MERGE.out.bam)

    //
    // MODULE: Run Deepvariant
    //
    bam_merge_ch.link
    .map { meta, bam -> [ meta, *bam ]} // the spread operator (*) flattens the bam list
    .join(SAMTOOLS_INDEX_FILTER.out.bai, by:0)
    .mix(SAMTOOLS_MERGE.out.bam
        .join(SAMTOOLS_INDEX_MERGE.out.bai, by:0)
    )
    .join(bam_bed_ch
    .map { meta, bam, bed -> [meta, bed]}
    .unique())
    .set {deepvariant_ch}

    DEEPVARIANT(
        deepvariant_ch,
        ch_assemblies,
        SAMTOOLS_FAIDX.out.fai,
        [[],[]]     // tuple val(meta4), path(gzi)
    )

    DEEPVARIANT.out.vcf
    .join(DEEPVARIANT.out.vcf_tbi, by:0)
    .set { bcftools_view_ch }

    bcftools_view_ch.view()

    //
    // MODULE: bcftools view [predefined filter options TODO: needs tests] 
    //
    BCFTOOLS_VIEW (
        bcftools_view_ch,
        [], // path(regions)
        [], // path(targets)
        [] // path(samples)
    )

    //
    // MODULE: tabix
    //
    TABIX_TABIX(
        BCFTOOLS_VIEW.out.vcf
    )

    // in case of multiple vcf files, merge them prior the consenus step
    BCFTOOLS_VIEW.out.vcf
    .map { meta, vcf -> [ meta.subMap('id', 'single_end'), vcf ] }
    .groupTuple(by:0)
    .set { filt_vcf_list_ch }

    TABIX_TABIX.out.tbi
    .map { meta, vcf -> [ meta.subMap('id', 'single_end'), vcf ] }
    .groupTuple(by:0)
    .set { filt_tbi_list_ch }

    filt_vcf_list_ch
    .join(filt_tbi_list_ch, by:0)
    .branch { meta, vcf_list, vcf_index_list ->
        merge: vcf_list.size() > 1
        other: true
    }
    .set { vcf_merge_ch }

    BCFTOOLS_MERGE(
        vcf_merge_ch.merge,
        ch_assemblies,
        SAMTOOLS_FAIDX.out.fai,
        [] // path(bed)
    )

    TABIX_TABIX_MERGED(
        BCFTOOLS_MERGE.out.merged_variants
    )

    vcf_merge_ch.other
    .map { meta, vcf, idx  -> [ meta, *vcf, *idx ] } // the spread operator (*) flattens the bam list
    .mix(BCFTOOLS_MERGE.out.merged_variants
        .join(TABIX_TABIX_MERGED.out.tbi)
    )
    .join(ch_assemblies)
    .set { bcftools_consensus_ch }

    vcf_merge_ch.other.view { "vcf_merge_ch.other: " + it}
    bcftools_consensus_ch.view { "bcftools_consensus_ch: " + it }
    ch_assemblies.view { "input.assembly_ch: " + it}    

    BCFTOOLS_CONSENSUS(
        bcftools_consensus_ch
    )

    ch_polished_assemblies = BCFTOOLS_CONSENSUS.out.fasta

    emit:
    assemblies = ch_polished_assemblies
}
