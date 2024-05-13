/*
 * Workflow based around the DeepVariant tool to polish homozygous variants.
 * https://git.mpi-cbg.de/assembly/programs/polishing
 */
include { getPrimaryAssembly                      } from "$projectDir/modules/local/functions"
include { constructAssemblyRecord                 } from "$projectDir/modules/local/functions"
include { joinByMetaKeys                          } from "$projectDir/modules/local/functions"
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

    reads_plus_assembly_ch = combineByMetaKeys (
            ch_hifi,
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
            assembly_ch: [meta, assembly ]
        }
        .set { input }

    uniq_assembly_ch = getPrimaryAssembly(ch_assemblies)

    // index assembly file(s)
    SAMTOOLS_FAIDX (
        uniq_assembly_ch,
        [[],[]]
    )
    

    // split assembly into smaller chunks, this step just creates bed files 
    // that represent the assembly chunks, no sequence is split
    DVPOLISH_CHUNKFA (
        SAMTOOLS_FAIDX.out.fai
    )

    // create minimap2 index for assemblies
    DVPOLISH_PBMM2_INDEX (
        uniq_assembly_ch
    )

    // map reads with pbmm2 to complete assemblies (chunks are not used in that step)
    DVPOLISH_PBMM2_ALIGN (
        input.assembly_ch,
        input.reads_ch
    )

    def path_closure = {meta, files -> files.collect(){[meta, it ]}}

    DVPOLISH_PBMM2_ALIGN.out.bam
    .flatMap(path_closure)
    .combine(DVPOLISH_CHUNKFA.out.bed.flatMap(path_closure), by:0)
    .combine(DVPOLISH_PBMM2_ALIGN.out.bai.flatMap(path_closure), by:0)
    .multiMap { meta, bam, bai, bed ->
            bam_bai_ch: [ meta + [mergeID: bed.baseName], bam, bai]
            bed_ch: bed
        }
    .set { bam_bai_bed_ch }
    // split bam files according to bed file chunks 
    SAMTOOLS_VIEW (bam_bai_bed_ch.bam_bai_ch,   //  val(meta), path(input), path(index)
    [[],[]],                                    //  val(meta2), path(fasta)  NOT USED
    bam_bai_bed_ch.bed_ch)                      //  path qname IS REUSED, i.e. qname is removed(+patched) from the SAMTOOLS_VIEW process 
//                                                  and a bed file is provided, which needs to be liked into the CWD
//                                                  but its actually only be used via the ext.args = "-L ${index}"

    // index the splitted bam files 
    SAMTOOLS_INDEX_FILTER(SAMTOOLS_VIEW.out.bam)

    SAMTOOLS_VIEW.out.bam
    .groupTuple(by:0)
    .branch { meta, bam_list ->
        multiples: bam_list.size() > 1
        singleton: true
    }
    .set { bam_merge_ch }

    // in case multiple reads files are present, all corresping bam files 
    // the where splitted in the previous step need to be merged. key:bed file ID
    SAMTOOLS_MERGE(
        bam_merge_ch.multiples,
        [[],[]],
        [[],[]]
    )
    // index merged bam files
    SAMTOOLS_INDEX_MERGE(SAMTOOLS_MERGE.out.bam)

    bam_merge_ch.singleton
    .map { meta, bam -> [ meta, *bam ]} // the spread operator (*) flattens the bam list
    .join(SAMTOOLS_INDEX_FILTER.out.bai, by:0)
    .mix(SAMTOOLS_MERGE.out.bam
        .join(SAMTOOLS_INDEX_MERGE.out.bai, by:0)
    )
    .join(bam_bai_bed_ch
    .map { meta, bam, bai, bed -> [meta, bed]}
    .unique())
    .set {deepvariant_ch}
    // run deepvariant and the chunked bam files 
    DEEPVARIANT(
        deepvariant_ch,
        uniq_assembly_ch,
        SAMTOOLS_FAIDX.out.fai,
        [[],[]]     // tuple val(meta4), path(gzi)
    )

    DEEPVARIANT.out.vcf
    .join(DEEPVARIANT.out.vcf_tbi, by:0)
    .set { bcftools_view_ch }
    // filter vcf files for PASS and homozygous varinats
    // TODO add a minimim and maximum coverage filter ??? Needs to be tested
    BCFTOOLS_VIEW (
        bcftools_view_ch,
        [], // path(regions)
        [], // path(targets)
        [] // path(samples)
    )

    // index vcf file 
    TABIX_TABIX(
        BCFTOOLS_VIEW.out.vcf
    )

    // in case of multiple vcf files, merge them prior the consenus step
    BCFTOOLS_VIEW.out.vcf
    .map { meta, vcf -> [ meta - meta.subMap('mergeID'), vcf ] }
    .groupTuple(by:0)
    .set { filt_vcf_list_ch }

    TABIX_TABIX.out.tbi
    .map { meta, tbi -> [ meta - meta.subMap('mergeID'), tbi ] }
    .groupTuple(by:0)
    .set { filt_tbi_list_ch }

    filt_vcf_list_ch
    .join(filt_tbi_list_ch, by:0)
    .branch { meta, vcf_list, vcf_index_list ->
        multiples: vcf_list.size() > 1
        singleton: true
    }
    .set { vcf_merge_ch }

    // merge all vcf files 
    BCFTOOLS_MERGE(
        vcf_merge_ch.multiples,
        uniq_assembly_ch,
        SAMTOOLS_FAIDX.out.fai,
        [] // path(bed)
    )

    // index merged vcf file
    TABIX_TABIX_MERGED(
        BCFTOOLS_MERGE.out.merged_variants
    )

    vcf_plus_index_ch = vcf_merge_ch.singleton
    .map { meta, vcf, idx  -> [ meta, *vcf, *idx ] } // the spread operator (*) flattens the bam list
    .mix(BCFTOOLS_MERGE.out.merged_variants
        .join(TABIX_TABIX_MERGED.out.tbi)
    )

    vcf_plus_index_plus_assembly_ch = joinByMetaKeys (
        vcf_plus_index_ch,
        uniq_assembly_ch,
        keySet: ['id','single_end'],
        meta: 'lhs'
    )

    // create consensus sequence 
    BCFTOOLS_CONSENSUS(
        vcf_plus_index_plus_assembly_ch
    )

    ch_polished_assemblies = constructAssemblyRecord(
    BCFTOOLS_CONSENSUS.out.fasta
    )

    emit:
    assemblies = ch_polished_assemblies
}
