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
include { MERQURY as MERQURY_INPUT_ASM            } from "$projectDir/modules/nf-core/merqury/main"
include { MERQURY as MERQURY_POLISHED_ASM         } from "$projectDir/modules/nf-core/merqury/main"
include { DVPOLISH_CREATE_FINALASM                } from "$projectDir/modules/local/dvpolish/createFinalAsm"

/*
outline: 

|- create pbmm2 index for assembly (1)
|- create bed chunks for given assembly (1..n)
|- align all read files to full assembly (1..m)
    \- split each alignment file to contig according to bed chunks files (n*m)
    |- in case of multiple read files (therefore optional) merge all read files that belong to same assembly chunk (1..n)
    |- index merged alignment files (1..n)
    |- call variants with DeepVartiant (1..n)
    |- filter variants (PASS + homozygous) (1..n)
    |- create tabix index files  (1..n)
    |- merge all variants (1)
    |- create consensus sequence (1)
*/

workflow DVPOLISH {

    take:
    ch_assemblies // [ meta, assembly ]
    ch_hifi       // [ meta, hifi ]
    ch_meryl_hifi // [ meta, union.meryldb ] 

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
            reads_ch: [ meta + [ readID: reads.baseName ], reads ]
            assembly_ch: [ meta, assembly ]
        }
        .set { input }

    uniq_assembly_ch = getPrimaryAssembly(ch_assemblies)

    assembly_plus_meryl_ch = combineByMetaKeys (
            ch_meryl_hifi,
            uniq_assembly_ch,
            keySet: ['id','sample'],
            meta: 'rhs'
        )

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
        input.reads_ch,
        input.assembly_ch
    )

    combineByMetaKeys (
        DVPOLISH_PBMM2_ALIGN.out.bam_bai,
        DVPOLISH_CHUNKFA.out.bed.transpose(),
        keySet: ['sample','assembly'],
        meta: 'rhs'
    )
    .multiMap { meta, bam, bai, bed ->
        meta_bam_bai_ch:  [ meta + [ mergeID: bed.baseName ], bam, bai ]
        meta_bed_ch:      [ meta + [ mergeID: bed.baseName ], bed ]
        bed_ch:             bed
    }
    .set { alignment }
    
    // split bam files according to bed file chunks 
    SAMTOOLS_VIEW (alignment.meta_bam_bai_ch,
    [[],[]],                            
    alignment.bed_ch)                   

    // index the splitted bam files 
    SAMTOOLS_INDEX_FILTER(SAMTOOLS_VIEW.out.bam)

    SAMTOOLS_VIEW.out.bam
    .groupTuple(by:0)
    .branch { meta, bam_list ->
        multiples: bam_list.size() > 1
        singleton: true
    }
    .set { bam_merge_ch }

    // in case multiple reads files are present, all corresponding bam files 
    // that were splitted in the previous step need to be merged. key:bed file ID
    SAMTOOLS_MERGE(
        bam_merge_ch.multiples,
        [[],[]],
        [[],[]]
    )
    // index merged bam files
    SAMTOOLS_INDEX_MERGE(SAMTOOLS_MERGE.out.bam)

    bam_merge_ch.singleton
    .map { meta, bam -> [ meta, *bam ] } // the spread operator (*) flattens the bam list
    .join(SAMTOOLS_INDEX_FILTER.out.bai)
    .mix(SAMTOOLS_MERGE.out.bam
        .join(SAMTOOLS_INDEX_MERGE.out.bai)
    )
    .join(alignment.meta_bed_ch)
    .set { dv_bam_bai_bed_ch }

    asm_fai_ch = joinByMetaKeys (
        uniq_assembly_ch,
        SAMTOOLS_FAIDX.out.fai,
        keySet: ['sample','assembly'],
        meta: 'lhs'
    )

    combineByMetaKeys (
        dv_bam_bai_bed_ch,
        asm_fai_ch,
        keySet: ['sample','assembly'],
        meta: 'lhs'
    )
    .multiMap { meta, bam, bai, bed, fasta, fai ->
        bam_bai_bed_ch: [ meta, bam, bai, bed ]
        fasta_ch:       [ meta, fasta ]
        fai_ch:         [ meta, fai ]
    }
    .set { dv_input }

    // run deepvariant and the chunked bam files 
    DEEPVARIANT(
        dv_input.bam_bai_bed_ch,    // tuple val(meta), path(input), path(index), path(intervals)
        dv_input.fasta_ch,          // tuple val(meta2), path(fasta)
        dv_input.fai_ch,            // tuple val(meta3), path(fai)
        [[],[]]                     // tuple val(meta4), path(gzi)
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

    joinByMetaKeys (
        vcf_merge_ch.multiples,
        asm_fai_ch,
        keySet: ['sample','assembly'],
        meta: 'lhs'
    )
    .multiMap { meta, vcfs, tbis, fasta, fai ->
        vcf_tbis_ch:    [ meta, vcfs, tbis ]
        fasta_ch:       [ meta, fasta ]
        fai_ch:         [ meta, fai ]
    }
    .set { bcf_input }

    // merge all vcf files 
    BCFTOOLS_MERGE(
        bcf_input.vcf_tbis_ch,
        bcf_input.fasta_ch,
        bcf_input.fai_ch,
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
        keySet: ['sample','assembly'],
        meta: 'lhs'
    )

    // create consensus sequence 
    BCFTOOLS_CONSENSUS(
        vcf_plus_index_plus_assembly_ch
    )

    // run merqury on input assembly 
    MERQURY_INPUT_ASM(assembly_plus_meryl_ch)

    // run merqury on polished assembly
    polishedASM_plus_meryl_ch = combineByMetaKeys (
        ch_meryl_hifi,
        BCFTOOLS_CONSENSUS.out.fasta,
        keySet: ['id','sample'],
        meta: 'rhs'
    )
    MERQURY_POLISHED_ASM(polishedASM_plus_meryl_ch)

    unpolASM_merqQV_ch = combineByMetaKeys(
        uniq_assembly_ch,
        MERQURY_INPUT_ASM.out.scaffold_qv,
        keySet: ['id','assembly'],
        meta: 'rhs'
    )

    polASM_merqQV_ch = combineByMetaKeys(
        BCFTOOLS_CONSENSUS.out.fasta,
        MERQURY_POLISHED_ASM.out.scaffold_qv,
        keySet: ['id','assembly'],
        meta: 'rhs'
    )

    combineByMetaKeys(
        unpolASM_merqQV_ch,
        polASM_merqQV_ch,
        keySet: ['id','assembly'],
        meta: 'rhs'
    )
    .multiMap { meta, unpol_asm, unpol_qv, pol_asm, pol_qv ->
        unpolASM_qv_ch: [ meta, unpol_asm, unpol_qv ]
        polASM_qv_ch: [ meta, pol_asm, pol_qv ]
    }
    .set { createFinalAsm }

    DVPOLISH_CREATE_FINALASM(
        createFinalAsm.unpolASM_qv_ch,
        createFinalAsm.polASM_qv_ch,
    )

    // TODO 
    // 3. publish result

    ch_polished_assemblies = constructAssemblyRecord(
    DVPOLISH_CREATE_FINALASM.out.fasta_gz
    )

    emit:
    assemblies = ch_polished_assemblies
}
