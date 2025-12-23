/*
 * Workflow based around the DeepVariant tool to polish homozygous variants.
 * https://git.mpi-cbg.de/assembly/programs/polishing
 */
include { getPrimaryAssembly                } from "../../../modules/local/functions"
include { constructAssemblyRecord           } from "../../../modules/local/functions"
include { joinByMetaKeys                    } from "../../../modules/local/functions"
include { combineByMetaKeys                 } from "../../../modules/local/functions"
include { DVPOLISH_CHUNKFA                  } from "../../../modules/local/dvpolish/chunkfa"
include { DVPOLISH_PBMM2_INDEX              } from "../../../modules/local/dvpolish/pbmm2_index"
include { DVPOLISH_PBMM2_ALIGN              } from "../../../modules/local/dvpolish/pbmm2_align"
include { SAMTOOLS_FAIDX                    } from "../../../modules/nf-core/samtools/faidx/main"
include { SAMTOOLS_VIEW                     } from "../../../modules/nf-core/samtools/view/main"
include { SAMTOOLS_INDEX                    } from "../../../modules/nf-core/samtools/index/main"
include { SAMTOOLS_MERGE                    } from "../../../modules/nf-core/samtools/merge/main"
include { DEEPVARIANT                       } from "../../../modules/nf-core/deepvariant/rundeepvariant/main"
include { BCFTOOLS_VIEW                     } from "../../../modules/nf-core/bcftools/view/main"
include { TABIX_TABIX as TABIX_TABIX        } from "../../../modules/nf-core/tabix/tabix/main"
include { TABIX_TABIX as TABIX_TABIX_MERGED } from "../../../modules/nf-core/tabix/tabix/main"
include { BCFTOOLS_MERGE                    } from "../../../modules/nf-core/bcftools/merge/main"
include { BCFTOOLS_CONSENSUS                } from "../../../modules/nf-core/bcftools/consensus/main"
include { MERQURY as MERQURY_INPUT_ASM      } from "../../../modules/nf-core/merqury/main"
include { MERQURY as MERQURY_POLISHED_ASM   } from "../../../modules/nf-core/merqury/main"
include { DVPOLISH_CREATE_FINALASM          } from "../../../modules/local/dvpolish/createFinalAsm"

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
    ch_meryl_hifi // [ meta, union.meryl ]

    main:
    ch_logs     = channel.empty()
    ch_versions = channel.empty()

    // Create input channel for assembly-level tasks
    uniq_assembly_ch = getPrimaryAssembly(ch_assemblies)

    // Index assembly file(s)
    SAMTOOLS_FAIDX (
        uniq_assembly_ch,
        [ [] , [] ]
    )

    // Generate BED intervals partitioning the assembly into fixed-size regions
    DVPOLISH_CHUNKFA (
        SAMTOOLS_FAIDX.out.fai
    )

    // create minimap2 index for assemblies
    DVPOLISH_PBMM2_INDEX (
        uniq_assembly_ch
    )

    // Create input channel for read-level tasks
    reads_plus_assembly_ch = combineByMetaKeys (
        ch_hifi,
        ch_assemblies,
        keySet: [ 'id', 'sample' ],
        meta: 'rhs'
    )
    reads_plus_assembly_plus_index_ch = combineByMetaKeys (
        reads_plus_assembly_ch,
        DVPOLISH_PBMM2_INDEX.out.index,
        keySet: [ 'id', 'sample' ],
        meta: 'lhs'
    )
    reads_plus_assembly_plus_index_ch
        .flatMap { meta, reads, assembly, index -> reads instanceof List ?
            reads.collect{ [ meta, it, assembly.pri_fasta, index ] }
            : [ [ meta, reads, assembly.pri_fasta, index ] ] }
        .multiMap { meta, reads, assembly, index ->
            reads_ch: [ meta + [ readID: reads.baseName ], reads ]
            assembly_ch: [ meta, assembly ]
            index_ch: [ meta, index ]
        }
        .set { input }

    // Map reads with pbmm2 to complete assemblies
    DVPOLISH_PBMM2_ALIGN (
        input.reads_ch,
        input.index_ch
    )

    // Create channel combining whole genome alignment with bed chunk files
    combineByMetaKeys (
        DVPOLISH_PBMM2_ALIGN.out.bam_bai,
        DVPOLISH_CHUNKFA.out.bed.transpose(),
        keySet: [ 'sample', 'assembly' ],
        meta: 'rhs'
    )
    .multiMap { meta, bam, bai, bed ->
        meta_bam_bai_ch:  [ meta + [ mergeID: bed.baseName ], bam, bai ]
        meta_bed_ch:      [ meta + [ mergeID: bed.baseName ], bed ]
        bed_ch:             bed
    }
    .set { alignment } // Note on meta: readID dropped, mergeID added

    // Split alignment bam by BED coordinates
    SAMTOOLS_VIEW (
        alignment.meta_bam_bai_ch,
        [ [],[] ],
        alignment.bed_ch
    )

    // Groups BAMs on meta (BAMs from different readsets covering the same BED region are grouped)
    SAMTOOLS_VIEW.out.bam
        .groupTuple(by:0)
        .branch { _meta, bam_list ->
            multiples: bam_list.size() > 1 // Multiple read files aligned to same region
            singleton: true                // Single read file aligned to region
        }
        .set { bam_merge_ch }

    // Merge BAMs aligned to the same region (originating from different input readsets). key:bed file ID
    SAMTOOLS_MERGE (
        bam_merge_ch.multiples,
        [ [], [] ],
        [ [], [] ]
    )

    // Mix BAMs for indexing
    bam_merge_ch.singleton
        .map { meta, bam -> [ meta, bam[0] ] }
        .mix( SAMTOOLS_MERGE.out.bam )
        .set { bams_to_index }

    // Index BAMS
    SAMTOOLS_INDEX(
        bams_to_index
    )

    // Join indexed bams with corresponding BED files
    bams_to_index
        .join( SAMTOOLS_INDEX.out.bai )
        .join( alignment.meta_bed_ch )
        .set { dv_bam_bai_bed_ch }

    // Join assembly with samtools index
    asm_fai_ch = joinByMetaKeys (
        uniq_assembly_ch,
        SAMTOOLS_FAIDX.out.fai,
        keySet: [ 'sample', 'assembly' ],
        meta: 'lhs'
    )
    // Create DeepVariant input channel
    combineByMetaKeys (
        dv_bam_bai_bed_ch,
        asm_fai_ch,
        keySet: ['sample', 'assembly' ],
        meta: 'lhs'
    )
    .multiMap { meta, bam, bai, bed, fasta, fai ->
        bam_bai_bed_ch: [ meta, bam, bai, bed ]
        fasta_ch:       [ meta, fasta ]
        fai_ch:         [ meta, fai ]
    }
    .set { dv_input }

    // run deepvariant on the chunked bam files
    DEEPVARIANT(
        dv_input.bam_bai_bed_ch,    // tuple val(meta), path(input), path(index), path(intervals)
        dv_input.fasta_ch,          // tuple val(meta2), path(fasta)
        dv_input.fai_ch,            // tuple val(meta3), path(fai)
        [ [], [] ]                  // tuple val(meta4), path(gzi)
    )

    // Prepare channel with indexed vcf
    DEEPVARIANT.out.vcf
        .join(DEEPVARIANT.out.vcf_tbi, by:0)
        .set { bcftools_view_ch }

    // filter vcf files for PASS and homozygous variants
    // TODO add a minimim and maximum coverage filter ??? Needs to be tested
    BCFTOOLS_VIEW (
        bcftools_view_ch,
        [], // path(regions)
        [], // path(targets)
        [] // path(samples)
    )

    // Index vcf files
    TABIX_TABIX(
        BCFTOOLS_VIEW.out.vcf
    )

    // in case of multiple vcf files, merge them prior to the consensus step
    BCFTOOLS_VIEW.out.vcf
        .map { meta, vcf -> [ meta - meta.subMap( 'mergeID' ), vcf ] }
        .groupTuple(by:0)
        .set { filt_vcf_list_ch }
    TABIX_TABIX.out.tbi
        .map { meta, tbi -> [ meta - meta.subMap( 'mergeID' ), tbi ] }
        .groupTuple( by:0 )
        .set { filt_tbi_list_ch }
    filt_vcf_list_ch
        .join( filt_tbi_list_ch, by:0 )
        .branch { _meta, vcf_list, _vcf_index_list ->
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

    // Prepare input for consensus step
    vcf_plus_index_ch = vcf_merge_ch.singleton
        .map { meta, vcf, idx  -> [ meta, vcf[0], idx[0] ] }
        .mix(BCFTOOLS_MERGE.out.merged_variants
        .join(TABIX_TABIX_MERGED.out.tbi)
    )
    vcf_plus_index_plus_assembly_ch = joinByMetaKeys (
        vcf_plus_index_ch,
        uniq_assembly_ch,
        keySet: [ 'sample', 'assembly' ],
        meta: 'lhs'
    )

    // create consensus sequence
    BCFTOOLS_CONSENSUS(
        vcf_plus_index_plus_assembly_ch
    )

    // run merqury on input assembly
    MERQURY_INPUT_ASM (
        combineByMetaKeys (
            ch_meryl_hifi,
            uniq_assembly_ch,
            keySet: [ 'id', 'sample' ],
            meta: 'rhs'
        )
    )

    // run merqury on polished assembly
    MERQURY_POLISHED_ASM (
        combineByMetaKeys (
            ch_meryl_hifi,
            BCFTOOLS_CONSENSUS.out.fasta,
            keySet: [ 'id', 'sample' ],
            meta: 'rhs'
        )
    )

    // combine merqury results & assemblies
    unpolASM_merqQV_ch = combineByMetaKeys(
        uniq_assembly_ch,
        MERQURY_INPUT_ASM.out.scaffold_qv,
        keySet: [ 'id', 'assembly' ],
        meta: 'rhs'
    )
    polASM_merqQV_ch = combineByMetaKeys(
        BCFTOOLS_CONSENSUS.out.fasta,
        MERQURY_POLISHED_ASM.out.scaffold_qv,
        keySet: [ 'id', 'assembly' ],
        meta: 'rhs'
    )
    combineByMetaKeys(
        unpolASM_merqQV_ch,
        polASM_merqQV_ch,
        keySet: [ 'id', 'assembly' ],
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

    ch_polished_assemblies = constructAssemblyRecord(
        DVPOLISH_CREATE_FINALASM.out.fasta_gz,
        params.use_phased
    )

    ch_logs = ch_logs.mix(
        DVPOLISH_CREATE_FINALASM.out.dvpolish_selection_tsv
        )
        .map { _meta, file -> file }

    ch_versions = ch_versions.mix(
        SAMTOOLS_FAIDX.out.versions.first(),
        DVPOLISH_CHUNKFA.out.versions.first(),
        DVPOLISH_PBMM2_INDEX.out.versions.first(),
        DVPOLISH_PBMM2_ALIGN.out.versions.first(),
        SAMTOOLS_VIEW.out.versions.first(),
        SAMTOOLS_MERGE.out.versions.first(),
        SAMTOOLS_INDEX.out.versions.first(),
        DEEPVARIANT.out.versions.first(),
        BCFTOOLS_VIEW.out.versions.first(),
        TABIX_TABIX.out.versions.first(),
        BCFTOOLS_MERGE.out.versions.first(),
        TABIX_TABIX_MERGED.out.versions.first(),
        BCFTOOLS_CONSENSUS.out.versions.first(),
        MERQURY_INPUT_ASM.out.versions.first(),
        MERQURY_POLISHED_ASM.out.versions.first(),
        DVPOLISH_CREATE_FINALASM.out.versions.first()
    )

    emit:
    assemblies = ch_polished_assemblies
    logs       = ch_logs
    versions   = ch_versions
}
