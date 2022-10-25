/*
 * Workflow based around the purge_dups tool
 * https://github.com/dfguan/purge_dups
 */

include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_READS    } from "$projectDir/modules/nf-core/modules/minimap2/align/main"
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ASSEMBLY } from "$projectDir/modules/nf-core/modules/minimap2/align/main"
include { PURGEDUPS_CALCUTS                         } from "$projectDir/modules/local/purgedups/calcuts"
include { PURGEDUPS_GETSEQS                         } from "$projectDir/modules/local/purgedups/getseqs"
include { PURGEDUPS_PBCSTAT                         } from "$projectDir/modules/local/purgedups/pbcstat"
include { PURGEDUPS_SPLITFA                         } from "$projectDir/modules/local/purgedups/splitfa"
include { PURGEDUPS_PURGEDUPS                       } from "$projectDir/modules/local/purgedups/purgedups"

workflow PURGE_DUPLICATES {

    take:
    reads_plus_assembly_ch     // [ meta, [reads], [assembly] ], where reads are the pacbio files, and assembly is the primary and alternate asms

    main:
    reads_plus_assembly_ch
        .flatMap { meta, reads, assembly -> reads instanceof List ? reads.collect{ [ meta + [ single_end: true ], it, assembly.pri_fasta ] } : [ [ meta + [ single_end: true ], reads, assembly.pri_fasta ] ] }
        .multiMap { meta, reads, assembly -> 
            reads_ch: [ meta, reads ]
            assembly_ch: assembly
        }
        .set { input }
    reads_plus_assembly_ch
        .map { meta, reads, assembly -> [ meta, assembly.pri_fasta ] }
        .set { assembly_ch }
    /*
    # Map Pacbio CSS reads
    for i in $pb_list
    do
        minimap2 -xasm20 $pri_asm $i | gzip -c - > $i.paf.gz
    done
    bin/pbcstat *.paf.gz (produces PB.base.cov and PB.stat files)
    bin/calcuts PB.stat > cutoffs 2>calcults.log

    # Split assembly and do self alignment
    bin/split_fa $pri_asm > $pri_asm.split
    minimap2 -xasm5 -DP $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz

    # Purge haplotigs
    bin/purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log

    # purged primary and haplotig sequences from draft assembly
    bin/get_seqs -e dups.bed $pri_asm 
    */
    // Map pacbio reads
    MINIMAP2_ALIGN_READS(
        input.reads_ch,
        input.assembly_ch,
        false, // bam output
        false, // cigar in paf file
        false  // cigar in bam file
    )
    PURGEDUPS_PBCSTAT( MINIMAP2_ALIGN_READS.out.paf.collect() )
    PURGEDUPS_CALCUTS( PURGEDUPS_PBCSTAT.out.stat )

    // Split assembly and do self alignment
    PURGEDUPS_SPLITFA( assembly_ch )
    MINIMAP2_ALIGN_ASSEMBLY (
        PURGEDUPS_SPLITFA.out.split_fasta,
        [],    // Trigger read to read alignment
        false, // bam output
        false, // cigar in paf file
        false  // cigar in bam file
    )

    PURGEDUPS_PURGEDUPS(
        PURGEDUPS_PBCSTAT.out.basecov
            .join( PURGEDUPS_CALCUTS.out.cutoff )
            .join( MINIMAP2_ALIGN_ASSEMBLY.out.paf )
    )

    PURGEDUPS_GETSEQS( assembly_ch.join( PURGEDUPS_PURGEDUPS.out.bed ) )

    emit:
    assembly = PURGEDUPS_GETSEQS.out.assembly
    coverage = PURGEDUPS_PBCSTAT.out.basecov
}