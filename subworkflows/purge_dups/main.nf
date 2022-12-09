/*
 * Workflow based around the purge_dups tool
 * https://github.com/dfguan/purge_dups
 */

// TODO:: purgedups are no on nf-core
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_READS    } from "$projectDir/modules/nf-core/minimap2/align/main"
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ASSEMBLY } from "$projectDir/modules/nf-core/minimap2/align/main"
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
        .flatMap { meta, reads, assembly -> reads instanceof List ? reads.collect{ [ meta + [ single_end: true, build: assembly.id ], it, assembly.pri_fasta ] } : [ [ meta + [ single_end: true, build: assembly.id ], reads, assembly.pri_fasta ] ] }
        .multiMap { meta, reads, assembly -> 
            reads_ch: [ meta, reads ]
            assembly_ch: assembly
        }
        .set { input }
    reads_plus_assembly_ch
        .map { meta, reads, assembly -> [ meta + [ build: assembly.id ], assembly.pri_fasta ] }
        .set { assembly_ch }
    // Map pacbio reads
    MINIMAP2_ALIGN_READS(
        input.reads_ch,
        input.assembly_ch,
        false, // bam output
        false, // cigar in paf file
        false  // cigar in bam file
    )
    PURGEDUPS_PBCSTAT( MINIMAP2_ALIGN_READS.out.paf.groupTuple() )
    PURGEDUPS_CALCUTS( PURGEDUPS_PBCSTAT.out.stat )
    // TODO:: Cutoffs can likely be estimated from genescope model output.

    // Split assembly and do self alignment
    PURGEDUPS_SPLITFA( assembly_ch )
    MINIMAP2_ALIGN_ASSEMBLY (
        PURGEDUPS_SPLITFA.out.split_fasta,
        [],    // Trigger read to read alignment
        false, // bam output
        false, // cigar in paf file
        false  // cigar in bam file
    )

    // TODO: Check the output from here
    PURGEDUPS_PURGEDUPS(
        PURGEDUPS_PBCSTAT.out.basecov
            .join( PURGEDUPS_CALCUTS.out.cutoff )
            .map { meta, cov, cutoff -> [ meta.findAll { !(it.key in [ 'single_end' ]) }, cov, cutoff ] }
            .join( MINIMAP2_ALIGN_ASSEMBLY.out.paf )
    ) 

    PURGEDUPS_GETSEQS( assembly_ch.join( PURGEDUPS_PURGEDUPS.out.bed ) )

    // TODO: Mix haplotigs back into haplotig set / Verify alternate contigs.

    emit:
    assembly = PURGEDUPS_GETSEQS.out.purged
    coverage = PURGEDUPS_PBCSTAT.out.basecov
}