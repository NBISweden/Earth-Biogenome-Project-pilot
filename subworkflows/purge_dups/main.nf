/*
 * Workflow based around the purge_dups tool
 * https://github.com/dfguan/purge_dups
 */

include { MINIMAP2_ALIGN } from "$projectDir/modules/nf-core/modules/minimap2/align/main"

workflow PURGE_DUPS {

    take:
    reads_plus_assembly_ch     // [ meta, [reads], [assembly] ], where reads are the pacbio files, and assembly is the primary and alternate asms

    main:
    MINIMAP2_ALIGN(
        reads_plus_assembly_ch.flatMap{ meta, reads, assembly -> reads.collect{ [ meta, it, assembly.pri_asm ] } }
            .multiMap { meta, reads, reference -> 
                reads_ch: [meta, reads]
                ref_ch: reference
            }
        false, // bam output
        false, // cigar in paf file
        false  // cigar in bam file
    )
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

    emit:
    assembly
    coverage
}