/*
 * Workflow based around the purge_dups tool
 * https://github.com/dfguan/purge_dups
 */

include { joinByMetaKeys                                       } from "$projectDir/modules/local/functions"
include { combineByMetaKeys                                    } from "$projectDir/modules/local/functions"
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_READS               } from "$projectDir/modules/nf-core/minimap2/align/main"
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ASSEMBLY_PRIMARY    } from "$projectDir/modules/nf-core/minimap2/align/main"
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ASSEMBLY_ALTERNATE  } from "$projectDir/modules/nf-core/minimap2/align/main"
include { PURGEDUPS_PBCSTAT                                    } from "$projectDir/modules/nf-core/purgedups/pbcstat"
include { PURGEDUPS_CALCUTS                                    } from "$projectDir/modules/nf-core/purgedups/calcuts"
include { PURGEDUPS_HISTPLOT                                   } from "$projectDir/modules/nf-core/purgedups/histplot"
include { PURGEDUPS_SPLITFA as PURGEDUPS_SPLITFA_PRIMARY       } from "$projectDir/modules/nf-core/purgedups/splitfa"
include { PURGEDUPS_SPLITFA as PURGEDUPS_SPLITFA_ALTERNATE     } from "$projectDir/modules/nf-core/purgedups/splitfa"
include { PURGEDUPS_PURGEDUPS as PURGEDUPS_PURGEDUPS_PRIMARY   } from "$projectDir/modules/nf-core/purgedups/purgedups"
include { PURGEDUPS_PURGEDUPS as PURGEDUPS_PURGEDUPS_ALTERNATE } from "$projectDir/modules/nf-core/purgedups/purgedups"
include { PURGEDUPS_GETSEQS as PURGEDUPS_GETSEQS_PRIMARY       } from "$projectDir/modules/nf-core/purgedups/getseqs"
include { PURGEDUPS_GETSEQS as PURGEDUPS_GETSEQS_ALTERNATE     } from "$projectDir/modules/nf-core/purgedups/getseqs"

workflow PURGE_DUPLICATES {

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
            assembly_ch: assembly
        }
        .set { input }
    reads_plus_assembly_ch
        .map { meta, reads, assembly -> [ meta, assembly.pri_fasta ] }
        .set { primary_assembly_ch }
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
    PURGEDUPS_HISTPLOT( PURGEDUPS_PBCSTAT.out.stat.join( PURGEDUPS_CALCUTS.out.cutoff ) )

    // Purge primary assembly
    PURGEDUPS_SPLITFA_PRIMARY( primary_assembly_ch )
    MINIMAP2_ALIGN_ASSEMBLY_PRIMARY(
        PURGEDUPS_SPLITFA_PRIMARY.out.split_fasta,
        [],    // Trigger read to read alignment
        false, // bam output
        false, // cigar in paf file
        false  // cigar in bam file
    )
    PURGEDUPS_PURGEDUPS_PRIMARY(
        joinByMetaKeys(
            PURGEDUPS_PBCSTAT.out.basecov.join( PURGEDUPS_CALCUTS.out.cutoff ),
            MINIMAP2_ALIGN_ASSEMBLY_PRIMARY.out.paf,
            keySet: ['sample','assembly'],
            meta: 'rhs'
        )
    )
    PURGEDUPS_GETSEQS_PRIMARY(
        PURGEDUPS_SPLITFA_PRIMARY.out.merged_fasta
            .join( PURGEDUPS_PURGEDUPS_PRIMARY.out.bed )
    )

    // Purge alternate contigs. // TODO: Skip when using consensus
    reads_plus_assembly_ch
        .filter { meta, reads, assembly -> assembly.alt_fasta != null }
        .map { meta, reads, assembly -> [ meta, assembly.alt_fasta ] }
        // Purges only primary haplotigs when using consensus
        .mix( PURGEDUPS_GETSEQS_PRIMARY.out.haplotigs )
        .groupTuple()  // TODO Find size to prevent blocking
        .set { alternate_assembly_ch }
    PURGEDUPS_SPLITFA_ALTERNATE( alternate_assembly_ch )
    MINIMAP2_ALIGN_ASSEMBLY_ALTERNATE(
        PURGEDUPS_SPLITFA_ALTERNATE.out.split_fasta,
        [],    // Trigger read to read alignment
        false, // bam output
        false, // cigar in paf file
        false  // cigar in bam file
    )
    PURGEDUPS_PURGEDUPS_ALTERNATE(
        joinByMetaKeys(
            PURGEDUPS_PBCSTAT.out.basecov.join( PURGEDUPS_CALCUTS.out.cutoff ),
            MINIMAP2_ALIGN_ASSEMBLY_ALTERNATE.out.paf,
            keySet: ['sample','assembly'],
            meta: 'rhs'
        )
    )
    PURGEDUPS_GETSEQS_ALTERNATE(
        PURGEDUPS_SPLITFA_ALTERNATE.out.merged_fasta
            .join( PURGEDUPS_PURGEDUPS_ALTERNATE.out.bed )
    )
    ch_purged_assemblies = PURGEDUPS_GETSEQS_PRIMARY.out.purged
        .mix(PURGEDUPS_GETSEQS_ALTERNATE.out.purged)
        .groupTuple( sort: { a, b -> a.name <=> b.name } )
        .map { meta, fasta ->
            def asm_meta = meta.assembly.subMap(['assembler','stage','id','build'])
            [ meta, asm_meta + (fasta.size() == 1 ? [ pri_fasta: fasta[0] ] : [ pri_fasta: fasta[0], alt_fasta: fasta[1] ] ) ]
        }

    emit:
    assemblies = ch_purged_assemblies
    coverage   = PURGEDUPS_PBCSTAT.out.basecov
}