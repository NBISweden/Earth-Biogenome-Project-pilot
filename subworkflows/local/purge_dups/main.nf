/*
 * Workflow based around the purge_dups tool
 * https://github.com/dfguan/purge_dups
 */

include { joinByMetaKeys                                       } from "../../../modules/local/functions"
include { combineByMetaKeys                                    } from "../../../modules/local/functions"
include { constructAssemblyRecord                              } from "../../../modules/local/functions"
include { getPrimaryAssembly                                   } from "../../../modules/local/functions"
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_READS               } from "../../../modules/nf-core/minimap2/align/main"
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ASSEMBLY_PRIMARY    } from "../../../modules/nf-core/minimap2/align/main"
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ASSEMBLY_ALTERNATE  } from "../../../modules/nf-core/minimap2/align/main"
include { PURGEDUPS_PBCSTAT                                    } from "../../../modules/nf-core/purgedups/pbcstat"
include { PURGEDUPS_CALCUTS                                    } from "../../../modules/nf-core/purgedups/calcuts"
include { PURGEDUPS_HISTPLOT                                   } from "../../../modules/nf-core/purgedups/histplot"
include { PURGEDUPS_SPLITFA as PURGEDUPS_SPLITFA_PRIMARY       } from "../../../modules/nf-core/purgedups/splitfa"
include { PURGEDUPS_SPLITFA as PURGEDUPS_SPLITFA_ALTERNATE     } from "../../../modules/nf-core/purgedups/splitfa"
include { PURGEDUPS_PURGEDUPS as PURGEDUPS_PURGEDUPS_PRIMARY   } from "../../../modules/nf-core/purgedups/purgedups"
include { PURGEDUPS_PURGEDUPS as PURGEDUPS_PURGEDUPS_ALTERNATE } from "../../../modules/nf-core/purgedups/purgedups"
include { PURGEDUPS_GETSEQS as PURGEDUPS_GETSEQS_PRIMARY       } from "../../../modules/nf-core/purgedups/getseqs"
include { PURGEDUPS_GETSEQS as PURGEDUPS_GETSEQS_ALTERNATE     } from "../../../modules/nf-core/purgedups/getseqs"
include { SEQKIT_SEQ                                           } from "../../../modules/nf-core/seqkit/seq/main"

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
    PURGEDUPS_SPLITFA_PRIMARY( getPrimaryAssembly( ch_assemblies ) )
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
            keySet: [ 'sample', 'assembly' ],
            meta: 'rhs'
        )
    )
    PURGEDUPS_GETSEQS_PRIMARY(
        PURGEDUPS_SPLITFA_PRIMARY.out.merged_fasta
            .join( PURGEDUPS_PURGEDUPS_PRIMARY.out.bed )
    )
    def ch_to_format = PURGEDUPS_GETSEQS_PRIMARY.out.purged
        .mix( PURGEDUPS_GETSEQS_PRIMARY.out.haplotigs )

    if( params.use_phased ){
        // Purge alternate contigs.
        reads_plus_assembly_ch
            .filter { _meta, _reads, assembly -> assembly.alt_fasta != null }
            .map { meta, _reads, assembly -> [ meta, assembly.alt_fasta ] }
            //! WARN Purges only primary haplotigs when using consensus
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
                keySet: [ 'sample', 'assembly' ],
                meta: 'rhs'
            )
        )
        PURGEDUPS_GETSEQS_ALTERNATE(
            PURGEDUPS_SPLITFA_ALTERNATE.out.merged_fasta
                .join( PURGEDUPS_PURGEDUPS_ALTERNATE.out.bed )
        )
        ch_to_format = PURGEDUPS_GETSEQS_PRIMARY.out.purged.mix( PURGEDUPS_GETSEQS_ALTERNATE.out.purged )
    }
    // Format sequences and enforce line breaks
    SEQKIT_SEQ( ch_to_format )
    // If phased [ meta, [*hap.fa, *purged.fa] ] // don't sort by name
    // else [ meta, [*hap1.fa, *hap2.fa] ]       // sort by name
    ch_purged_assemblies = constructAssemblyRecord( SEQKIT_SEQ.out.fastx, params.use_phased )

    PURGEDUPS_HISTPLOT.out.png
        .mix( PURGEDUPS_PURGEDUPS_PRIMARY.out.bed )
        .map { _meta, file -> file }
        .set { logs }

    MINIMAP2_ALIGN_READS.out.versions.first().mix(
        PURGEDUPS_PBCSTAT.out.versions.first(),
        PURGEDUPS_CALCUTS.out.versions.first(),
        PURGEDUPS_HISTPLOT.out.versions.first(),
        PURGEDUPS_SPLITFA_PRIMARY.out.versions.first(),
        MINIMAP2_ALIGN_ASSEMBLY_PRIMARY.out.versions.first(),
        PURGEDUPS_PURGEDUPS_PRIMARY.out.versions.first(),
        PURGEDUPS_GETSEQS_PRIMARY.out.versions.first(),
        SEQKIT_SEQ.out.versions.first(),
    ).set { versions }

    emit:
    assemblies = ch_purged_assemblies
    coverage   = PURGEDUPS_PBCSTAT.out.basecov
    logs
    versions
}