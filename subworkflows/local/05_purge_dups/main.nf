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
    // Combine reads with assemblies and prepare for alignment
    ch_reads_with_assemblies = combineByMetaKeys(
        ch_hifi,
        ch_assemblies,
        keySet: ['id', 'sample'],
        meta: 'rhs'
    )

    ch_reads_and_assembly = ch_reads_with_assemblies
        // Flatten read sets and add single_end flag for minimap2
        .flatMap { meta, reads, assembly ->
            def read_list = reads instanceof List ? reads : [reads]
            if( params.use_phased ){
                read_list.collect { read_set ->
                    [meta + [single_end: true, haplotype: 'hap1'], read_set, assembly.pri_fasta]
                }
                + read_list.collect { read_set ->
                    [meta + [single_end: true, haplotype: 'hap2'], read_set, assembly.alt_fasta]
                }
            } else {
                read_list.collect { read_set ->
                    [meta + [single_end: true], read_set, assembly.pri_fasta]
                }
            }
        }
        .multiMap { meta, reads, assembly ->
            reads: [meta, reads]
            assembly: assembly
            meta_assembly : [meta, assembly]
        }

    MINIMAP2_ALIGN_READS(
        ch_reads_and_assembly.reads,
        ch_reads_and_assembly.assembly,
        false, // bam output
        false, // cigar in paf
        false  // cigar in bam
    )

    PURGEDUPS_PBCSTAT(MINIMAP2_ALIGN_READS.out.paf.groupTuple())
    PURGEDUPS_CALCUTS(PURGEDUPS_PBCSTAT.out.stat)
    PURGEDUPS_HISTPLOT(
        PURGEDUPS_PBCSTAT.out.stat.join(PURGEDUPS_CALCUTS.out.cutoff)
    )
    ch_coverage_and_cutoffs = PURGEDUPS_PBCSTAT.out.basecov
        .join(PURGEDUPS_CALCUTS.out.cutoff)

    PURGEDUPS_SPLITFA_PRIMARY(
        ch_reads_and_assembly.meta_assembly
            .filter { meta, _paf -> meta.haplotype ? meta.haplotype == 'hap1' : true }
    )

    MINIMAP2_ALIGN_ASSEMBLY_PRIMARY(
        PURGEDUPS_SPLITFA_PRIMARY.out.split_fasta,
        [],    // Empty reference triggers self-alignment
        false, // bam output
        false, // cigar in paf
        false  // cigar in bam
    )

    PURGEDUPS_PURGEDUPS_PRIMARY(
        joinByMetaKeys(
            ch_coverage_and_cutoffs,
            MINIMAP2_ALIGN_ASSEMBLY_PRIMARY.out.paf,
            keySet: ['sample', 'assembly', 'haplotype'],
            meta: 'rhs'
        )
    )

    PURGEDUPS_GETSEQS_PRIMARY(
        PURGEDUPS_SPLITFA_PRIMARY.out.merged_fasta
            .join(PURGEDUPS_PURGEDUPS_PRIMARY.out.bed)
    )

    ch_assemblies_to_format = PURGEDUPS_GETSEQS_PRIMARY.out.purged
        .mix(PURGEDUPS_GETSEQS_PRIMARY.out.haplotigs)

    if (params.use_phased) {
        // Combine alternate assembly with primary haplotigs for purging
        ch_alternate_with_haplotigs = ch_reads_and_assembly.meta_assembly
            .filter { meta, _assembly -> meta.haplotype == 'hap2' }
            .mix( PURGEDUPS_GETSEQS_PRIMARY.out.haplotigs
                .map{ meta, haptigs -> tuple( meta + [ haplotype: 'hap2'], haptigs ) }
            )
            .groupTuple()

        PURGEDUPS_SPLITFA_ALTERNATE(ch_alternate_with_haplotigs)

        MINIMAP2_ALIGN_ASSEMBLY_ALTERNATE(
            PURGEDUPS_SPLITFA_ALTERNATE.out.split_fasta,
            [],    // Empty reference triggers self-alignment
            false, // bam output
            false, // cigar in paf
            false  // cigar in bam
        )

        PURGEDUPS_PURGEDUPS_ALTERNATE(
            joinByMetaKeys(
                ch_coverage_and_cutoffs,
                MINIMAP2_ALIGN_ASSEMBLY_ALTERNATE.out.paf,
                keySet: ['sample', 'assembly', 'haplotype'],
                meta: 'rhs'
            )
        )

        PURGEDUPS_GETSEQS_ALTERNATE(
            PURGEDUPS_SPLITFA_ALTERNATE.out.merged_fasta
                .join(PURGEDUPS_PURGEDUPS_ALTERNATE.out.bed)
        )

        // Override output channel with both primary and alternate purged assemblies
        ch_assemblies_to_format = PURGEDUPS_GETSEQS_PRIMARY.out.purged
            .mix(PURGEDUPS_GETSEQS_ALTERNATE.out.purged)
    }

    // Standardize sequence formatting and line breaks
    SEQKIT_SEQ(ch_assemblies_to_format)

    // Construct final assembly record
    // Phased: [meta, [hap0.purged.fa, hap1.purged.fa]] - sort by name
    // Unphased: [meta, [*purged.fa, *hap.fa]] - reverse sort by name
    ch_purged_assemblies = constructAssemblyRecord(
        SEQKIT_SEQ.out.fastx,
        params.use_phased
    )

    logs = PURGEDUPS_HISTPLOT.out.png
        .mix(PURGEDUPS_PURGEDUPS_PRIMARY.out.bed)
        .map { _meta, file -> file }

    versions = MINIMAP2_ALIGN_READS.out.versions.first()
        .mix(
            PURGEDUPS_PBCSTAT.out.versions.first(),
            PURGEDUPS_CALCUTS.out.versions.first(),
            PURGEDUPS_HISTPLOT.out.versions.first(),
            PURGEDUPS_SPLITFA_PRIMARY.out.versions.first(),
            MINIMAP2_ALIGN_ASSEMBLY_PRIMARY.out.versions.first(),
            PURGEDUPS_PURGEDUPS_PRIMARY.out.versions.first(),
            PURGEDUPS_GETSEQS_PRIMARY.out.versions.first(),
            SEQKIT_SEQ.out.versions.first()
        )

    emit:
    assemblies = ch_purged_assemblies
    coverage   = PURGEDUPS_PBCSTAT.out.basecov
    logs
    versions
}