include { constructAssemblyRecord                 } from "../../../modules/local/functions"
include { getEachAssembly                         } from "../../../modules/local/functions"
include { getPrimaryAssembly                      } from "../../../modules/local/functions"
include { joinByMetaKeys                          } from "../../../modules/local/functions"
include { combineByMetaKeys                       } from "../../../modules/local/functions"
include { BWAMEM2_INDEX as BWAMEM2_INDEX_SCAFFOLD } from "../../../modules/nf-core/bwamem2/index/main"
include { BWAMEM2_MEM as BWAMEM2_MEM_SCAFFOLD     } from "../../../modules/nf-core/bwamem2/mem/main"
include { SAMTOOLS_FAIDX                          } from "../../../modules/nf-core/samtools/faidx/main"
include { PAIRTOOLS                               } from "../../../modules/local/pairtools/main"
include { YAHS                                    } from "../../../modules/nf-core/yahs/main.nf"

/*
    PURPOSE: Performs scaffolding of haplotype/consensus assemblies with a scaffolding tool.

    SCIENTIFIC RATIONALE:
    Each haplotype assembly is scaffolded separately as that is the underlying assumption of Yahs.
    The assembly should only be merged for curation.
 */

workflow SCAFFOLD {
    take:
    ch_assemblies // [ meta, assembly ]
    ch_hic        // [ meta, hic-pairs ]

    main:
    ch_to_scaffold = params.use_phased ?
        getEachAssembly(ch_assemblies)
            .flatMap { meta, hap1_hap2 -> hap1_hap2.withIndex()
                .collect { hap_asm, idx -> tuple( meta + [ haplotype: "hap${idx+1}" ], hap_asm ) }
            }
        : getPrimaryAssembly(ch_assemblies)
    BWAMEM2_INDEX_SCAFFOLD ( ch_to_scaffold )

    SAMTOOLS_FAIDX (
        ch_to_scaffold,
        [ [ ], [ ] ]
    )

    combineByMetaKeys( // Combine (Hi-C + index) with (BWA_INDEX + Assembly)
        ch_hic,
        joinByMetaKeys( // Join BWA index with Assembly
            BWAMEM2_INDEX_SCAFFOLD.out.index,
            ch_to_scaffold,
            keySet: ['id','sample','assembly', 'haplotype'],
            meta: 'rhs'
        ),
        keySet: ['id','sample'],
        meta: 'rhs'
    ).multiMap{ meta, hic_reads, index, fasta ->
        reads: [ meta, hic_reads ]
        index: [ meta, index ]
        fasta: [ meta, fasta ]
    }.set{ bwamem2_input }

    BWAMEM2_MEM_SCAFFOLD (
        bwamem2_input.reads,
        bwamem2_input.index,
        bwamem2_input.fasta,
        false
    )

    SAMTOOLS_FAIDX.out.fai
        .map{ _meta, fai ->
            fai.splitCsv( sep: '\t', header: false )
                .collect{ row -> row[ 0..1 ].join('\t') }
                .join('\n')
        }.collectFile()
        .set { chrom_sizes }

    // For each assembly, group hi-c alignment bams and combine with chrom sizes
    joinByMetaKeys(
        BWAMEM2_MEM_SCAFFOLD.out.bam.groupTuple()
            .combine(chrom_sizes)
            .map{ meta, hic_bam, chr_lengths -> [ meta, hic_bam, chr_lengths ] },
        ch_to_scaffold,
        keySet: ['id','sample','assembly','haplotype'],
        meta: 'rhs'
    ).multiMap{ meta, hic_bam, chr_lengths, fasta ->
            hicbams: [ meta, hic_bam ]
            lengths: chr_lengths
            fasta:  fasta
    }.set{ pairtools_input }

    PAIRTOOLS (
        pairtools_input.hicbams,
        pairtools_input.lengths,
        pairtools_input.fasta,
        true    // sort_bam
    )

    combineByMetaKeys( // Combine bam with (assembly + fai)
        PAIRTOOLS.out.bam, // FIXME: note this expects BAM, while module can output CRAM
        joinByMetaKeys( // join assembly + fai
            ch_to_scaffold,
            SAMTOOLS_FAIDX.out.fai,
            keySet: ['id','sample','assembly','haplotype'],
            meta: 'rhs'
        ),
        keySet: ['id','sample','assembly'],
        meta: 'rhs'
    ).multiMap{ meta, bam, fasta, fai ->
        bam:   [ meta, bam ]
        fasta: [ fasta ]
        fai:   [ fai ]
    }.set{ yahs_input }

    YAHS(
        yahs_input.bam,
        yahs_input.fasta,
        yahs_input.fai
    )

    // // Consensus case:
    // // Preserve haplotigs from purge dups
    // ch_scaff_and_alt = ch_assemblies
    //     .filter { _meta, assembly -> assembly.alt_fasta }
    //     .map { meta, assembly -> [ meta, assembly.alt_fasta ] }
    //     .mix( YAHS.out.scaffolds_fasta )
    ch_scaffolded_assemblies = constructAssemblyRecord(
        YAHS.out.scaffolds_fasta
            .map{ meta, fasta -> tuple( meta.subMap( meta.keySet() - ['haplotype'] ), fasta) },
        params.use_phased
    )

    logs = PAIRTOOLS.out.stat
        .flatten()

    versions = BWAMEM2_INDEX_SCAFFOLD.out.versions.first().mix(
        SAMTOOLS_FAIDX.out.versions.first(),
        BWAMEM2_MEM_SCAFFOLD.out.versions.first(),
        PAIRTOOLS.out.versions.first(),
        YAHS.out.versions.first()
    )

    emit:
    assemblies = ch_scaffolded_assemblies
    logs
    versions
}