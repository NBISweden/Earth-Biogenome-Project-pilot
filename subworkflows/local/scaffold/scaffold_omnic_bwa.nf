include { constructAssemblyRecord                 } from "../../../modules/local/functions"
include { getPrimaryAssembly                      } from "../../../modules/local/functions"
include { joinByMetaKeys                          } from "../../../modules/local/functions"
include { combineByMetaKeys                       } from "../../../modules/local/functions"
include { BWAMEM2_INDEX                           } from "../../../modules/nf-core/bwamem2/index/main"
include { BWAMEM2_MEM                             } from "../../../modules/nf-core/bwamem2/mem/main"
include { SAMTOOLS_FAIDX                          } from "../../../modules/nf-core/samtools/faidx/main"
include { PAIRTOOLS_PARSE                         } from "../../../modules/nf-core/pairtools/parse/main"
include { PAIRTOOLS_SORT                          } from "../../../modules/nf-core/pairtools/sort/main"
include { PAIRTOOLS_MERGE                         } from "../../../modules/nf-core/pairtools/merge/main"
include { PAIRTOOLS_DEDUP                         } from "../../../modules/nf-core/pairtools/dedup/main"
include { PAIRTOOLS_SPLIT                         } from "../../../modules/nf-core/pairtools/split/main"
include { YAHS                                    } from "../../../modules/nf-core/yahs/main.nf"


workflow SCAFFOLD_OMNIC_BWA {
    take:
    ch_assemblies // [ meta, assembly ]
    ch_hic        // [ meta, hic-pairs ]

    main:

    ch_versions = Channel.empty()

    BWAMEM2_INDEX(getPrimaryAssembly(ch_assemblies))
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())


    SAMTOOLS_FAIDX(
        getPrimaryAssembly(ch_assemblies),
        [[], []],
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    combineByMetaKeys(
        combineByMetaKeys(
            ch_hic,
            BWAMEM2_INDEX.out.index,
            keySet: ['id', 'sample'],
            meta: 'rhs',
        ),
        getPrimaryAssembly(ch_assemblies),
        keySet: ['id', 'sample'],
        meta: 'rhs',
    ).multiMap { meta, hic_reads, index, fasta ->
        reads: [meta, hic_reads]
        index: [meta, index]
        fasta: [meta, fasta]
    }.set { bwamem2_input }

    BWAMEM2_MEM(
        bwamem2_input.reads,
        bwamem2_input.index,
        bwamem2_input.fasta,
        false,
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    SAMTOOLS_FAIDX.out.fai
        .map { meta, fai ->
            fai
                .splitCsv(sep: '\t', header: false)
                .collect { row -> row[0..1].join('\t') }
                .join('\n')
        }
        .collectFile()
        .set { chrom_sizes }

    // Combine hi-c alignment with chrom sizes
    BWAMEM2_MEM.out.bam
        .combine(chrom_sizes)
        .multiMap { meta, hic_bam, chr_lengths ->
            hicbams: [meta, hic_bam]
            lengths: chr_lengths
        }
        .set { pairtools_parse_input }

    PAIRTOOLS_PARSE(
        pairtools_parse_input.hicbams,
        pairtools_parse_input.lengths,
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_PARSE.out.versions.first())

    PAIRTOOLS_SORT(PAIRTOOLS_PARSE.out.pairsam)
    ch_versions = ch_versions.mix(PAIRTOOLS_SORT.out.versions.first())

    PAIRTOOLS_SORT.out.sorted
        .groupTuple()
        .branch { meta, pairsam ->
            single: pairsam.size() == 1
            return [meta, pairsam.first()]
            multi: true
            return [meta, pairsam]
        }
        .set { pairtools_sorted }

    PAIRTOOLS_MERGE(pairtools_sorted.multi)
    ch_versions = ch_versions.mix(PAIRTOOLS_MERGE.out.versions.first())

    PAIRTOOLS_DEDUP(
        PAIRTOOLS_MERGE.out.pairs.mix(
            pairtools_sorted.single
        )
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_DEDUP.out.versions.first())

    joinByMetaKeys(
        PAIRTOOLS_DEDUP.out.pairs,
        getPrimaryAssembly(ch_assemblies),
        keySet: ['id', 'sample'],
        meta: 'rhs',
    ).multiMap { meta, pairs, fasta ->
        pairs: [meta, pairs]
        fasta: [meta, fasta]
    }.set { pairtools_split_input }

    PAIRTOOLS_SPLIT(
        pairtools_split_input.pairs,
        pairtools_split_input.fasta,
        true,
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_SPLIT.out.versions.first())

    combineByMetaKeys(
        combineByMetaKeys(
            PAIRTOOLS_SPLIT.out.bam,
            getPrimaryAssembly(ch_assemblies),
            keySet: ['id', 'sample'],
            meta: 'rhs',
        ),
        SAMTOOLS_FAIDX.out.fai,
        keySet: ['id', 'sample'],
        meta: 'rhs',
    ).multiMap { meta, bam, fasta, fai ->
        bam: [meta, bam]
        fasta: [fasta]
        fai: [fai]
    }.set { yahs_input }

    YAHS(
        yahs_input.bam,
        yahs_input.fasta,
        yahs_input.fai,
    )
    ch_versions = ch_versions.mix(YAHS.out.versions.first())

    // Consensus case:
    // Preserve haplotigs from purge dups
    ch_scaff_and_alt = ch_assemblies
        .map { meta, assembly -> [meta, assembly.alt_fasta] }
        .mix(YAHS.out.scaffolds_fasta)
    ch_scaffolded_assemblies = constructAssemblyRecord(ch_scaff_and_alt, false)

    PAIRTOOLS_PARSE.out.stat
        .mix(
            PAIRTOOLS_DEDUP.out.stat
        )
        .map { meta, stats -> stats }
        .set { logs }

    emit:
    assemblies = ch_scaffolded_assemblies
    logs
    versions   = ch_versions
}
