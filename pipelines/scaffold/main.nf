#!/usr/bin/env nextflow

/*
 * Scaffolding of haplotype/consensus assemblies with YaHS and Hi-C reads, followed by assembly QC.
 *
 * Each haplotype assembly is scaffolded separately as that is the underlying assumption of YaHS.
 * The haplotypes are only merged for curation.
 *
 * Takes one wide assembly record per assembly build (see types/records.nf). The Hi-C read pairs
 * are read from `hic`. Returns the scaffolded assemblies with the same fields.
 */

nextflow.enable.types = true

include { BWAMEM2_INDEX as BWAMEM2_INDEX_SCAFFOLD } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as BWAMEM2_MEM_SCAFFOLD     } from '../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_FAIDX                          } from '../../modules/nf-core/samtools/faidx/main'
include { EXTRACT_CHROMSIZES                      } from '../../modules/local/extract/chromsizes/main'
include { PAIRTOOLS                               } from '../../modules/local/pairtools/main'
include { YAHS                                    } from '../../modules/nf-core/yahs/main'
include { SEQKIT_MERGE as MERGE_HAPLOTYPES        } from '../../modules/local/seqkit/merge/main'
include { workflow as ASSEMBLY_QC_SCAFFOLDED      } from '../assembly_qc'
include { Assembly; KmerParams; BuscoParams; assemblyWith } from '../../types/records'

params {
    assemblies: Channel<Assembly>  // Assemblies to scaffold
    provided  : Channel<Assembly>? // Assemblies of this stage that were supplied. Evaluated, not scaffolded.

    // Scaffolding
    use_phased     : Boolean = false      // Scaffold the primary and alternate assembly separately
    hic_type       : String  = 'arima_v2' // Hi-C library type: arima_v2, arima_v1, omni-c or ''
    hic_bam_sorting: String  = 'by_name'  // Sort of the Hi-C bam file: by_coordinate or by_name

    // Assembly QC
    busco     : BuscoParams = record(lineages: null, lineages_db_path: null, max_lineages: null)
    fastk     : KmerParams  = record(kmer_size: null)
    meryl     : KmerParams  = record(kmer_size: null)
    reference : Path?
}

// Identifies the sample, assembly build and haplotype (if any)
def haplotypeKey(meta: Record) -> String {
    return "${meta.id}:${meta.assembly.build}:${meta.haplotype ?: ''}"
}

// File name prefix for the scaffolding of one haplotype
def scaffoldPrefix(meta: Record) -> String {
    return [ meta.id, meta.assembly.build, meta.haplotype ?: '' ].findAll { p -> p }.join('_')
}

workflow {

    main:
    // Scaffold each haplotype separately when phased
    to_scaffold = params.assemblies.flatMap { a ->
        def hic = a.hic ?: []
        def metas = params.use_phased
            ? [ a.meta + record(haplotype: 'hap1'), a.meta + record(haplotype: 'hap2') ]
            : [ a.meta ]
        def fastas = params.use_phased ? [ a.pri_fasta, a.alt_fasta ] : [ a.pri_fasta ]
        (0..<metas.size()).collect { i ->
            record(id: haplotypeKey(metas[i]), meta: metas[i], fasta: fastas[i], hic: hic)
        }.toList()
    }

    index = BWAMEM2_INDEX_SCAFFOLD(
        to_scaffold.map { h -> record(id: h.id, meta: h.meta, fasta: h.fasta) }
    )

    faidx = SAMTOOLS_FAIDX(
        to_scaffold.map { h -> record(id: h.id, meta: h.meta, fasta: h.fasta, fai: null, get_sizes: false) }
    )

    // Align each Hi-C read pair to the BWA index of the haplotype
    aligned = BWAMEM2_MEM_SCAFFOLD(
        to_scaffold
            .join(index.map { i -> record(id: i.id, index: i.index) }, by: 'id')
            .flatMap { x ->
                x.hic.collect { pair ->
                    record(
                        id: x.id,
                        meta: x.meta,
                        reads: pair,
                        index: x.index,
                        fasta: x.fasta,
                        sort_bam: false,
                        bwamem2_args: '-SP -T0',
                        bwamem2_args2: '',
                        output_prefix: "${scaffoldPrefix(x.meta)}_${pair[0].name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')}"
                    )
                }.toList()
            }
    )

    chromsizes = EXTRACT_CHROMSIZES(
        faidx.map { f -> record(id: f.id, meta: f.meta, faidx: f.fai) }
    )

    // For each assembly, group hi-c alignment bams and combine with chrom sizes
    pairs = PAIRTOOLS(
        aligned
            .map { a -> tuple(haplotypeKey(a.meta), a) }
            .groupBy()
            .map { key, alignments ->
                record(id: key, meta: alignments.toList()[0].meta, bam: alignments.collect { a -> a.bam }.toList().toSorted())
            }
            .join(chromsizes.map { c -> record(id: c.id, chromsizes: c.chr_sizes) }, by: 'id')
            .join(to_scaffold.map { h -> record(id: h.id, fasta: h.fasta) }, by: 'id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    bam: x.bam,
                    chromsizes: x.chromsizes,
                    fasta: x.fasta,
                    sort_bam: true,
                    pairtools_args: "--min-mapq 20 --walks-policy 5unique --max-inter-align-gap 30",
                    pairtools_args2: '--tmpdir=\$PWD',
                    pairtools_args3: null,
                    pairtools_args4: '--mark-dups',
                    pairtools_args5: null,
                    pairtools_args6: params.hic_bam_sorting == "by_coordinate" ? "--write-index" : "-n",
                    output_prefix: scaffoldPrefix(x.meta)
                )
            }
    )

    // FIXME: note this expects BAM, while the module can output CRAM
    scaffolds = YAHS(
        pairs
            .map { p -> record(id: p.id, meta: p.meta, hic_map: p.bam) }
            .join(to_scaffold.map { h -> record(id: h.id, fasta: h.fasta) }, by: 'id')
            .join(faidx.map { f -> record(id: f.id, fai: f.fai) }, by: 'id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    fasta: x.fasta,
                    fai: x.fai,
                    hic_map: x.hic_map,
                    agp: null,
                    yahs_args: params.hic_type == "arima_v2" ? "-e GATC,GANTC,CTNAG,TTAA" : params.hic_type == "arima_v1" ? "-e GATC,GANTC" : "",
                    output_prefix: scaffoldPrefix(x.meta)
                )
            }
    )

    // Recombine the haplotypes of each assembly, primary first when phased
    scaffolded = scaffolds
        .filter { s -> s.scaffolds_fasta != null }
        .map { s -> tuple("${s.meta.id}:${s.meta.assembly.build}", s.scaffolds_fasta) }
        .groupBy()
        .map { key, fastas ->
            def sorted = fastas.toSorted { a, b -> params.use_phased ? a.name <=> b.name : b.name <=> a.name }
            record(assembly_key: key, fastas: sorted)
        }
        .join(params.assemblies.map { a -> record(assembly_key: "${a.meta.id}:${a.meta.assembly.build}", assembly: a) }, by: 'assembly_key')

    // The wide record, with the scaffolded assembly files
    scaffolded_assemblies = scaffolded.map { x ->
        assemblyWith(
            x.assembly,
            [
                pri_fasta: x.fastas[0],
                alt_fasta: x.fastas.size() > 1 ? x.fastas[1] : null,
                pri_gfa: null,
                alt_gfa: null
            ]
        )
    }
    assemblies = scaffolded_assemblies.mix(params.provided ?: channel.empty())

    merged_haplotypes = MERGE_HAPLOTYPES(
        scaffolded
            .filter { x -> x.fastas.size() > 1 }
            .map { x ->
                record(
                    id: x.assembly.id,
                    meta: x.assembly.meta,
                    fastx: x.fastas.toSorted { a, b -> a.name <=> b.name },
                    out_ext: "fasta.gz",
                    seqkit_args: "-w 100",
                    output_prefix: "${x.assembly.meta.id}_${x.assembly.meta.assembly.build}"
                )
            }
    )

    qc = ASSEMBLY_QC_SCAFFOLDED(
        record(
            assemblies: assemblies,
            busco: params.busco,
            fastk: params.fastk,
            meryl: params.meryl,
            reference: params.reference,
            use_phased: params.use_phased
        )
    )

    publish:
    assemblies        = assemblies
    pairtools         = pairs
    yahs              = scaffolds
    merged_haplotypes = merged_haplotypes
    merquryfk         = qc.merquryfk
    merqury           = qc.merqury
    busco             = qc.busco
    gfastats          = qc.gfastats
    quast             = qc.quast
    logs              = pairs.flatMap { p -> p.stat.toList() }.mix(qc.logs)
}

output {
    assemblies: Channel<Assembly> {
        path { a ->
            a.pri_fasta >> ("assemblies/${a.id}/")
            a.alt_fasta >> ("assemblies/${a.id}/")
        }
    }
    pairtools: Channel<Record> {
        path { r -> r.stat >> ("pairtools/${r.meta.assembly.build}/") }
    }
    yahs: Channel<Record> {
        path { r -> "yahs/${r.meta.assembly.build}" }
    }
    merged_haplotypes: Channel<Record> {
        path { r -> "merged_haps/${r.meta.assembly.build}" }
    }
    merquryfk: Channel<Record> {
        path { r -> "merquryfk/${r.meta.assembly.build}" }
    }
    merqury: Channel<Record> {
        path { r -> "merqury/${r.meta.assembly.build}" }
    }
    busco: Channel<Record> {
        path { r ->
            r.short_summaries_txt >> ("busco/${r.meta.assembly.build}/")
            r.short_summaries_json >> ("busco/${r.meta.assembly.build}/")
        }
    }
    gfastats: Channel<Record> {
        path { r -> "gfastats/${r.meta.assembly.build}" }
    }
    quast: Channel<Record> {
        path 'quast'
    }
    logs: Channel<Path> {
        path 'logs'
    }
}
