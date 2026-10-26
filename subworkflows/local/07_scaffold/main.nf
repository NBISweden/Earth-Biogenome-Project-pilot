nextflow.enable.types = true

include { BWAMEM2_INDEX as BWAMEM2_INDEX_SCAFFOLD } from "../../../modules/nf-core/bwamem2/index/main"
include { BWAMEM2_MEM as BWAMEM2_MEM_SCAFFOLD     } from "../../../modules/nf-core/bwamem2/mem/main"
include { SAMTOOLS_FAIDX                          } from "../../../modules/nf-core/samtools/faidx/main"
include { EXTRACT_CHROMSIZES                      } from "../../../modules/local/extract/chromsizes/main"
include { PAIRTOOLS                               } from "../../../modules/local/pairtools/main"
include { YAHS                                    } from "../../../modules/nf-core/yahs/main"
include { SEQKIT_MERGE as MERGE_HAPLOTYPES        } from "../../../modules/local/seqkit/merge/main"

/*
    PURPOSE: Performs scaffolding of haplotype/consensus assemblies with a scaffolding tool.

    SCIENTIFIC RATIONALE:
    Each haplotype assembly is scaffolded separately as that is the underlying assumption of Yahs.
    The assembly should only be merged for curation.
 */

record Reads {
    id   : String
    meta : Record
    reads: List<Path>
}

record Assembly {
    id       : String  // meta.assembly.build
    meta     : Record  // sample meta data, including meta.assembly
    pri_fasta: Path
    alt_fasta: Path?
    pri_gfa  : Path?
    alt_gfa  : Path?
}

// Identifies the sample, assembly build and haplotype (if any)
def haplotypeKey(meta: Record) -> String {
    return "${meta.id}:${meta.assembly.build}:${meta.haplotype ?: ''}"
}

// File name prefix for the scaffolding of one haplotype
def scaffoldPrefix(meta: Record) -> String {
    return [ meta.id, meta.assembly.build, meta.haplotype ?: '' ].findAll { p -> p }.join('_')
}

def scaffoldInput(meta: Record, fasta: Path) -> Record {
    return record(id: haplotypeKey(meta), meta: meta, fasta: fasta)
}

workflow SCAFFOLD {
    take:
    assemblies  : Channel<Assembly>
    hic         : Channel<Reads> // Hi-C read pairs
    use_phased  : Boolean
    hic_type    : String
    hic_bam_sorting: String

    main:
    // Scaffold each haplotype separately when phased
    to_scaffold = assemblies.flatMap { a ->
        def inputs = use_phased
            ? [ scaffoldInput(a.meta + record(haplotype: 'hap1'), a.pri_fasta), scaffoldInput(a.meta + record(haplotype: 'hap2'), a.alt_fasta) ]
            : [ scaffoldInput(a.meta, a.pri_fasta) ]
        inputs.toList()
    }

    index = BWAMEM2_INDEX_SCAFFOLD(
        to_scaffold.map { h -> record(id: h.id, meta: h.meta, fasta: h.fasta) }
    )

    faidx = SAMTOOLS_FAIDX(
        to_scaffold.map { h -> record(id: h.id, meta: h.meta, fasta: h.fasta, fai: null, get_sizes: false) }
    )

    // Combine Hi-C reads with the BWA index and assembly
    aligned = BWAMEM2_MEM_SCAFFOLD(
        to_scaffold
            .map { h -> record(id: h.id, sample_id: h.meta.id, meta: h.meta, fasta: h.fasta) }
            .join(index.map { i -> record(id: i.id, index: i.index) }, by: 'id')
            .join(hic.map { r -> record(sample_id: r.id, reads: r.reads) }, by: 'sample_id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    reads: x.reads,
                    index: x.index,
                    fasta: x.fasta,
                    sort_bam: false,
                    bwamem2_args: '-SP -T0',
                    bwamem2_args2: '',
                    output_prefix: "${scaffoldPrefix(x.meta)}_${x.reads[0].name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')}"
                )
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
                    pairtools_args6: hic_bam_sorting == "by_coordinate" ? "--write-index" : "-n",
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
                    yahs_args: hic_type == "arima_v2" ? "-e GATC,GANTC,CTNAG,TTAA" : hic_type == "arima_v1" ? "-e GATC,GANTC" : "",
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
            def sorted = fastas.toSorted { a, b -> use_phased ? a.name <=> b.name : b.name <=> a.name }
            record(assembly_key: key, fastas: sorted)
        }
        .join(assemblies.map { a -> record(assembly_key: "${a.meta.id}:${a.meta.assembly.build}", id: a.id, meta: a.meta) }, by: 'assembly_key')

    assemblies_out = scaffolded.map { x ->
        record(id: x.id, meta: x.meta, pri_fasta: x.fastas[0], alt_fasta: x.fastas.size() > 1 ? x.fastas[1] : null, pri_gfa: null, alt_gfa: null)
    }

    MERGE_HAPLOTYPES(
        scaffolded
            .filter { x -> x.fastas.size() > 1 }
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    fastx: x.fastas.toSorted { a, b -> a.name <=> b.name },
                    out_ext: "fasta.gz",
                    seqkit_args: "-w 100",
                    output_prefix: "${x.meta.id}_${x.meta.assembly.build}"
                )
            }
    )

    emit:
    assemblies: Channel<Assembly> = assemblies_out
    logs      : Channel<Path>     = pairs.flatMap { p -> p.stat.toList() }
}
