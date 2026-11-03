#!/usr/bin/env nextflow

/*
 * Polishing of homozygous variants with DeepVariant, followed by assembly QC.
 * https://git.mpi-cbg.de/assembly/programs/polishing
 *
 * Takes one wide assembly record per assembly build (see types/records.nf). The HiFi reads are
 * read from `hifi_merged` and the HiFi Meryl database from `meryl_hifi_db`. Returns the polished
 * assemblies with the same fields.
 *
 * outline:
 *
 * |- create pbmm2 index for assembly (1)
 * |- create bed chunks for given assembly (1..n)
 * |- align reads to full assembly (1)
 *     \- split alignment by bed chunk regions (1..n)
 *     |- index chunk alignment files (1..n)
 *     |- call variants with DeepVariant (1..n)
 *     |- filter variants (PASS + homozygous) (1..n)
 *     |- create tabix index files  (1..n)
 *     |- merge all variants (1)
 *     |- create consensus sequence (1)
 */

nextflow.enable.types = true

include { DVPOLISH_CHUNKFA                        } from "../../modules/local/dvpolish/chunkfa"
include { DVPOLISH_PBMM2_INDEX                    } from "../../modules/local/dvpolish/pbmm2_index"
include { DVPOLISH_PBMM2_ALIGN                    } from "../../modules/local/dvpolish/pbmm2_align"
include { SAMTOOLS_FAIDX                          } from "../../modules/nf-core/samtools/faidx/main"
include { SAMTOOLS_VIEW                           } from "../../modules/nf-core/samtools/view/main"
include { DEEPVARIANT_RUNDEEPVARIANT              } from "../../modules/nf-core/deepvariant/rundeepvariant/main"
include { BCFTOOLS_VIEW                           } from "../../modules/nf-core/bcftools/view/main"
include { TABIX_TABIX as TABIX_TABIX              } from "../../modules/nf-core/tabix/tabix/main"
include { TABIX_TABIX as TABIX_TABIX_MERGED       } from "../../modules/nf-core/tabix/tabix/main"
include { BCFTOOLS_MERGE                          } from "../../modules/nf-core/bcftools/merge/main"
include { BCFTOOLS_CONSENSUS                      } from "../../modules/nf-core/bcftools/consensus/main"
include { MERQURY_MERQURY as MERQURY_INPUT_ASM    } from "../../modules/nf-core/merqury/merqury/main"
include { MERQURY_MERQURY as MERQURY_POLISHED_ASM } from "../../modules/nf-core/merqury/merqury/main"
include { DVPOLISH_CREATE_FINALASM                } from "../../modules/local/dvpolish/createFinalAsm"
include { workflow as ASSEMBLY_QC_POLISHED        } from "../assembly_qc"
include { Assembly; QcParams; defaultQcParams; withFastas } from "../../types/records"

// One haplotype of an assembly. The id identifies the sample, the assembly build and the haplotype.
record Haplotype {
    id      : String
    meta    : Record  // includes meta.assembly and meta.haplotype
    fasta   : Path
    reads   : Path    // All HiFi reads of the sample in one file
    meryl_db: Path    // HiFi Meryl database of the sample
}

params {
    assemblies: Channel<Assembly>  // Assemblies to polish
    provided  : Channel<Assembly>? // Assemblies of this stage that were supplied. Evaluated, not polished.

    chunk_size: MemoryUnit = 100.MB // Size of the regions that are called separately, e.g. 1.5 GB

    // Assembly QC
    qc        : QcParams    = defaultQcParams() // Options of the assembly QC. use_phased: also process and evaluate the alternate assembly
}

def haplotypeKey(meta: Record) -> String {
    return "${meta.id}:${meta.assembly.build}:${meta.haplotype}"
}

def haplotype(a: Record, hap: String, fasta: Path) -> Haplotype {
    def meta = a.meta + record(haplotype: hap)
    return record(id: haplotypeKey(meta), meta: meta, fasta: fasta, reads: a.hifi_merged, meryl_db: a.meryl_hifi_db)
}

workflow {

    main:
    def chunk_size_bytes = params.chunk_size.toBytes()

    // Create tagged channels for each available haplotype
    haplotypes = params.assemblies.flatMap { a ->
        // A primary-only assembly has no second haplotype, also when phased
        def haps = params.qc.use_phased && a.alt_fasta
            ? [ haplotype(a, 'hap1', a.pri_fasta), haplotype(a, 'hap2', a.alt_fasta) ]
            : [ haplotype(a, 'hap1', a.pri_fasta) ]
        haps.toList()
    }

    // Index assembly files
    faidx = SAMTOOLS_FAIDX(
        haplotypes.map { h -> record(id: h.id, meta: h.meta, fasta: h.fasta, fai: null, get_sizes: false) }
    )
    fai = faidx.map { f -> record(id: f.id, fai: f.fai) }

    // Generate BED intervals partitioning each assembly into fixed-size regions
    chunks = DVPOLISH_CHUNKFA(
        faidx.map { f -> record(id: f.id, meta: f.meta, fai: f.fai, chunk_size_bytes: chunk_size_bytes, output_prefix: f.meta.id) }
    )
    beds = chunks.flatMap { c -> c.bed.toList().collect { bed -> record(id: c.id, bed: bed) }.toList() }

    // Create minimap2 index for each assembly
    index = DVPOLISH_PBMM2_INDEX(
        haplotypes.map { h -> record(id: h.id, meta: h.meta, fasta: h.fasta) }
    )

    // Map reads with pbmm2
    aligned = DVPOLISH_PBMM2_ALIGN(
        haplotypes
            .map { h -> record(id: h.id, meta: h.meta, reads: h.reads) }
            .join(index.map { i -> record(id: i.id, index: i.index) }, by: 'id')
            .map { x -> record(id: x.id, meta: x.meta, reads: x.reads, index: x.index, pbmm2_args: '-N 1 -l 4000') }
    )

    // Combine indexed bam with bed chunk file, tag by region (mergeID) and emit each region separately
    chunk_inputs = aligned
        .map { a -> record(id: a.id, meta: a.meta, bam: a.bam, bai: a.bai) }
        .join(beds, by: 'id')
        .join(haplotypes.map { h -> record(id: h.id, fasta: h.fasta) }, by: 'id')
        .join(fai, by: 'id')
        .map { x ->
            record(
                id: "${x.id}:${x.bed.baseName}",
                meta: x.meta + record(mergeID: x.bed.baseName),
                bam: x.bam,
                bai: x.bai,
                bed: x.bed,
                fasta: x.fasta,
                fai: x.fai
            )
        }

    // Split alignment bam by BED coordinates
    chunk_bams = SAMTOOLS_VIEW(
        chunk_inputs.map { c ->
            record(
                id: c.id,
                meta: c.meta,
                input: c.bam,
                index: c.bai,
                fasta: null,
                fai: null,
                bed: c.bed,
                index_format: 'bai',
                samtools_args: "-bh -F 2308 -M",
                output_prefix: "${c.bam.baseName}_${c.meta.mergeID}"
            )
        }
    )

    // run deepvariant on the chunked bam files
    calls = DEEPVARIANT_RUNDEEPVARIANT(
        chunk_bams
            .map { v -> record(id: v.id, meta: v.meta, bam: v.bam, bai: v.bai) }
            .join(chunk_inputs.map { c -> record(id: c.id, bed: c.bed, fasta: c.fasta, fai: c.fai) }, by: 'id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    input: x.bam,
                    index: x.bai,
                    intervals: x.bed,
                    fasta: x.fasta,
                    fai: x.fai,
                    deepvariant_args: "--model_type=PACBIO",
                    output_prefix: x.bam.baseName
                )
            }
    )

    // filter vcf files for PASS and homozygous variants
    // TODO add a minimim and maximum coverage filter ??? Needs to be tested
    filtered = BCFTOOLS_VIEW(
        calls.map { c ->
            record(
                id: c.id,
                meta: c.meta,
                vcf: c.vcf,
                index: c.vcf_tbi,
                samples: null,
                bcftools_args: "-f 'PASS' -i 'GT=\"1/1\"' --no-version -Oz",
                output_prefix: "${c.vcf.name.replaceAll(/\.vcf\.gz$/, '')}_filt"
            )
        }
    )

    // Index vcf files
    indexed = TABIX_TABIX(
        filtered.map { f -> record(id: f.id, meta: f.meta, tab: f.vcf) }
    )

    // in case of multiple vcf files, merge them prior to the consensus step
    vcf_groups = filtered
        .map { f -> tuple(haplotypeKey(f.meta), f.vcf) }
        .groupBy()
        .map { key, vcfs -> record(id: key, vcfs: vcfs.toSorted()) }
        .join(
            indexed
                .map { t -> tuple(haplotypeKey(t.meta), t.index) }
                .groupBy()
                .map { key, tbis -> record(id: key, tbis: tbis.toSorted()) },
            by: 'id'
        )

    merged = BCFTOOLS_MERGE(
        vcf_groups
            .filter { g -> g.vcfs.size() > 1 }
            .join(haplotypes.map { h -> record(id: h.id, meta: h.meta, fasta: h.fasta) }, by: 'id')
            .join(fai, by: 'id')
            .map { g ->
                record(
                    id: g.id,
                    meta: g.meta,
                    vcfs: g.vcfs,
                    tbis: g.tbis,
                    bed: null,
                    fasta: g.fasta,
                    fai: g.fai,
                    bcftools_args: "--force-samples -Oz",
                    output_prefix: g.meta.id
                )
            }
    )

    // index merged vcf file
    merged_index = TABIX_TABIX_MERGED(
        merged.map { m -> record(id: m.id, meta: m.meta, tab: m.vcf) }
    )

    // Prepare input for consensus step
    vcf_plus_index = vcf_groups
        .filter { g -> g.vcfs.size() == 1 }
        .map { g -> record(id: g.id, vcf: g.vcfs[0], tbi: g.tbis[0]) }
        .mix(
            merged
                .map { m -> record(id: m.id, vcf: m.vcf) }
                .join(merged_index.map { t -> record(id: t.id, tbi: t.index) }, by: 'id')
        )

    // create consensus sequence
    consensus = BCFTOOLS_CONSENSUS(
        vcf_plus_index
            .join(haplotypes.map { h -> record(id: h.id, meta: h.meta, fasta: h.fasta) }, by: 'id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    vcf: x.vcf,
                    tbi: x.tbi,
                    fasta: x.fasta,
                    output_prefix: "${x.meta.id}_consensus"
                )
            }
    )

    // run merqury on input assembly
    input_qv = MERQURY_INPUT_ASM(
        haplotypes.map { h -> record(id: h.id, meta: h.meta, meryl_db: h.meryl_db, assembly: [ h.fasta ], output_prefix: "${h.meta.assembly.build}_merqury") }
    )

    // run merqury on polished assembly
    polished_qv = MERQURY_POLISHED_ASM(
        consensus
            .map { c -> record(id: c.id, meta: c.meta, fasta: c.fasta) }
            .join(haplotypes.map { h -> record(id: h.id, meryl_db: h.meryl_db) }, by: 'id')
            .map { x -> record(id: x.id, meta: x.meta, meryl_db: x.meryl_db, assembly: [ x.fasta ], output_prefix: "${x.meta.assembly.build}_merqury") }
    )

    // select the best contigs
    final_asm = DVPOLISH_CREATE_FINALASM(
        haplotypes
            .map { h -> record(id: h.id, meta: h.meta, unpol_fasta: h.fasta) }
            .join(input_qv.map { q -> record(id: q.id, unpol_merqury_csv: q.scaffold_qv.toList()[0]) }, by: 'id')
            .join(consensus.map { c -> record(id: c.id, pol_fasta: c.fasta) }, by: 'id')
            .join(polished_qv.map { q -> record(id: q.id, pol_merqury_csv: q.scaffold_qv.toList()[0]) }, by: 'id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    unpol_fasta: x.unpol_fasta,
                    unpol_merqury_csv: x.unpol_merqury_csv,
                    pol_fasta: x.pol_fasta,
                    pol_merqury_csv: x.pol_merqury_csv,
                    output_prefix: "${x.unpol_fasta.name.replaceAll(/\.[^.]+\.[^.]+$/, '')}_dvpolish"
                )
            }
    )

    // Recombine the haplotypes of each assembly into the wide record
    polished_assemblies = final_asm
        .map { f -> tuple("${f.meta.id}:${f.meta.assembly.build}", f) }
        .groupBy()
        .map { key, results ->
            def sorted = results.toSorted { a, b -> params.qc.use_phased ? a.fasta_gz.name <=> b.fasta_gz.name : b.fasta_gz.name <=> a.fasta_gz.name }
            record(assembly_key: key, fastas: sorted.collect { f -> f.fasta_gz }.toList())
        }
        .join(params.assemblies.map { a -> record(assembly_key: "${a.meta.id}:${a.meta.assembly.build}", assembly: a) }, by: 'assembly_key')
        .map { x -> withFastas(x.assembly, x.fastas) }
    assemblies = polished_assemblies.mix(params.provided ?: channel.empty())

    qc = ASSEMBLY_QC_POLISHED(
        record(assemblies: assemblies) + params.qc
    )

    publish:
    assemblies = assemblies
    dvpolish   = final_asm
    merqury_input    = input_qv
    merqury_polished = polished_qv
    merquryfk  = qc.merquryfk
    merqury    = qc.merqury
    busco      = qc.busco
    gfastats   = qc.gfastats
    logs       = final_asm.map { f -> f.dvpolish_selection_tsv }.mix(qc.logs)
}

output {
    assemblies: Channel<Assembly> {
        path { a ->
            a.pri_fasta >> ("assemblies/${a.id}/")
            a.alt_fasta >> ("assemblies/${a.id}/")
        }
    }
    dvpolish: Channel<Record> {
        path { r -> "dvpolish/${r.meta.assembly.build}" }
    }
    // Merqury of the haplotypes before polishing and of the consensus, which decide the contig selection
    merqury_input: Channel<Record> {
        path { r -> "dvpolish/${r.meta.assembly.build}/merqury_input" }
    }
    merqury_polished: Channel<Record> {
        path { r -> "dvpolish/${r.meta.assembly.build}/merqury_polished" }
    }
    merquryfk: Channel<Record> {
        path { r -> "merquryfk/${r.meta.assembly.build}" }
    }
    merqury: Channel<Record> {
        path { r -> "merqury/${r.meta.assembly.build}" }
    }
    busco: Channel<Record> {
        path { r -> "busco/${r.meta.assembly.build}" }
    }
    gfastats: Channel<Record> {
        path { r -> "gfastats/${r.meta.assembly.build}" }
    }
    logs: Channel<Path> {
        path 'logs'
    }
}
