nextflow.enable.types = true

include { BWAMEM2_INDEX                         } from "../../../modules/nf-core/bwamem2/index/main"
include { BWAMEM2_MEM as BWAMEM2_MEM_CURATION   } from "../../../modules/nf-core/bwamem2/mem/main"
include { FILTER_FIVE_END                       } from "../../../modules/local/hic_curation/filter_five_end"
include { TWOREADCOMBINER_FIXMATE_SORT          } from "../../../modules/local/hic_curation/tworeadcombiner_fixmate_sort"
include { BIOBAMBAM_BAMMARKDUPLICATES2          } from "../../../modules/nf-core/biobambam/bammarkduplicates2/main"
include { SAMTOOLS_FAIDX                        } from "../../../modules/nf-core/samtools/faidx/main"
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_HIC  } from "../../../modules/nf-core/samtools/merge/main"
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_HIFI } from "../../../modules/nf-core/samtools/merge/main"
include { BAM2BED_SORT                          } from "../../../modules/local/hic_curation/bam2bed_sort"
include { COOLER_CLOAD                          } from "../../../modules/nf-core/cooler/cload/main"
include { COOLER_ZOOMIFY                        } from "../../../modules/nf-core/cooler/zoomify/main"
include { CREATE_CHROMOSOME_SIZES_FILE          } from "../../../modules/local/hic_curation/create_chromosome_sizes"
include { PRETEXTMAP                            } from "../../../modules/nf-core/pretextmap/main"
include { MINIMAP2_INDEX                        } from "../../../modules/nf-core/minimap2/index/main"
include { MINIMAP2_ALIGN                        } from "../../../modules/nf-core/minimap2/align/main"
include { BAM2COVERAGE_TRACKS                   } from "../../../modules/local/hic_curation/bam2coverageTracks"
include { SEQTK_CUTN                            } from "../../../modules/nf-core/seqtk/cutn/main"
include { CREATE_GAP_TRACKS                     } from "../../../modules/local/hic_curation/create_gap_tracks"
include { TIDK_SEARCH as TIDK_SEARCH_BEDGRAPH   } from "../../../modules/nf-core/tidk/search/main"
include { TIDK_SEARCH as TIDK_SEARCH_TSV        } from "../../../modules/nf-core/tidk/search/main"
include { TIDK_PLOT                             } from "../../../modules/nf-core/tidk/plot/main"
include { CREATE_TELOMER_BIGWIG_TRACK           } from "../../../modules/local/hic_curation/create_telomer_bigwig_track"
include { PRETEXT_TRACKS_INGESTION              } from "../../../modules/local/hic_curation/pretext_tracks_ingestion"

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

// Identifies the sample and the assembly build
def assemblyKey(meta: Record) -> String {
    return "${meta.id}:${meta.assembly.build}"
}

// File name without the fastq or fasta extensions
def readsName(reads: Path) -> String {
    return reads.name.replaceAll(/\.(fastq|fq|fasta|fa)(\.gz)?$/, '')
}

workflow SCAFFOLD_CURATION {
    take:
    assemblies      : Channel<Assembly>
    hic             : Channel<Reads>  // Hi-C read pairs
    hifi            : Channel<Reads>  // HiFi reads
    hic_map_sort_by : String
    hic_map_qv      : Integer
    cooler_bin_size : Integer
    hifi_coverage_cap: Integer
    telomer_motif   : String

    main:
    primary = assemblies.map { a -> record(id: assemblyKey(a.meta), sample_id: a.meta.id, meta: a.meta, fasta: a.pri_fasta) }

    index = BWAMEM2_INDEX(
        primary.map { p -> record(id: p.id, meta: p.meta, fasta: p.fasta) }
    )

    faidx = SAMTOOLS_FAIDX(
        primary.map { p -> record(id: p.id, meta: p.meta, fasta: p.fasta, fai: null, get_sizes: false) }
    )
    fai = faidx.map { f -> record(id: f.id, fai: f.fai) }

    // Align each Hi-C read separately
    hic_reads = hic.flatMap { r ->
        r.reads.collect { f -> record(sample_id: r.id, pair_id: r.meta.pair_id, read: f) }.toList()
    }
    aligned = BWAMEM2_MEM_CURATION(
        primary
            .map { p -> record(id: p.id, sample_id: p.sample_id, meta: p.meta, fasta: p.fasta) }
            .join(index.map { i -> record(id: i.id, index: i.index) }, by: 'id')
            .join(hic_reads, by: 'sample_id')
            .map { x ->
                record(
                    id: "${x.id}:${x.pair_id}:${x.read.name}",
                    meta: x.meta + record(pair_id: x.pair_id),
                    reads: [ x.read ],
                    index: x.index,
                    fasta: x.fasta,
                    sort_bam: false,
                    bwamem2_args: '-SP -T0',
                    bwamem2_args2: '',
                    output_prefix: "${x.meta.id}_${x.meta.assembly.build}_${readsName(x.read)}"
                )
            }
    )

    // filter alignments
    filtered = FILTER_FIVE_END(
        aligned.map { a -> record(id: a.id, meta: a.meta, bam: a.bam, output_prefix: "${a.bam.baseName}_filt5prime") }
    )

    // combine the alignments of the two reads of each pair
    combined = TWOREADCOMBINER_FIXMATE_SORT(
        filtered
            .map { f -> tuple("${assemblyKey(f.meta)}:${f.meta.pair_id}", f) }
            .groupBy()
            .map { key, results ->
                def sorted = results.toSorted { a, b -> a.bam.name <=> b.bam.name }
                record(
                    id: key,
                    meta: sorted[0].meta,
                    bam: sorted.collect { s -> s.bam }.toList(),
                    combiner_args: "${hic_map_qv}",
                    combiner_args2: '--write-index -l1',
                    output_prefix: sorted[0].bam.baseName
                )
            }
    )

    // merge bam files in case multiple HIC paired-end libraries are present
    hic_groups = combined
        .map { c -> tuple(assemblyKey(c.meta), c) }
        .groupBy()
        .map { key, results -> record(id: key, meta: results.toList()[0].meta, bams: results.collect { r -> r.bam }.toList().toSorted()) }

    hic_merged = SAMTOOLS_MERGE_HIC(
        hic_groups
            .filter { g -> g.bams.size() > 1 }
            .map { g ->
                record(
                    id: g.id,
                    meta: g.meta,
                    input_files: g.bams,
                    index_files: [],
                    fasta: null,
                    fai: null,
                    gzi: null,
                    samtools_args: '--write-index',
                    output_prefix: "${g.meta.assembly.build}_merged"
                )
            }
    )

    dedup_bam = hic_groups
        .filter { g -> g.bams.size() == 1 }
        .map { g -> record(id: g.id, meta: g.meta, bam: g.bams[0]) }
        .mix(hic_merged.map { m -> record(id: m.id, meta: m.meta, bam: m.bam) })

    // dedupliucate bam file
    marked = BIOBAMBAM_BAMMARKDUPLICATES2(
        dedup_bam.map { d -> record(id: d.id, meta: d.meta, bam: d.bam, output_prefix: "${d.bam.baseName}_dupMarked") }
    )

    // convert bam to sorted bed file
    pairs = BAM2BED_SORT(
        marked.map { m ->
            record(
                id: m.id,
                meta: m.meta,
                bam: m.bam,
                bam2bed_args: '-u -F0x400 -h',          // arguments for samtools view
                bam2bed_args2: '-T \$PWD',              // arguments for sort
                bam2bed_args3: "${hic_map_qv}",         // min hic map quality
                output_prefix: m.bam.baseName
            )
        }
    )

    // create chrome sizes file from fai file
    chrom_sizes = CREATE_CHROMOSOME_SIZES_FILE(
        faidx.map { f -> record(id: f.id, meta: f.meta, fai: f.fai, hic_map_sort_by: hic_map_sort_by, output_prefix: f.meta.id) }
    )
    sizes = chrom_sizes.map { c -> record(id: c.id, chrom_sizes: c.sizes) }

    // create cooler files
    cool = COOLER_CLOAD(
        pairs
            .join(sizes, by: 'id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    contacts: x.pairs,
                    index: null,
                    chromsizes: x.chrom_sizes,
                    mode: 'pairs',
                    cool_bin: cooler_bin_size,
                    cooler_args: '-0 -c1 2 -p1 3 -c2 4 -p2 5',
                    output_prefix: x.pairs.baseName
                )
            }
    )

    zoomify = COOLER_ZOOMIFY(
        cool.map { c ->
            record(
                id: c.id,
                meta: c.meta,
                cool: c.cool,
                cooler_args: '--resolutions 1000,5000,10000,20000,40000,60000,80000,100000,120000,150000,200000,300000,400000,500000',
                output_prefix: c.cool.baseName
            )
        }
    )

    // create pretext maps
    pretext = PRETEXTMAP(
        marked
            .map { m -> record(id: m.id, meta: m.meta, input: m.bam) }
            .join(primary.map { p -> record(id: p.id, fasta: p.fasta) }, by: 'id')
            .join(fai, by: 'id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    input: x.input,
                    fasta: x.fasta,
                    fai: x.fai,
                    pretextmap_args: "--sortby ${hic_map_sort_by} --mapq ${hic_map_qv}",
                    pretextmap_args2: null,
                    output_prefix: x.input.baseName
                )
            }
    )

    // create tracks for PretextMap:
    // coverage, gap, telomer

    // Coverage-track: create minimap2 index
    MINIMAP2_INDEX(
        primary.map { p -> record(id: p.id, meta: p.meta, fasta: p.fasta, minimap2_args: '-x map-hifi') }
    )

    // Coverage-track: align HiFi reads
    hifi_aligned = MINIMAP2_ALIGN(
        primary
            .map { p -> record(id: p.id, sample_id: p.sample_id, meta: p.meta, fasta: p.fasta) }
            .join(hifi.map { r -> record(sample_id: r.id, reads: r.reads) }, by: 'sample_id')
            .map { x ->
                record(
                    id: "${x.id}:${x.reads[0].name}",
                    meta: x.meta,
                    reads: x.reads[0],
                    reference: x.fasta,
                    bam_format: true,
                    bam_index_extension: null,
                    cigar_paf_format: false,
                    cigar_bam: false,
                    minimap2_args: "--MD -x map-hifi",
                    minimap2_args2: null,
                    minimap2_args3: null,
                    minimap2_args4: null,
                    output_prefix: "${x.meta.assembly.build}_${readsName(x.reads[0])}"
                )
            }
    )

    // Coverage-track: merge hifi-bam files (if multiple hifi reads are present)
    hifi_groups = hifi_aligned
        .map { h -> tuple(assemblyKey(h.meta), h) }
        .groupBy()
        .map { key, results -> record(id: key, meta: results.toList()[0].meta, bams: results.collect { r -> r.bam }.toList().toSorted()) }

    hifi_merged = SAMTOOLS_MERGE_HIFI(
        hifi_groups
            .filter { g -> g.bams.size() > 1 }
            .map { g ->
                record(
                    id: g.id,
                    meta: g.meta,
                    input_files: g.bams,
                    index_files: [],
                    fasta: null,
                    fai: null,
                    gzi: null,
                    samtools_args: '--write-index',
                    output_prefix: "${g.meta.assembly.build}_merged"
                )
            }
    )

    hifi_bam = hifi_groups
        .filter { g -> g.bams.size() == 1 }
        .map { g -> record(id: g.id, meta: g.meta, bam: g.bams[0]) }
        .mix(hifi_merged.map { m -> record(id: m.id, meta: m.meta, bam: m.bam) })

    // Coverage-track: create hitile- and bed-coverage tracks
    coverage = BAM2COVERAGE_TRACKS(
        hifi_bam
            .join(sizes, by: 'id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    bam: x.bam,
                    chrom_sizes: x.chrom_sizes,
                    hifi_coverage_cap: hifi_coverage_cap,
                    coverage_args: '',               // samtools sort
                    coverage_args2: '-b -h -F 256',  // samtools view
                    coverage_args3: '-bga -split',   // genomeCoverageBed
                    coverage_args4: '-T \$PWD',      // sort
                    output_prefix: x.bam.baseName
                )
            }
    )

    // Gap-track: create bedtrack
    cutn = SEQTK_CUTN(
        primary.map { p -> record(id: p.id, meta: p.meta, fasta: p.fasta, seqtk_args: '-n 1', output_prefix: p.meta.assembly.build) }
    )

    // Gap-track: create beddb- and bed-gap tracks
    gaps = CREATE_GAP_TRACKS(
        cutn
            .join(sizes, by: 'id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    bed: x.bed,
                    chrom_sizes: x.chrom_sizes,
                    gap_args: null,
                    gap_args2: '--no-header',        // clodius aggregate
                    output_prefix: x.meta.assembly.build
                )
            }
    )

    // Telomer-track: create track
    telomer_bedgraph = TIDK_SEARCH_BEDGRAPH(
        primary.map { p ->
            record(id: p.id, meta: p.meta, fasta: p.fasta, string: telomer_motif, tidk_args: "--window 10000 --extension bedgraph", output_prefix: "${p.meta.assembly.build}_telomer")
        }
    )

    telomer_tsv = TIDK_SEARCH_TSV(
        primary.map { p ->
            record(id: p.id, meta: p.meta, fasta: p.fasta, string: telomer_motif, tidk_args: "--window 10000 --extension tsv", output_prefix: "${p.meta.assembly.build}_telomer")
        }
    )

    // Telomer-track: create plot - thats not really necessary, but nice to have
    TIDK_PLOT(
        telomer_tsv
            .filter { t -> t.tsv != null }
            .map { t -> record(id: t.id, meta: t.meta, tsv: t.tsv, output_prefix: t.meta.id) }
    )

    // Telomer-track: convert telomer bedgraph into beddb file that can be ingested into HiGlass
    telomer_bigwig = CREATE_TELOMER_BIGWIG_TRACK(
        telomer_bedgraph
            .filter { t -> t.bedgraph != null }
            .join(sizes, by: 'id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    bedgraph: x.bedgraph,
                    chrom_sizes: x.chrom_sizes,
                    telomer_args: '',                // sort
                    output_prefix: x.meta.assembly.build
                )
            }
    )

    // ingest coverage, gap and telomer track into Pretext
    // for each assembly always join:
    // 1-PretextMap with 1-CoverageTrack with 1-TelomereTrack with 1-GapTrack
    pretext_tracks = PRETEXT_TRACKS_INGESTION(
        pretext
            .map { p -> record(id: p.id, meta: p.meta, pretext_in: p.pretext) }
            .join(coverage.map { c -> record(id: c.id, cov_bedgraph: c.capped_bed) }, by: 'id')
            .join(telomer_bedgraph.filter { t -> t.bedgraph != null }.map { t -> record(id: t.id, telomer_bedgraph: t.bedgraph) }, by: 'id')
            .join(gaps.map { g -> record(id: g.id, gap_bedgraph: g.bedgraph) }, by: 'id')
            .map { x ->
                record(
                    id: x.id,
                    meta: x.meta,
                    pretext_in: x.pretext_in,
                    cov_bedgraph: x.cov_bedgraph,
                    telomer_bedgraph: x.telomer_bedgraph,
                    gap_bedgraph: x.gap_bedgraph,
                    output_prefix: x.meta.assembly.build
                )
            }
    )

    emit:
    zoomify       : Channel<Record> = zoomify
    coverage      : Channel<Record> = coverage
    gaps          : Channel<Record> = gaps
    telomer_bigwig: Channel<Record> = telomer_bigwig
    pretext_tracks: Channel<Record> = pretext_tracks
}
