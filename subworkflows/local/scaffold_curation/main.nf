include { constructAssemblyRecord               } from "$projectDir/modules/local/functions"
include { getPrimaryAssembly                    } from "$projectDir/modules/local/functions"
include { combineByMetaKeys                     } from "$projectDir/modules/local/functions"
include { joinByMetaKeys                        } from "$projectDir/modules/local/functions"
include { BWAMEM2_INDEX                         } from "$projectDir/modules/nf-core/bwamem2/index/main"
include { BWAMEM2_MEM as BWAMEM2_MEM_CURATION   } from "$projectDir/modules/nf-core/bwamem2/mem/main"
include { FILTER_FIVE_END                       } from "$projectDir/modules/local/hic_curation/filter_five_end"
include { TWOREADCOMBINER_FIXMATE_SORT          } from "$projectDir/modules/local/hic_curation/tworeadcombiner_fixmate_sort"
include { BIOBAMBAM_BAMMARKDUPLICATES2          } from "$projectDir/modules/nf-core/biobambam/bammarkduplicates2/main"
include { SAMTOOLS_FAIDX                        } from "$projectDir/modules/nf-core/samtools/faidx/main"
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_HIC  } from "$projectDir/modules/nf-core/samtools/merge/main"
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_HIFI } from "$projectDir/modules/nf-core/samtools/merge/main"
include { SAMTOOLS_SORT                         } from "$projectDir/modules/nf-core/samtools/sort/main"
include { BAM2BED_SORT                          } from "$projectDir/modules/local/hic_curation/bam2bed_sort"
include { COOLER_CLOAD                          } from "$projectDir/modules/nf-core/cooler/cload/main"
include { COOLER_ZOOMIFY                        } from "$projectDir/modules/nf-core/cooler/zoomify/main"
include { CREATE_CHROMOSOME_SIZES_FILE          } from "$projectDir/modules/local/hic_curation/create_chromosome_sizes"
include { PRETEXTMAP                            } from "$projectDir/modules/nf-core/pretextmap/main"   
include { MINIMAP2_INDEX                        } from "$projectDir/modules/nf-core/minimap2/index/main"
include { MINIMAP2_ALIGN                        } from "$projectDir/modules/nf-core/minimap2/align/main"
include { BAM2COVERAGE_TRACKS                   } from "$projectDir/modules/local/hic_curation/bam2coverageTracks"
include { SEQTK_CUTN                            } from "$projectDir/modules/nf-core/seqtk/cutn/main"
include { CREATE_GAP_TRACKS                     } from "$projectDir/modules/local/hic_curation/create_gap_tracks"
include { TIDK_SEARCH as TIDK_SEARCH_BEDGRAPH   } from "$projectDir/modules/nf-core/tidk/search/main"
include { TIDK_SEARCH as TIDK_SEARCH_TSV        } from "$projectDir/modules/nf-core/tidk/search/main"
include { TIDK_PLOT                             } from "$projectDir/modules/nf-core/tidk/plot/main"
include { CREATE_TELOMER_HITILE_TRACK           } from "$projectDir/modules/local/hic_curation/create_telomer_hitile_track"
include { PRETEXT_TRACKS_INGESTION              } from "$projectDir/modules/local/hic_curation/pretext_tracks_ingestion"

workflow SCAFFOLD_CURATION {
    take:
    ch_assemblies // [ meta, assembly ]
    ch_hic        // [ meta, hic-pairs ]
    ch_hifi       // [ meta, hifi-reads ]

    main:

    ch_versions  = Channel.empty()

    BWAMEM2_INDEX ( getPrimaryAssembly(ch_assemblies) )
    ch_versions  = ch_versions.mix( BWAMEM2_INDEX.out.versions )

    SAMTOOLS_FAIDX (
        getPrimaryAssembly(ch_assemblies), // [meta, fasta]
        [ [ ], [ ] ]                       // [meta2, fai] 
    )
    ch_versions  = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )

    combineByMetaKeys( // Combine (Hi-C + index) with Assembly
        combineByMetaKeys( // Combine Hi-C reads with BWA index
            ch_hic, BWAMEM2_INDEX.out.index,
            keySet: ['id','sample'],
            meta: 'merge'
        ),
        getPrimaryAssembly( ch_assemblies ),
        keySet: ['id','sample'],
        meta: 'merge'
    ).transpose(by:1)       // by meta info: [id, sample, settings, single_end, pair_id, assembly]
    .multiMap{ meta, hic_reads, index, fasta ->
        reads: [ meta, hic_reads ]
        index: [ meta, index ]
        fasta: [ meta, fasta ]
    }.set{ bwamem2_input }

    BWAMEM2_MEM_CURATION (bwamem2_input.reads, bwamem2_input.index, bwamem2_input.fasta, false)
    ch_versions  = ch_versions.mix( BWAMEM2_MEM_CURATION.out.versions )

    // filter alignments 
    FILTER_FIVE_END(BWAMEM2_MEM_CURATION.out.bam)
    ch_versions  = ch_versions.mix( FILTER_FIVE_END.out.versions )

    // combine reads 
    FILTER_FIVE_END.out.bam
    .map { meta, bam -> [ meta.subMap('id', 'sample', 'assembly', 'pair_id'), bam ] }
    .groupTuple( sort: {a, b -> a.name <=> b.name } )
    .set { combine_input }

    TWOREADCOMBINER_FIXMATE_SORT(combine_input)
    ch_versions  = ch_versions.mix( TWOREADCOMBINER_FIXMATE_SORT.out.versions )

    // merge bam files in case multiple HIC paired-end libraries are present 
    TWOREADCOMBINER_FIXMATE_SORT.out.bam
        .map { meta, bam_list -> [ meta - meta.subMap( 'pair_id' ), bam_list ] }
        .groupTuple()
        .branch { meta, bam_list ->
            multiples: bam_list.size() > 1
            singleton: true
        }
    .set { merge_bam }

    SAMTOOLS_MERGE_HIC(
        merge_bam.multiples,
        [ [], [] ],             // [meta2, fasta]
        [ [], [] ]              // [meta2, fai]
    )
    ch_versions  = ch_versions.mix( SAMTOOLS_MERGE_HIC.out.versions )

    merge_bam.singleton
        .map { meta, bam -> [ meta, *bam ] } // the spread operator (*) flattens the bam list
        .mix( SAMTOOLS_MERGE_HIC.out.bam)
        .set { dedup_bam }

    // dedupliucate bam file 
    BIOBAMBAM_BAMMARKDUPLICATES2(dedup_bam)
    ch_versions  = ch_versions.mix( BIOBAMBAM_BAMMARKDUPLICATES2.out.versions )

    // convert bam to sorted bed file
    BAM2BED_SORT(BIOBAMBAM_BAMMARKDUPLICATES2.out.bam)
    ch_versions  = ch_versions.mix( BAM2BED_SORT.out.versions )

    // create chrome sizes file from fai file
    CREATE_CHROMOSOME_SIZES_FILE(SAMTOOLS_FAIDX.out.fai)
    ch_versions  = ch_versions.mix( CREATE_CHROMOSOME_SIZES_FILE.out.versions )

    // create cooler files 
    // tuple val(meta), path(pairs), path(index), val(cool_bin)
    // path chromsizes
    BAM2BED_SORT.out.pairs
        .map { meta, pairs  -> [ meta, pairs, [ ] ] }
        .combine(
            Channel.value(params.cooler_bin_size)
        )
        .set { pairs_idx_binsize_ch }
    
    combineByMetaKeys( // Combine Hi-C reads with BWA index
        pairs_idx_binsize_ch, CREATE_CHROMOSOME_SIZES_FILE.out.sizes,
        keySet: ['id','sample'],
        meta: 'rhs'
    )
    .multiMap { meta, pairs, fake_index, cool_bin, chrom_sizes -> 
        cload_in     : [ meta, pairs, fake_index, cool_bin ] 
        chrom_sizes  : chrom_sizes
    }
    .set { cooler_cload_ch }

    COOLER_CLOAD(
        cooler_cload_ch.cload_in,
        cooler_cload_ch.chrom_sizes
    )
    ch_versions  = ch_versions.mix( COOLER_CLOAD.out.versions )

    COOLER_CLOAD.out.cool.map { meta, cool, cool_bin ->  [ meta, cool ] }
    .set { cooler_zoomify_ch }
    
    COOLER_ZOOMIFY(cooler_zoomify_ch)
    ch_versions  = ch_versions.mix( COOLER_ZOOMIFY.out.versions )

    // create pretext maps 
    joinByMetaKeys(joinByMetaKeys(
            getPrimaryAssembly(ch_assemblies),
            SAMTOOLS_FAIDX.out.fai,
            keySet: ['id','sample'],
            meta: 'rhs'
        ), BIOBAMBAM_BAMMARKDUPLICATES2.out.bam,
        keySet: ['id','sample'],
        meta: 'rhs'
    )
    .multiMap { meta, fasta, fai, bam ->
        bam       : [ meta, bam ]
        fasta_fai : [ meta, fasta, fai ]
    }
    .set { pretext_ch }
    
    PRETEXTMAP(
        pretext_ch.bam,
        pretext_ch.fasta_fai
    )
    ch_versions  = ch_versions.mix( PRETEXTMAP.out.versions )

    // create tracks for PretextMap:
    // coverage, gap, telomer
    
    joinByMetaKeys( 
        ch_hifi,
        getPrimaryAssembly( ch_assemblies ),
        keySet: ['id','sample'],
        meta: 'rhs'
    )
    .multiMap { meta, hifi_reads, fasta -> 
        hifi_reads: [ meta, hifi_reads ]
        reference : fasta 
    }
    .set{ asm_hifi_ch }

    // Coverage-track: create minimap2 index 
    MINIMAP2_INDEX(getPrimaryAssembly( ch_assemblies ))
    ch_versions  = ch_versions.mix( MINIMAP2_INDEX.out.versions )

    // Coverage-track: align HiFi reads 
    MINIMAP2_ALIGN(
        asm_hifi_ch.hifi_reads, // [ meta, reads ]
        asm_hifi_ch.reference,  // reference 
        1,                      // bam_format 
        0,                      // cigar_paf_format
        0                       // cigar_ba, 
    )
    ch_versions  = ch_versions.mix( MINIMAP2_ALIGN.out.versions )

    // Coverage-track: merge hifi-bam files (if multiple hifi reads are present)
    MINIMAP2_ALIGN.out.bam
        .branch { meta, bam_list ->
            multiples: bam_list.size() > 1
            singleton: true
        }
    .set { merge_hifi_bam }

    SAMTOOLS_MERGE_HIFI(
        merge_hifi_bam.multiples,
        [ [], [] ],
        [ [], [] ]
    )
    ch_versions  = ch_versions.mix( SAMTOOLS_MERGE_HIFI.out.versions )

    joinByMetaKeys( 
        ch_hifi,
        getPrimaryAssembly( ch_assemblies ),
        keySet: ['id','sample'],
        meta: 'rhs'
    )
    .multiMap { meta, hifi_reads, fasta -> 
        hifi_reads: [ meta, hifi_reads ]
        reference : fasta 
    }
    .set{ asm_hifi_ch }

    joinByMetaKeys( 
        merge_hifi_bam.singleton
            .map { meta, bam -> [ meta, *bam ] } // the spread operator (*) flattens the bam list
            .mix( SAMTOOLS_MERGE_HIFI.out.bam),
        CREATE_CHROMOSOME_SIZES_FILE.out.sizes,
        keySet: ['id','sample'],
        meta: 'rhs'
    )
    .multiMap { meta, bam, chrom_sizes -> 
        merged_bam   : [ meta, bam ] 
        chrom_sizes  : chrom_sizes
    }
    .set { bam2coverage_ch }

    // Coverage-track: create hitile- and bed-coverage tracks 
    BAM2COVERAGE_TRACKS(
        bam2coverage_ch.merged_bam,
        bam2coverage_ch.chrom_sizes
    )
    ch_versions  = ch_versions.mix( BAM2COVERAGE_TRACKS.out.versions )

    // Gap-track: create bedtrack 
    SEQTK_CUTN(getPrimaryAssembly( ch_assemblies ))
    ch_versions  = ch_versions.mix( SEQTK_CUTN.out.versions )

    // Gap-track: create beddb- and bed-gap tracks 
    joinByMetaKeys( 
        SEQTK_CUTN.out.bed,
        CREATE_CHROMOSOME_SIZES_FILE.out.sizes,
        keySet: ['id','sample'],
        meta: 'rhs'
    )
    .multiMap { meta, bed, chrom_sizes -> 
        bed          : [ meta, bed ] 
        chrom_sizes  : chrom_sizes
    }
    .set { bed2gap_ch }

    CREATE_GAP_TRACKS(
        bed2gap_ch.bed,
        bed2gap_ch.chrom_sizes
    ) 
    ch_versions  = ch_versions.mix( CREATE_GAP_TRACKS.out.versions )

    // Telomer-track: create track
    TIDK_SEARCH_BEDGRAPH(
        getPrimaryAssembly( ch_assemblies ),
        params.telomer_motif
    )
    ch_versions  = ch_versions.mix( TIDK_SEARCH_BEDGRAPH.out.versions )

    TIDK_SEARCH_TSV(
        getPrimaryAssembly( ch_assemblies ),
        params.telomer_motif
    )
    ch_versions  = ch_versions.mix( TIDK_SEARCH_TSV.out.versions )

    // Telomer-track: create plot - thats not really necessary, but nice to have 
    TIDK_PLOT(TIDK_SEARCH_TSV.out.tsv)
    ch_versions  = ch_versions.mix( TIDK_PLOT.out.versions )

    // Telomer-track: convert telomer bedgraph into beddb file that can be ingested into HiGlass 
    joinByMetaKeys( 
        TIDK_SEARCH_BEDGRAPH.out.bedgraph,
        CREATE_CHROMOSOME_SIZES_FILE.out.sizes,
        keySet: ['id','sample'],
        meta: 'rhs'
    )
    .multiMap { meta, bedgraph, chrom_sizes -> 
        bedgraph          : [ meta, bedgraph ] 
        chrom_sizes  : chrom_sizes
    }
    .set { bedgraph_telomer_ch }

    CREATE_TELOMER_HITILE_TRACK(
        bedgraph_telomer_ch.bedgraph,
        bedgraph_telomer_ch.chrom_sizes
    )
    ch_versions  = ch_versions.mix( CREATE_TELOMER_HITILE_TRACK.out.versions )

    // ingest coverage, gap and telomer track into Pretext
    joinByMetaKeys(
        joinByMetaKeys(
            joinByMetaKeys(
                PRETEXTMAP.out.pretext,
                BAM2COVERAGE_TRACKS.out.capped_bed,
                keySet: ['id','sample', 'assembly'],
                meta: 'rhs'
            ), 
            TIDK_SEARCH_BEDGRAPH.out.bedgraph,
            keySet: ['id','sample','assembly'],
            meta: 'rhs'
        ),
        CREATE_GAP_TRACKS.out.bedgraph,
        keySet: ['id','sample','assembly'],
        meta: 'rhs'
    )    
    .set { pretext_tracks_ch }

    PRETEXT_TRACKS_INGESTION(pretext_tracks_ch)    
    ch_versions  = ch_versions.mix( PRETEXT_TRACKS_INGESTION.out.versions )

    emit:
    assemblies = Channel.empty()    
    versions   = ch_versions
}