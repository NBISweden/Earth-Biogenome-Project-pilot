include { BLOBTOOLKIT_ADD    } from "$projectDir/modules/local/blobtoolkit/add"
include { BLOBTOOLKIT_CREATE } from "$projectDir/modules/local/blobtoolkit/create"
include { BLOBTOOLKIT_VIEW   } from "$projectDir/modules/local/blobtoolkit/view"
include { BLAST_BLASTN       } from "$projectDir/modules/nf-core/modules/blast/blastn/main"
include { DIAMOND_BLASTX     } from "$projectDir/modules/nf-core/modules/diamond/blastx/main"
include { MINIMAP2_ALIGN     } from "$projectDir/modules/nf-core/modules/minimap2/align/main"
include { BUSCO              } from "$projectDir/modules/nf-core/modules/busco/main"
include { SEQKIT_SPLIT2      } from "$projectDir/modules/nf-core/modules/seqkit/split2/main"

workflow BLOBTOOLKIT {
    take:
    read_assembly_ch   // input type: [ [id: 'name'], [ file(read1), file(read2) ], file(assembly) ]
    busco_lineages     // Busco lineages to check against
    busco_lineage_path // Path to Busco lineage files
    uniprot_db         // Path to Uniprot database
    ncbi_nt_db         // Path the ncbi nt database
    ncbi_taxonomy      // Path to ncbi taxonomy database

    main:
    // Blobtoolkit workflow : See README.md
    input = read_assembly_ch.multiMap { sample, reads, assembly -> 
        read_ch: [ sample, reads ] 
        asm_ch:  [ sample, assembly ] 
    }

    // Split assembly for faster taxonomic assignment
    SEQKIT_SPLIT2 ( input.asm_ch )
    versions_ch = SEQKIT_SPLIT2.out.versions
    
    // Generate blob DB
    BLOBTOOLKIT_CREATE ( 
        input.asm_ch, 
        [] // ignore blobtools meta file for now 
    )
    versions_ch = versions_ch.mix(BLOBTOOLKIT_CREATE.out.versions)

    // Taxonomic hits
    BLAST_BLASTN ( 
        SEQKIT_SPLIT2.out.reads.transpose(), 
        ncbi_nt_db 
    )
    DIAMOND_BLASTX ( 
        SEQKIT_SPLIT2.out.reads.transpose(),
        uniprot_db,
        "txt",
        'qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    )
    versions_ch = versions_ch.mix( BLAST_BLASTN.out.versions.first(), DIAMOND_BLASTX.out.versions.first() )

    // Coverage
    MINIMAP2_ALIGN ( 
        read_assembly_ch.multiMap { sample, reads, assembly -> 
                fastx_ch : [ sample, reads ]
                reference: assembly
            }, // supplies reads channel and reference channel
        true, // use bam format
        false, // use cigar_paf_format
        false // use cigar_bam
    )
    versions_ch = versions_ch.mix( MINIMAP2_ALIGN.out.versions )

    // Busco
    BUSCO (
        input.asm_ch,
        busco_lineages,
        busco_lineage_path,
        []
    )

    // Plot
    BLOBTOOLKIT_ADD (
        BLOBTOOLKIT_CREATE.out.blobdb,
        BLAST_BLASTN.out.txt.mix(DIAMOND_BLASTX.out.txt).collect(),
        MINIMAP2_ALIGN.out.bam,
        BUSCO.out.busco_dir,
        [], // bed - ignore
        [], // beddir - ignore
        [], // bedtsv - ignore
        ncbi_taxonomy
    )
    BLOBTOOLKIT_VIEW (
        BLOBTOOLKIT_ADD.out.blobdb,
        [ 'blob', 'cumulative', 'snail' ]
    )
    versions_ch = versions_ch.mix( BLOBTOOLKIT_ADD.out.versions, BLOBTOOLKIT_VIEW.out.versions )

    emit:
    blobplots = BLOBTOOLKIT_VIEW.out.png.mix(BLOBTOOLKIT_VIEW.out.svg)
    versions = versions_ch
}
