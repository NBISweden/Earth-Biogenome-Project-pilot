// include { SEQKIT_SPLIT2      } from "$projectDir/modules/local/seqkit/seqkit_split2"
include { BLAST_BLASTN       } from "$projectDir/modules/nf-core/modules/blast/blastn/main"
include { DIAMOND_BLASTX     } from "$projectDir/modules/nf-core/modules/diamond/blastx/main"

workflow ASSEMBLY_CONTAMINATION_SCREEN {
    take:
    assembly_ch
    blast_db         // Paths to blast databases
    diamond_db       // Paths to diamond databases

    main:
    // Split assembly for faster taxonomic assignment
    // SEQKIT_SPLIT2 ( assembly_ch )
    // versions_ch = SEQKIT_SPLIT2.out.versions
    versions_ch = Channel.empty()

    // Taxonomic hits
    BLAST_BLASTN ( 
        // SEQKIT_SPLIT2.out.reads
            // .transpose()
        assembly_ch
            .combine( blast_db )
            .multiMap { meta, assembly, database ->
                queries: [ meta, assembly.primary_asm_path ]
                db: database
            }
    )
    DIAMOND_BLASTX ( 
        // SEQKIT_SPLIT2.out.reads
            // .transpose()
        assembly_ch
            .combine( diamond_db )
            .multiMap { meta, assembly, database ->
                queries: [ meta, assembly.primary_asm_path ]
                db: database
            },
        "txt",
        'qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    )
    versions_ch = versions_ch.mix( BLAST_BLASTN.out.versions.first(), DIAMOND_BLASTX.out.versions.first() )
    
    emit:
    blast_hits   = BLAST_BLASTN.out.txt.groupTuple()
    diamond_hits = DIAMOND_BLASTX.out.txt.groupTuple()
    versions     = versions_ch
}