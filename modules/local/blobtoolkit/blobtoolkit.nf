process BLOBTOOLKIT {

    label 'process_high'

    conda "${params.enable_conda ? 'https://github.com/blobtoolkit/pipeline/blob/master/env.yaml' : '' }"
    container 'genomehubs/blobtoolkit:2.6.5'

    input:
    tuple val(meta), path(reads, stageAs: 'data/*'), path(assembly, stageAs: 'data/*')
    val busco_lineages
    // path blobtoolkit_config
    path busco_lineage_path, stageAs: 'databases/busco'
    path uniprot_db        , stageAs: 'databases/uniprot_db'
    path ncbi_nt_db        , stageAs: 'databases/ncbi_db'
    path ncbi_taxonomy     , stageAs: 'databases/ncbi_taxdump'

    output:
    tuple val(sample), path("*.snakemake.stats"), emit: snakemake_stats

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: ''
    """
    cat << EOF > ${prefix}.yaml
    assembly:
      accession: draft
      level: scaffold
      prefix: $prefix
    busco:
      lineages:
    ${busco_lineages.tokenize(',').collect{ "    - $it" }.join("\n")}
      lineage_dir: databases/busco
    reads:
      single:
    ${reads.collect{ "    -\n      - ${it.baseName}\n      - PACBIO_SMRT" }.join("\n")}
    settings:
      blobtools2_path: /blobtoolkit/blobtools2
      taxonomy: databases/ncbi_taxdump
      tmp: /tmp
      blast_chunk: 100000
      blast_max_chunks: 10
      blast_overlap: 500
      chunk: 1000000
    similarity:
      defaults:
        evalue: 1e-25
        max_target_seqs: 10
        root: 1
        mask_ids:
          - 7215
      databases:
        -
          local: databases/ncbi_db
          name: nt
          source: ncbi
          tool: blast
          type: nucl
        -
          local: databases/uniprot_db
          max_target_seqs: 1
          name: reference_proteomes
          source: uniprot
          tool: diamond
          type: prot
      taxrule: bestsumorder
    keep_intermediates: true
    EOF
    snakemake -p \\
         --use-conda \\
         --conda-prefix .conda \\
         --directory data \\
         --configfile ${prefix}.yaml \\
         --stats ${prefix}.snakemake.stats \\
         -j $task.cpus \\
         -s /blobtoolkit/insdc-pipeline/Snakefile \\
         --resources btk=1
    """
}