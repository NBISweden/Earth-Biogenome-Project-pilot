nextflow.enable.types = true

process OATK {
    tag "${id}"
    label 'process_medium'

    ext [
        args: { oatk_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/oatk:1.0--h577a1d6_1'
        : 'biocontainers/oatk:1.0--h577a1d6_1'}"

    input:
    record(
        id: String,
        meta: Record?,
        reads: Path,
        mito_hmm: List<Path>,
        pltd_hmm: List<Path>,
        oatk_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        mito_fasta: file("*mito.ctg.fasta", optional: true),
        pltd_fasta: file("*pltd.ctg.fasta", optional: true),
        mito_bed: file("*mito.ctg.bed", optional: true),
        pltd_bed: file("*pltd.ctg.bed", optional: true),
        mito_gfa: file("*mito.gfa", optional: true),
        pltd_gfa: file("*pltd.gfa", optional: true),
        annot_mito_txt: file("*annot_mito.txt", optional: true),
        annot_pltd_txt: file("*annot_pltd.txt", optional: true),
        clean_gfa: file("*utg.clean.gfa", optional: true),
        final_gfa: file("*utg.final.gfa", optional: true),
        initial_gfa: file("*utg.gfa", optional: true),
        multiplex_gfa: file("*utg.multiplex.gfa", optional: true),
        unzip_gfa: file("*utg.unzip.gfa", optional: true)
    )

    topic:
    tuple(task.process, 'oatk', eval("oatk --version 2>&1")) >> 'versions'

    script:
    def prefix = task.ext.prefix
    def mito_fam = mito_hmm.find { f -> f.name.endsWith('.fam') }
    def pltd_fam = pltd_hmm.find { f -> f.name.endsWith('.fam') }
    def mito_hmm_arg = mito_fam ? "-m ${mito_fam}" : ''
    def pltd_hmm_arg = pltd_fam ? "-p ${pltd_fam}" : ''
    """
    oatk \\
        ${task.ext.args} \\
        $mito_hmm_arg \\
        $pltd_hmm_arg \\
        -t $task.cpus \\
        -o ${prefix} \\
        $reads
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}.mito.ctg.fasta ${prefix}.mito.ctg.bed ${prefix}.mito.gfa
    touch ${prefix}.annot_mito.txt ${prefix}.utg.final.gfa ${prefix}.utg.clean.gfa
    """
}
