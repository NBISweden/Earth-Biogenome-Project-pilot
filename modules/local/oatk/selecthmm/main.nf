process OATK_SELECTHMM {
    tag "${species}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/oatk:1.0--h577a1d6_1'
        : 'biocontainers/oatk:1.0--h577a1d6_1'}"

    input:
    tuple val(meta), val(species), val(lineage)
    val oatkdb

    output:
    tuple path("mito/*_mito.fam"), path("mito/*_mito.fam.h3f"), path("mito/*_mito.fam.h3i"), path("mito/*_mito.fam.h3m"), path("mito/*_mito.fam.h3p"), emit: mito_hmm, optional: true
    tuple path("pltd/*_pltd.fam"), path("pltd/*_pltd.fam.h3f"), path("pltd/*_pltd.fam.h3i"), path("pltd/*_pltd.fam.h3m"), path("pltd/*_pltd.fam.h3p"), emit: pltd_hmm, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Search for most specific lineage match
    taxa=""
    echo "${lineage.toLowerCase()}" \
        | tr -s " ;" "\\n" \
        > lineages.txt
    while IFS= read -r lineage; do
        if grep -q "\$lineage" ${oatkdb}/TAXID; then
            taxa="\$lineage"
        fi
    done < lineages.txt
    echo "Selected taxa: \$taxa"

    if [ -n "\$taxa" ]; then
        mkdir mito pltd
        find -L ${oatkdb} -type f -name "\${taxa}_mito.*" -exec cp "{}" mito \\;
        find -L ${oatkdb} -type f -name "\${taxa}_pltd.*" -exec cp "{}" pltd \\;
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oatk: \$(oatk --version)
    END_VERSIONS
    """
}
