nextflow.enable.types = true

process OATK_SELECTHMM {
    tag "${species}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/oatk:1.0--h577a1d6_1'
        : 'biocontainers/oatk:1.0--h577a1d6_1'}"

    input:
    record(
        id: String,
        meta: Record?,
        species: String,
        lineage: String,
        oatkdb: Path
    )

    output:
    record(
        id: id,
        meta: meta,
        mito_hmm: files("mito/*_mito.fam*", optional: true),
        pltd_hmm: files("pltd/*_pltd.fam*", optional: true)
    )

    topic:
    tuple(task.process, 'oatk', eval("oatk --version")) >> 'versions'

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
    """

    stub:
    """
    mkdir mito pltd
    for EXT in fam fam.h3f fam.h3i fam.h3m fam.h3p; do
        touch mito/stub_mito.\$EXT pltd/stub_pltd.\$EXT
    done
    """
}
