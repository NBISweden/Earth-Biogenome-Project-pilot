// Portions of this code adapted from Sanger Institute's "curationpretext" (https://github.com/sanger-tol/curationpretext)
// and from https://github.com/MartinPippel/DAMARVEL/blob/master-v1/scripts/createHiCPlans.sh
process BAM2BED_SORT {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/bedtools_samtools_pip_awk:61be4cfccef27592"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.pairs") , emit: pairs
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: '1'
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools view -@$task.cpus $args ${bam} | \\
    bamToBed | \\
    sort -k4 --parallel=$task.cpus $args2 | \\
    paste -d '\t' - - | \\
    awk -v q=$args3 'BEGIN {FS=\"\t\"; OFS=\"\t\"} { if(int(\$5) >= int(q) && int(\$11) >= int(q)) { if (\$1 > \$7) { print substr(\$4,1,length(\$4)-2),\$7,\$8,\$1,\$2,\$12,\$6,\"UU\"} else { print substr(\$4,1,length(\$4)-2),\$1,\$2,\$7,\$8,\$6,\$12,\"UU\"} } }' | \\
    sort -k2,2V -k4,4V -k3,3n -k5,5n --parallel=$task.cpus $args2 > ${prefix}.pairs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed '1!d; s/.* //')
        bedtools: \$(bedtools --version | sed 's/.*v//')
        sort: \$(sort --version | sed '1!d; s/.* //')
        awk: \$(awk --version | sed '1!d; s/mawk //; s/ .*//')
        paste: \$(paste --version | sed '1!d; s/.* //')
    END_VERSIONS
    """
}
