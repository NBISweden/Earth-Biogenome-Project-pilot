process REPORT_SOFTWAREVERSIONS {
    tag ""
    label 'process_single'

    input:
    val versions_list

    output:
    path "assembly_mqc_versions.yml", emit: yml

    when:
    task.ext.when == null || task.ext.when

    exec:
    def versions_yml = versions_list
        .collect { vf ->
            vf.text
        }
        .join('')
    file("${task.workDir}/assembly_mqc_versions.yml").text = versions_yml
}
