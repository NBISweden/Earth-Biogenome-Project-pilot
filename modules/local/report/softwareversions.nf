process REPORT_SOFTWAREVERSIONS {
    tag ""
    label 'process_single'

    input:
    val versions_list // List of filenames containing versions

    output:
    path "versions_mqc.yml", emit: yml

    when:
    task.ext.when == null || task.ext.when

    exec:
    def versions_yml = versions_list.collect { vf ->
        vf.text
    }.join('\n')
    file("$task.workDir/versions_mqc.yml").text = versions_yml
}
