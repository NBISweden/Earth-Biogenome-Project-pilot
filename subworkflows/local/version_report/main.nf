nextflow.enable.types = true

include { WRITE_FILE } from "../../../modules/local/report/write_file"

// Tool versions of the workflow as YAML, grouped by process
def versionsYaml(versions: List<List<String>>, workflow_version: String) -> String {
    def by_process = versions
        .collect { v -> v[0] }
        .toSet()
        .toList()
        .toSorted()
        .collect { process ->
            def tool_versions = versions
                .findAll { v -> v[0] == process }
                .collect { v -> "  ${v[1]}: ${v[2]}" }
                .toSet()
                .toList()
                .toSorted()
            "${process}:\n${tool_versions.join('\n')}"
        }
    return "${workflow_version}\n${by_process.join('\n')}"
}

workflow REPORT_VERSIONS {

    take:
    versions        : Channel<Tuple>   // ( process, tool, version ) from the versions topic
    workflow_version: String                               // YAML text with the workflow and Nextflow versions

    main:
    // Code adapted from the nf-core pipeline template: https://github.com/nf-core/tools/tree/main/nf_core/pipeline-template)
    yaml = versions
        .map { t -> [ t[0] as String, t[1] as String, t[2] as String ] }
        .collect()
        .map { all -> versionsYaml(all.toList(), workflow_version) }
    versions_file = WRITE_FILE(
        yaml.map { y -> record(id: 'versions', filename: 'assembly_mqc_versions.yml', content: y) }
    )

    emit:
    report: Value<Path> = versions_file.map { f -> f.file }
}
