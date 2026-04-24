include { softwareVersionsToYAML } from "../../../modules/local/nf-core-utilities.nf"

workflow REPORT_VERSIONS {

    take:
    versions_topic_ch
    versions_heredoc_ch

    main:
    // Code adapted from the nf-core pipeline template: https://github.com/nf-core/tools/tree/main/nf_core/pipeline-template)

    // Branch topic versions by data structure
    def topic_type = versions_topic_ch
        .distinct()
        .branch { entry ->
            file: entry instanceof Path
            tuple: true
        }

    // Process topic channel tuples into strings
    def topic_versions_string = topic_type.tuple
        .map { process, tool, version ->
            [ process, "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    // Mix all versions channels and convert to YAML
    softwareVersionsToYAML(versions_heredoc_ch.mix(topic_type.file))
        .mix(topic_versions_string)
        .collectFile(
            name: 'assembly_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { multiqc_versions_report_ch }

    emit:
    multiqc_versions_report_ch

}