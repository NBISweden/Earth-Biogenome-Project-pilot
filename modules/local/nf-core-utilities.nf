// These functions were adopted from the nf-core pipeline template: https://github.com/nf-core/tools/tree/main/nf_core/pipeline-template

// Get channel of software versions used in pipeline in YAML format
def softwareVersionsToYAML(ch_versions) {
    return ch_versions.unique().map { version -> processVersionsFromYAML(version) }.unique().mix(channel.of(workflowVersionToYAML()))
}

// Get software versions for pipeline
def processVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions = yaml.load(yaml_file)
    return yaml.dumpAsMap(versions).trim()
}

// Get workflow version for pipeline
def workflowVersionToYAML() {
    return """
    Workflow:
      ${workflow.manifest.name}: ${getWorkflowVersion()}
      Nextflow: ${workflow.nextflow.version}
    """.stripIndent().trim()
}

// Generate workflow version string
def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    return version_string
}
