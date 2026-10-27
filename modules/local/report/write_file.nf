nextflow.enable.types = true

// Write text to a file, for inputs that are built in the workflow
process WRITE_FILE {
    tag "${filename}"
    label 'process_single'

    input:
    record(
        id: String,
        filename: String,
        content: String
    )

    output:
    record(
        id: id,
        file: file("${filename}")
    )

    exec:
    file("${task.workDir}/${filename}").text = content
}
