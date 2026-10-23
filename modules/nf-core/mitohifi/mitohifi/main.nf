nextflow.enable.types = true

process MITOHIFI_MITOHIFI {
    tag "${id}"
    label 'process_high'

    ext [
        args: { mitohifi_args ?: '' },
        prefix: { output_prefix ?: id }
    ]

    // Docker image available at the project github repository
    container 'ghcr.io/marcelauliano/mitohifi:3.2.3'

    input:
    record(
        id: String,
        meta: Record?,
        input: Path,
        ref_fa: Path,
        ref_gb: Path,
        input_mode: String,
        mito_code: String,
        mitohifi_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        fasta: file("*final_mitogenome.fasta", optional: true),
        all_candidates_fa: files("potential_*/*/*.rotated.fa", optional: true),
        stats: file("*contigs_stats.tsv", optional: true),
        gb: file("*final_mitogenome.gb", optional: true),
        gff: file("*final_mitogenome.gff", optional: true),
        all_potential_contigs: file("*all_potential_contigs.fa", optional: true),
        contigs_annotations: file("*contigs_annotations.png", optional: true),
        contigs_circularization: file("contigs_circularization/", optional: true),
        contigs_filtering: file("contigs_filtering/", optional: true),
        coverage_mapping: file("coverage_mapping/", optional: true),
        coverage_plot: file("*coverage_plot.png", optional: true),
        final_mitogenome_annotation: file("*final_mitogenome.annotation.png", optional: true),
        final_mitogenome_choice: file("final_mitogenome_choice/", optional: true),
        final_mitogenome_coverage: file("*final_mitogenome.coverage.png", optional: true),
        potential_contigs: file("potential_contigs/", optional: true),
        reads_mapping_and_assembly: file("*reads_mapping_and_assembly/", optional: true),
        shared_genes: file("*shared_genes.tsv", optional: true)
    )

    topic:
    // WARN: Incorrect version information is provided by tool on CLI. Please update this string when bumping container versions.
    // old version command: \$(mitohifi.py -v | sed 's/.* //')
    tuple(task.process, 'mitohifi', eval('echo 3.2.3')) >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Error: MitoHiFi module does not support Conda. Please use Docker / Singularity instead."
    }
    if (! ["c", "r"].contains(input_mode)) {
        error "r for reads or c for contigs must be specified"
    }
    def prefix = task.ext.prefix
    def zipped = input.name.endsWith('.gz')
    // Mitohifi deletes the original file when renaming headers.
    def fasta = ( zipped ? input.name.replaceAll(/\.gz$/, '') : "${prefix}.${input.name}" )
    """
    ${zipped ? "gzip -dc $input >" : "cp ${input}"} ${fasta}

    # override workflow exit on mitohifi.py error (issues: #277, #86, #145)
    set +e
    mitohifi.py -${input_mode} ${fasta} \\
        -f ${ref_fa} \\
        -g ${ref_gb} \\
        -o ${mito_code} \\
        -t $task.cpus ${task.ext.args}
    set -e

    rm -f ${fasta} fixed_header_contigs.fasta

    # Rename files to include prefix
    find . -maxdepth 1 -type f ! -name '.*' -exec sh -c 'for f do mv "\$f" "${prefix}.\${f#./}"; done' sh {} +
    # Remove broken link that disrupts the publish operation
    find . -xtype l -delete

    # Test for mitohifi reference files:
    # *mitogenome.fasta && *contigs_stats.tsv: Mitohifi complete output, proceed with workflow
    # ! *mitogenome.fasta && ! *contigs_stats.tsv: Mitohifi did not complete, proceed with workflow
    # *mitogenome.fasta && ! *contigs_stats.tsv: Mitohifi partial output: program crashed, exit workflow

    FASTA=\$(find . -maxdepth 1 -name "*mitogenome.fasta" -type f | wc -l)
    STATS=\$(find . -maxdepth 1 -name "*contigs_stats.tsv" -type f | wc -l)

    if [[ \$FASTA -ne \$STATS ]]; then
        exit 1
    fi
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}.final_mitogenome.fasta
    touch ${prefix}.contigs_stats.tsv
    """
}
