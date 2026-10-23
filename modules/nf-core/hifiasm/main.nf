nextflow.enable.types = true

process HIFIASM {
    tag "${id}"
    label 'process_high'

    ext [
        args: { hifiasm_args ?: '' },
        prefix: { output_prefix ?: "${id}" }
    ]

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hifiasm:0.25.0--h5ca1c30_0' :
        'biocontainers/hifiasm:0.25.0--h5ca1c30_0' }"

    input:
    record(
        id: String,
        meta: Record?,
        long_reads: List<Path>,
        ul_reads: List<Path>,
        paternal_kmer_dump: Path?,
        maternal_kmer_dump: Path?,
        hic_read1: List<Path>,
        hic_read2: List<Path>,
        bin_files: List<Path>,
        hifiasm_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        raw_unitigs: file("*.r_utg.gfa"),
        bin_files: files("*.bin", optional: true),
        processed_unitigs: file("*.p_utg.gfa", optional: true),
        primary_contigs: file("${task.ext.prefix}.{p_ctg,bp.p_ctg,hic.p_ctg}.gfa", optional: true),
        alternate_contigs: file("${task.ext.prefix}.{a_ctg,hic.a_ctg}.gfa", optional: true),
        hap1_contigs: file("${task.ext.prefix}.*.hap1.p_ctg.gfa", optional: true),
        hap2_contigs: file("${task.ext.prefix}.*.hap2.p_ctg.gfa", optional: true),
        corrected_reads: file("*.ec.fa.gz", optional: true),
        read_overlaps: file("*.ovlp.paf.gz", optional: true),
        log: file("${task.ext.prefix}.stderr.log")
    )

    topic:
    tuple(task.process, 'hifiasm', eval('hifiasm --version 2>&1')) >> 'versions'

    script:
    def prefix = task.ext.prefix

    def long_reads_sorted = long_reads.toSorted { a, b -> a.name <=> b.name }.join(' ')
    def ultralong = ul_reads ? "--ul ${ul_reads.toSorted { a, b -> a.name <=> b.name }.join(',')}" : ""

    if (paternal_kmer_dump && hic_read1) {
        error "hifiasm trio binning mode and Hi-C phasing can not be used at the same time."
    }

    def input_trio = ""
    if (paternal_kmer_dump || maternal_kmer_dump) {
        if (!(paternal_kmer_dump && maternal_kmer_dump)) {
            error "Either the maternal or paternal kmer dump is missing!"
        }
        input_trio = "-1 ${paternal_kmer_dump} -2 ${maternal_kmer_dump}"
    }

    def input_hic = ""
    if (hic_read1 || hic_read2) {
        if (!(hic_read1 && hic_read2)) {
            error "Either the forward or reverse Hi-C reads are missing!"
        }
        input_hic = "--h1 ${hic_read1.join(",")} --h2 ${hic_read2.join(",")}"
    }
    """
    hifiasm \\
        ${task.ext.args} \\
        -t ${task.cpus} \\
        ${input_trio} \\
        ${input_hic} \\
        ${ultralong} \\
        -o ${prefix} \\
        ${long_reads_sorted} \\
        2>| >( tee ${prefix}.stderr.log >&2 )

    if [ -f ${prefix}.ec.fa ]; then
        gzip ${prefix}.ec.fa
    fi

    if [ -f ${prefix}.ovlp.paf ]; then
        gzip ${prefix}.ovlp.paf
    fi
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}.r_utg.gfa
    touch ${prefix}.ec.bin
    touch ${prefix}.ovlp.source.bin
    touch ${prefix}.ovlp.reverse.bin
    touch ${prefix}.hic.tlb.bin
    touch ${prefix}.hic.lk.bin
    touch ${prefix}.bp.p_ctg.gfa
    touch ${prefix}.p_utg.gfa
    touch ${prefix}.a_ctg.gfa
    touch ${prefix}.bp.hap1.p_ctg.gfa
    touch ${prefix}.bp.hap2.p_ctg.gfa
    echo "" | gzip > ${prefix}.ec.fa.gz
    echo "" | gzip > ${prefix}.ovlp.paf.gz
    touch ${prefix}.stderr.log
    """
}
