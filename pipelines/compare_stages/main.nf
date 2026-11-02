#!/usr/bin/env nextflow

/*
 * Compare the assemblies of a sample across the stages (raw, decontaminated, purged, polished,
 * scaffolded, curated) with QUAST, in one run per sample, so that the change of the contiguity
 * statistics from stage to stage is in one table and one set of plots.
 *
 * Included by the EBP assembly workflow with the assemblies of all stages:
 *
 *   include { workflow as COMPARE_STAGES } from './pipelines/compare_stages'
 *   comparison = COMPARE_STAGES( record(assemblies: all_stages, reference: params.reference) )
 *
 * On its own, it takes a list of assembly records in a -params-file, see tests/params.yml.
 */

nextflow.enable.types = true

include { COMPARE_ASSEMBLIES } from '../../subworkflows/local/compare_assemblies/main'
include { Assembly           } from '../../types/records'

params {
    assemblies: Channel<Assembly> // The assemblies of all stages of the samples
    reference : Path?             // Reference genome to compare to
}

workflow {

    main:
    comparison = COMPARE_ASSEMBLIES( params.assemblies, params.reference )

    publish:
    quast = comparison.quast
    logs  = comparison.logs
}

output {
    quast: Channel<Record> {
        path { r -> "quast/${r.meta.id}" }
    }
    logs: Channel<Path> {
        path 'logs'
    }
}
