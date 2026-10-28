#!/usr/bin/env nextflow

/*
 * Genome assembly of HiFi reads with hifiasm, organelle assembly with MitoHiFi and Oatk,
 * followed by assembly QC.
 *
 * Takes one wide sample record per sample (see types/records.nf). The HiFi reads are read from
 * `hifi_merged` and the Hi-C read pairs from `hic`. Returns the raw assemblies as wide assembly
 * records that carry the fields of their sample.
 */

nextflow.enable.types = true

include { ASSEMBLE_HIFI                 } from './assemble_hifi'
include { ASSEMBLE_ORGANELLES           } from './assemble_organelles'
include { workflow as ASSEMBLY_QC_RAW   } from '../assembly_qc'
include { Assembly; Sample; KmerParams; BuscoParams } from '../../types/records'

params {
    samples: Channel<Sample>       // Samples to assemble
    provided: Channel<Assembly>?   // Assemblies of this stage that were supplied. Evaluated, not assembled.

    // Assembly
    oatkdb                 : Path?                // Setup: git clone --depth 1 https://github.com/c-zhou/OatkDB.git . Then oatkdb = "/path/to/OatkDB/v20230921"
    nuclear_assembly_mode  : Boolean = true       // values: true, false
    organelle_assembly_mode: String  = 'contigs'  // values: 'reads', 'contigs', 'none'
    hifiasm                : Map<String,String>?  // Expects a Map of options [ opts01: "--opt1 --opt2 val", opts02: "--opt2 val --opt3" ]
    hifiasm_with_hic       : Boolean = false

    // Assembly QC
    busco     : BuscoParams = record(lineages: null, lineages_db_path: null, max_lineages: null)
    fastk     : KmerParams  = record(kmer_size: null)
    meryl     : KmerParams  = record(kmer_size: null)
    reference : Path?
    use_phased: Boolean     = false
}

workflow {

    main:
    if ( params.organelle_assembly_mode == 'contigs' && !params.nuclear_assembly_mode ) {
        error "Organelle assembly in 'contigs' mode requires 'nuclear_assembly_mode = true'"
    }
    def organelles = params.organelle_assembly_mode in ['contigs', 'reads']

    // Nuclear assembly
    nuclear = ASSEMBLE_HIFI(
        params.samples.filter { _s -> params.nuclear_assembly_mode },
        params.hifiasm ?: [:],
        params.hifiasm_with_hic,
        params.use_phased
    )

    // Organelle assembly
    organelle = ASSEMBLE_ORGANELLES(
        params.samples.filter { _s -> organelles },
        nuclear.assemblies.filter { _a -> organelles },
        params.organelle_assembly_mode,
        params.oatkdb
        // TODO: filter organelle contigs from primary assembly
    )

    assemblies = nuclear.assemblies.mix(params.provided ?: channel.empty())

    qc = ASSEMBLY_QC_RAW(
        record(
            assemblies: assemblies,
            busco: params.busco,
            fastk: params.fastk,
            meryl: params.meryl,
            reference: params.reference,
            use_phased: params.use_phased
        )
    )

    publish:
    assemblies     = assemblies
    hifiasm        = nuclear.hifiasm
    gfa2fa         = nuclear.gfa2fa
    mito_reference = organelle.reference
    mitohifi       = organelle.mitohifi
    oatk           = organelle.oatk
    dotplots       = organelle.dotplots
    merquryfk      = qc.merquryfk
    merqury        = qc.merqury
    busco          = qc.busco
    gfastats       = qc.gfastats
    quast          = qc.quast
    logs           = nuclear.logs.mix(organelle.logs).mix(qc.logs)
}

output {
    assemblies: Channel<Assembly> {
        path { a ->
            a.pri_fasta >> ("assemblies/${a.id}/")
            a.alt_fasta >> ("assemblies/${a.id}/")
        }
    }
    hifiasm: Channel<Record> {
        path { r -> "hifiasm/${r.meta.assembly.build}" }
    }
    gfa2fa: Channel<Record> {
        path { r -> "gfa2fa/${r.meta.assembly.build}" }
    }
    mito_reference: Channel<Record> {
        path "organelle/mitohifi/references"
    }
    mitohifi: Channel<Record> {
        path { r -> "organelle/mitohifi/${params.organelle_assembly_mode == 'contigs' ? r.meta.assembly.build : r.meta.id}" }
    }
    oatk: Channel<Record> {
        path "organelle/oatk"
    }
    dotplots: Channel<Record> {
        path "organelle/dotplots"
    }
    merquryfk: Channel<Record> {
        path { r -> "merquryfk/${r.meta.assembly.build}" }
    }
    merqury: Channel<Record> {
        path { r -> "merqury/${r.meta.assembly.build}" }
    }
    busco: Channel<Record> {
        path { r ->
            r.short_summaries_txt >> ("busco/${r.meta.assembly.build}/")
            r.short_summaries_json >> ("busco/${r.meta.assembly.build}/")
        }
    }
    gfastats: Channel<Record> {
        path { r -> "gfastats/${r.meta.assembly.build}" }
    }
    quast: Channel<Record> {
        path 'quast'
    }
    logs: Channel<Path> {
        path 'logs'
    }
}
