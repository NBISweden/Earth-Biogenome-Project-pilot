nextflow.enable.types = true

include { ASSEMBLE_HIFI       } from "./assemble_hifi"
include { ASSEMBLE_ORGANELLES } from "./assemble_organelles"

record Reads {
    id   : String
    meta : Record
    reads: List<Path>
}

record Assembly {
    id       : String  // meta.assembly.build
    meta     : Record  // sample meta data, including meta.assembly
    pri_fasta: Path
    alt_fasta: Path?
    pri_gfa  : Path?
    alt_gfa  : Path?
}

workflow ASSEMBLE {
    take:
    hifi_reads              : Channel<Reads> // One merged reads file per sample
    hic_reads               : Channel<Reads>
    oatkdb                  : Path?          // /path/to/oatkdb
    nuclear_assembly_mode   : Boolean        // assemble chromosomes
    organelle_assembly_mode : String         // enum: contigs, reads, or none
    hifiasm_configs         : Map<String,String>
    hifiasm_with_hic        : Boolean
    use_phased              : Boolean

    main:
    if ( organelle_assembly_mode == 'contigs' && !nuclear_assembly_mode ) {
        error "Organelle assembly in 'contigs' mode requires 'nuclear_assembly_mode = true'"
    }
    organelles = organelle_assembly_mode in ['contigs', 'reads']

    // Nuclear assembly
    nuclear = ASSEMBLE_HIFI(
        hifi_reads.filter { _r -> nuclear_assembly_mode },
        hic_reads,
        hifiasm_configs,
        hifiasm_with_hic,
        use_phased
    )

    // Organelle assembly
    organelle = ASSEMBLE_ORGANELLES(
        hifi_reads.filter { _r -> organelles },        // mitohifi + oatk: reads
        nuclear.assemblies.filter { _a -> organelles }, // mitohifi + oatk: assemblies
        organelle_assembly_mode,                        // mitohifi + oatk: mode
        oatkdb                                          // oatk: hmm database
        // TODO: filter organelle contigs from primary assembly
    )

    emit:
    raw_assemblies  : Channel<Assembly> = nuclear.assemblies
    logs            : Channel<Path>     = nuclear.logs.mix(organelle.logs)
    hifiasm         : Channel<Record>   = nuclear.hifiasm
    gfa2fa          : Channel<Record>   = nuclear.gfa2fa
    mito_reference  : Channel<Record>   = organelle.reference
    mitohifi        : Channel<Record>   = organelle.mitohifi
    oatk            : Channel<Record>   = organelle.oatk
    dotplots        : Channel<Record>   = organelle.dotplots
}
