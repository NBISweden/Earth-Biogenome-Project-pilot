include { ENTREZDIRECT_ESEARCH       } from "../../../modules/nf-core/entrezdirect/esearch/main"
include { MITOHIFI_FINDMITOREFERENCE } from "../../../modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI          } from "../../../modules/nf-core/mitohifi/mitohifi/main"

workflow ASSEMBLE_ORGANELLES {
    take:
    ch_assembly_data // Channel: tuple(meta, assembly / reads )
    ch_assembly_mode // String: enum('c','r')

    main:
    ch_versions = Channel.empty()
    // Lookup species name in NCBI taxonomy database, get count of hits
    ENTREZDIRECT_ESEARCH( ch_assembly_data.map { meta, _assembly_data -> [ meta, meta.sample.name ] }.unique(), 'taxonomy' )

    // Branch channel based on the number of hits in NCBI taxonomy (1 hit is expected, 0 or >1 requires an alternative strategy)
    ENTREZDIRECT_ESEARCH.out.count
        .map { meta, species, count -> [ meta, species, count.toInteger() ] }
        .branch { meta, species, count ->
            attempt_mitohifi: count == 1
            attempt_oatk: count != 1
        }
        .set { count_ch }

    // If single hit in taxonomy lookup, attempt MitoHiFi workflow: step one fetch mitochondrial reference sequence
    // TODO: Need to check options to mitohifi modules.
    MITOHIFI_FINDMITOREFERENCE( count_ch.attempt_mitohifi.map { meta, species, _count -> [ meta, species ] } )

    // Prepare mitohifi_ch based on input data type (mode c: contigs, mode r: reads)
    mitohifi_ch = ch_assembly_data
        .combine(
            MITOHIFI_FINDMITOREFERENCE.out.fasta
                .join(MITOHIFI_FINDMITOREFERENCE.out.gb),
            by: 0
        )
        .multiMap { meta, assembly_data, mitofa, mitogb ->
            input: [ meta, ch_assembly_mode == "c" ? assembly_data.pri_fasta : assembly_data ] // If c, use contigs. Else r, use reads.
            reference: mitofa
            genbank: mitogb
            mito_code: meta.sample.mito_code
        }

    // Run mitochondrial assembly
    MITOHIFI_MITOHIFI(
        mitohifi_ch.input,
        mitohifi_ch.reference,
        mitohifi_ch.genbank,
        ch_assembly_mode,
        mitohifi_ch.mito_code,
    )
    ch_versions = ch_versions.mix(
        ENTREZDIRECT_ESEARCH.out.versions.first(),
        MITOHIFI_FINDMITOREFERENCE.out.versions.first(),
        MITOHIFI_MITOHIFI.out.versions.first()
    )

    emit:
    // TODO: emit filtered assembly, or contig list to purge. Also emit from 03_assemble.
    versions = ch_versions
}