include { ENTREZDIRECT_ESEARCH       } from "../../../modules/nf-core/entrezdirect/esearch/main"
include { MITOHIFI_FINDMITOREFERENCE } from "../../../modules/nf-core/mitohifi/findmitoreference/main"
include { MITOHIFI_MITOHIFI          } from "../../../modules/nf-core/mitohifi/mitohifi/main"

workflow ASSEMBLE_ORGANELLES {
    take:
    ch_assembly_data // Channel: tuple(meta, assembly / reads )
    ch_assembly_mode // String: enum('c','r')

    main:
    ch_versions = Channel.empty()
    // Check species name is found in NCBI taxonomy database, exit if 0 or >1 found
    ENTREZDIRECT_ESEARCH( ch_assembly_data.map { meta, _assembly_data -> [ meta, meta.sample.name ] }.unique(), 'taxonomy' )

    // Fetch mitochondrial reference sequence
    // TODO: Need to check options to mitohifi modules.
    MITOHIFI_FINDMITOREFERENCE( ch_assembly_data.map { meta, _assembly_data -> [ meta, meta.sample.name ] }.unique(), ENTREZDIRECT_ESEARCH.out.success )

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