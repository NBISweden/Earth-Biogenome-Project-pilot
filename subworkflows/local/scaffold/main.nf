include { SCAFFOLD_OMNIC_BWA      } from "./scaffold_omnic_bwa"
include { SCAFFOLD_ARIMA_BWA      } from "./scaffold_arima_bwa"

workflow SCAFFOLD {
    take:
    ch_assemblies // [ meta, assembly ]
    ch_hic        // [ meta, hic-pairs ]

    main:

    ch_versions              = Channel.empty()
    ch_scaffolded_assemblies = Channel.empty()
    ch_logs                  = Channel.empty()

    if (params.hic_pipeline == 'omni-c')
    {
        SCAFFOLD_OMNIC_BWA(
            ch_assemblies,
            ch_hic
        )
        ch_versions = ch_versions.mix(SCAFFOLD_OMNIC_BWA.out.versions)
        ch_scaffolded_assemblies = ch_scaffolded_assemblies.mix(SCAFFOLD_OMNIC_BWA.out.assemblies)
        ch_logs = ch_logs.mix(SCAFFOLD_OMNIC_BWA.out.logs)
    }
    else if (params.hic_pipeline == 'arima')
    {
        SCAFFOLD_ARIMA_BWA(
            ch_assemblies,
            ch_hic
        )
        ch_versions = ch_versions.mix(SCAFFOLD_ARIMA_BWA.out.versions)
        ch_scaffolded_assemblies = ch_scaffolded_assemblies.mix(SCAFFOLD_ARIMA_BWA.out.assemblies)
        ch_logs = ch_logs.mix(SCAFFOLD_ARIMA_BWA.out.logs)
    }

    emit:
    assemblies = ch_scaffolded_assemblies
    logs       = ch_logs
    versions   = ch_versions
}
