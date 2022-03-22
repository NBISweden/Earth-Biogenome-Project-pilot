#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

import org.yaml.snakeyaml.Yaml

include { SAMTOOLS_FASTA } from "../../modules/local/samtools/fasta/main"

workflow {
    PREPARE_INPUT( params.input )
}

/* params.input example sample sheet (samplesheet.yml)
```yaml
sample:
  id: ERGA_Artic_Fox
assembly:
  - id: assemblerX_build1
    path: /path/to/assembly
  - id: assemblerY_build1
    path: /path/to/assembly
hic:
  - /path/to/reads
hifi:
  - /path/to/reads
rnaseq:
  - /path/to/reads
isoseq:
  - /path/to/reads
```
leads to the following YAML data structure
[
    sample:[
        id:ERGA_Artic_Fox
    ],
    assembly:[
        [
            id:assemblerX_build1,
            path:/path/to/assembly
        ],
        [
            id:assemblerY_build1,
            path:/path/to/assembly
        ]
    ],
    hic:[
        /path/to/reads
    ],
    hifi:[
        /path/to/reads
    ],
    rnaseq:[
        /path/to/reads
    ],
    isoseq:[
        /path/to/reads
    ]
]
*/

workflow PREPARE_INPUT {

    take:
    infile

    main:
    // Read in YAML
    Channel.fromPath( infile )
        .map { file -> readYAML( file ) }
        .multiMap { data ->
            assembly_ch : data.assembly ? [ data.sample, data.assembly ] : []
            hic_ch      : data.hic      ? [ data.sample, data.hic.collect { file( it, checkIfExists: true ) } ] : []
            hifi_ch     : data.hifi     ? [ data.sample, data.hifi.collect { file( it, checkIfExists: true ) } ] : []
            rnaseq_ch   : data.rnaseq   ? [ data.sample, data.rnaseq.collect { file( it, checkIfExists: true ) } ] : []
            isoseq_ch   : data.isoseq   ? [ data.sample, data.isoseq.collect { file( it, checkIfExists: true ) } ] : []
        }
        .set{ input }

    // Convert assembly filename to files for correct staging
    input.assembly_ch
        .filter { !it.isEmpty() }
        .transpose()     // Data is [ sample, [id:'assemblerX_build1', path:'/path/to/assembly']]
        .map { sample, assembly -> [ sample, [ id: assembly.id, path: file( assembly.path, checkIfExists: true ) ] ] }
        .set { assembly_ch }

    // Convert HiFi BAMS to FastQ
    input.hifi_ch
        .filter { !it.isEmpty() }
        .transpose()   // Transform to [ [ id: 'sample_name'], file('/path/to/read')  ]
        .branch { meta, filename ->
            bam_ch: filename.toString().endsWith(".bam")
            fastx_ch: true // assume everything else is fastx
        }.set { hifi }
    SAMTOOLS_FASTA ( hifi.bam_ch )
    hifi.fastx_ch.mix( SAMTOOLS_FASTA.out.fasta )
        .groupTuple()
        .set { hifi_fastx_ch }


    emit:
    assemblies = assembly_ch
    hic        = input.hic_ch.filter { !it.isEmpty() }
    hifi       = hifi_fastx_ch
    rnaseq     = input.rnaseq_ch.filter { !it.isEmpty() }
    isoseq     = input.isoseq_ch.filter { !it.isEmpty() }

}

def readYAML( yamlfile ) {
    // TODO: Validate sample file
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}
