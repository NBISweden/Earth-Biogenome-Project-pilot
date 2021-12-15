#! /usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

nextflow.enable.dsl = 2

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
    input.assembly_ch
        .filter { !it.isEmpty() }
        .transpose()     // Data is [ sample, [id:'assemblerX_build1', path:'/path/to/assembly']]
        .map { sample, assembly -> [ sample, [ id: assembly.id, path: file( assembly.path, checkIfExists: true ) ] ] }
        .set { assembly_ch }

    emit:
    assembly = assembly_ch
    hic      = input.hic_ch.filter { !it.isEmpty() }
    hifi     = input.hifi_ch.filter { !it.isEmpty() }
    rnaseq   = input.rnaseq_ch.filter { !it.isEmpty() }
    isoseq   = input.isoseq_ch.filter { !it.isEmpty() }

}

def readYAML( yamlfile ) {
    // TODO: Validate sample file
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}
