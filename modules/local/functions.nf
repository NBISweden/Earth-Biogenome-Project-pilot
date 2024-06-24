// A set of Nextflow channel helper functions

/* A function to enhance the channel operator `combine`,
to use specific keys of the meta map as keys for `combine`.

Example:

MYPROCESS (
    combineByMetaKeys (
        channel_one,    // [ [ id: 'sample', coverage: X ], bam ]
        channel_two,    // [ [ id: 'sample', group: Y    ], bed ]
        keySet: ['id'], // Use operator `combine` with `by` key meta.id
        meta: 'rhs'     // Output channel retains meta map from channel_two
    ) // Output: [ [ id: 'sample', group: Y ], bam, bed ]
)
*/
def combineByMetaKeys( Map args = [:], lhs, rhs ){
    assert args.keySet != null
    def combine_args = [ by: 0 ]
    switch(args.meta) {
        case 'lhs':
            // Return meta map from lhs
            return lhs.map{ [ it[0].subMap(args.keySet) ] + it }
                .combine( combine_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
                .map { it.tail() }
        case 'rhs':
            // Return meta map from rhs
            def restructure = { key, clip, ... tail -> def list = tail.toList(); [ list.removeAt(clip) ] + list }
            return lhs.map{ [ it[0].subMap(args.keySet), it.size()-1 ] + it.tail() }
                .combine( combine_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
                .map { restructure(*it) }
        case 'merge':
            // Return merged meta map
            def restructure = { key, clip, ... tail -> def list = tail.toList(); [ list.head().deepMerge(list.removeAt(clip)) ] + list.tail() }
            return lhs.map{ [ it[0].subMap(args.keySet), it.size() ] + it }
                .combine( combine_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
                .map { restructure(*it) }
        case 'subset':
            // Return meta keys subset
            return lhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() }
                .combine( combine_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
        default:
            assert args.meta in ['lhs','rhs','merge','subset']
    }
}

/* A function to enhance the channel operator `join`,
to use specific keys of the meta map as keys for `join`.

Example:

MYPROCESS (
    joinByMetaKeys (
        channel_one,    // [ [ id: 'sample', coverage: X ], bam ]
        channel_two,    // [ [ id: 'sample', group: Y    ], bed ]
        keySet: ['id'], // Use operator `join` with `by` key meta.id
        meta: 'rhs',    // Output channel retains meta map from channel_two
        remainder: true // additional flags valid for the channel operator can be set here
    ) // Output: [ [ id: 'sample', group: Y ], bam, bed ]
)
*/
def joinByMetaKeys( Map args = [:], lhs, rhs ) {
    assert args.keySet != null
    def join_args = args.subMap(['remainder', 'failOnMismatch', 'failOnDuplicate']) + [ by: 0 ]
    switch(args.meta) {
        case 'lhs':
            // Return meta map from lhs
            return lhs.map{ [ it[0].subMap(args.keySet) ] + it }
                .join( join_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
                .map { it.tail() }
        case 'rhs':
            // Return meta map from rhs
            def restructure = { key, clip, ... tail -> def list = tail.toList(); [ list.removeAt(clip) ] + list }
            return lhs.map{ [ it[0].subMap(args.keySet), it.size()-1 ] + it.tail() }
                .join( join_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
                .map { restructure(*it) }
        case 'merge':
            // Return merged meta map
            def restructure = { key, clip, ... tail -> def list = tail.toList(); [ list.head().deepMerge(list.removeAt(clip)) ] + list.tail() }
            return lhs.map{ [ it[0].subMap(args.keySet), it.size() ] + it }
                .join( join_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
                .map { restructure(*it) }
        case 'subset':
            // Return meta keys subset
            return lhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() }
                .join( join_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
        default:
            assert args.meta in ['lhs','rhs','merge','subset']
    }
}

/* Assemblies can be supplied as input.
This function filters assemblies for a particular stage of analysis.
*/
def assembliesFromStage( assemblies, String stage ) {
    assemblies.filter { meta, assembly -> meta.assembly?.stage == stage }
}

/* Assemblies need to have meta data updated to use in filenames.
This function updates the meta data so meta.assembly.build, etc can
be used for `ext.prefix` for example.
*/
def setAssemblyStage( assemblies, String stage ) {
    assemblies.map{ meta, assembly -> [ meta.deepMerge([ assembly: [ stage: stage, build: "${meta.assembly.assembler}-${stage}-${meta.assembly.id}" ] ]), assembly ] }
}

/* Assemblies are stored in a Map object.
This function extracts the primary assembly for analysis.
*/
def getPrimaryAssembly( assemblies ) {
    assemblies.map { meta, assembly -> [ meta, assembly.pri_fasta ] }
}

/* Assemblies are stored in a Map object.
This function extracts the primary and alternate assembly for analysis.
The channel operator `transpose` can be used on the output to emit each
fasta separately.
*/
def getEachAssembly( assemblies ) {
    assemblies.map { meta, assembly ->
        [ meta, ( assembly.alt_fasta ? [ assembly.pri_fasta, assembly.alt_fasta ] : assembly.pri_fasta ) ]
    }
}

/* Assemblies are stored in a Map object.
This function constructs that Map object from a mixed input channel,
relying on the name of assembly to denote which is primary and which is alternate.

Example:

PROCESSA.out.fasta // contains [ meta, 'assembler-stage-uuid-hap0.fasta' ] and [ meta, 'assembler-stage-uuid-hap1.fasta']
    .set { processed_assemblies }
output_assembly_ch = constructAssemblyRecord ( processed_assemblies )

The output is the map:
    [
        meta, // meta map
        [     // assembly map
            assembler: 'assembler',
            stage:     'stage',
            id:        'uuid',
            build:     'assembler-stage-uuid',
            pri_fasta: 'assembler-stage-uuid-hap0.fasta',
            alt_fasta: 'assembler-stage-uuid-hap1.fasta'
        ]
    ]
*/
def constructAssemblyRecord( assemblies ) {
    assemblies.groupTuple( sort: { a, b -> a.name <=> b.name } )
        .map { meta, fasta ->
            def asm_meta = meta.assembly.subMap(['assembler','stage','id','build'])
            [ meta, asm_meta + (fasta.size() == 1 ? [ pri_fasta: fasta[0] ] : [ pri_fasta: fasta[0], alt_fasta: fasta[1] ] ) ]
        }
}