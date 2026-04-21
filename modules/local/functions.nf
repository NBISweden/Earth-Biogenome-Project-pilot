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
    if ( args.meta == 'lhs') {
        // Return meta map from lhs
        return lhs.map{ [ it[0].subMap(args.keySet) ] + it }
            .combine( combine_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
            .map { it.tail() }
    } else if ( args.meta == 'rhs' ) {
        // Return meta map from rhs
        def restructure = { _key, clip, tail -> [ tail.removeAt(clip) ] + tail }
        return lhs.map{ [ it[0].subMap(args.keySet), it.size()-1 ] + it.tail() }
            .combine( combine_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
            .map { tuple -> restructure.call(tuple[0], tuple[1], tuple[2..-1]) }
    } else if ( args.meta == 'merge') {
        // Return merged meta map
        def restructure = { _key, clip, tail -> [ deepMergeMaps(tail.head(), tail.removeAt(clip)) ] + tail.tail() }
        return lhs.map{ [ it[0].subMap(args.keySet), it.size() ] + it }
            .combine( combine_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
            .map { tuple -> restructure.call(tuple[0], tuple[1], tuple[2..-1]) }
    } else if ( args.meta == 'subset') {
        // Return meta keys subset
        return lhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() }
            .combine( combine_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
    } else {
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
    if ( args.meta == 'lhs' ) {
        // Return meta map from lhs
        return lhs.map{ [ it[0].subMap(args.keySet) ] + it }
            .join( join_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
            .map { it.tail() }
    } else if ( args.meta == 'rhs' ) {
        // Return meta map from rhs
        def restructure = { _key, clip, tail -> [ tail.removeAt(clip) ] + tail }
        return lhs.map{ [ it[0].subMap(args.keySet), it.size()-1 ] + it.tail() }
            .join( join_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
            .map { tuple -> restructure.call(tuple[0], tuple[1], tuple[2..-1]) }
    } else if ( args.meta == 'merge' ){
        // Return merged meta map
        def restructure = { _key, clip, tail -> [ deepMergeMaps(tail.head(), tail.removeAt(clip)) ] + tail.tail() }
        return lhs.map{ [ it[0].subMap(args.keySet), it.size() ] + it }
            .join( join_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
            .map { tuple -> restructure.call(tuple[0], tuple[1], tuple[2..-1]) }
    } else if ( args.meta == 'subset' ) {
        // Return meta keys subset
        return lhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() }
            .join( join_args, rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
    } else {
        assert args.meta in ['lhs','rhs','merge','subset']
    }
}

/* Assemblies can be supplied as input.
This function filters assemblies for a particular stage of analysis.
*/
def assembliesFromStage( assemblies, String stage ) {
    assemblies.filter { meta, _assembly -> meta.assembly?.stage == stage }
}

/* Assemblies need to have meta data updated to use in filenames.
This function updates the meta data so meta.assembly.build, etc can
be used for `ext.prefix` for example.
*/
def setAssemblyStage( assemblies, String stage ) {
    assemblies.map{ meta, assembly -> [ deepMergeMaps(meta,[ assembly: [ stage: stage, build: "${meta.assembly.assembler}-${stage}-${meta.assembly.id}" ] ]), assembly ] }
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
output_assembly_ch = constructAssemblyRecord ( processed_assemblies, sort_by_name )

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
def constructAssemblyRecord( assemblies, Boolean byName ) {
    // assemblies: Channel: [ meta, fasta ]
    assemblies.groupTuple( sort: byName? { a, b -> a.name <=> b.name }: { a, b -> b.name <=> a.name } )
        .map { meta, fasta ->
            def asm_meta = meta.assembly.subMap(['assembler','stage','id','build'])
            [ meta, asm_meta + (fasta.size() == 1 ? [ pri_fasta: fasta[0] ] : [ pri_fasta: fasta[0], alt_fasta: fasta[1] ] ) ]
        }
}

/*
Returns a new Map with entries merged.
For each key, if the value is a map, deep merge the subMaps
If the key is not a map, see if it exists in rhs and replace current value in lhs (same behaviour as Map + Map )
Add the missing keys from rhs to lhs map.
*/

def deepMergeMaps(Map lhs, Map rhs) {
    lhs.collectEntries { k, v -> rhs[k] instanceof Map ? [ (k): deepMergeMaps(v, rhs[k] ) ] : [ (k): rhs[k]?:v ] } + rhs.subMap(rhs.keySet()-lhs.keySet())
}
