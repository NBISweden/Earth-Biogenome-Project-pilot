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

def assembliesFromStage( assemblies, String stage ) {
    assemblies.filter { meta, assembly -> meta.assembly?.stage == stage }
}

def setAssemblyStage( assemblies, String stage ) {
    assemblies.map{ meta, assemblies -> [ meta.deepMerge([ assembly: [ stage: stage, build: "${meta.assembly.assembler}-${stage}-${meta.assembly.id}" ] ]), assemblies ] }
}
