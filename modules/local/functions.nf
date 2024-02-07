def dropMeta( channel ){
    channel.map{ it[1] }
}

def combineByMetaKeys( Map args = [:], lhs, rhs ){
    assert args.keySet != null
    switch(args.meta) {
        case 'lhs':
            // Return meta map from lhs
            return lhs.map{ [ it[0].subMap(args.keySet) ] + it }
                .combine( args + [ by: 0 ], rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
                .map { it.tail() }
        case 'rhs':
            // Return meta map from rhs
            def restructure = { key, clip, ... tail -> def list = tail.toList(); [ list.removeAt(clip) ] + list }
            return lhs.map{ [ it[0].subMap(args.keySet), it.size()-1 ] + it.tail() }
                .combine( args + [ by: 0 ], rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
                .map { restructure(*it) }
        case 'merge':
            // Return merged meta map
            def restructure = { key, clip, ... tail -> def list = tail.toList(); [ list.head().deepMerge(list.removeAt(clip)) ] + list.tail() }
            return lhs.map{ [ it[0].subMap(args.keySet), it.size() ] + it }
                .combine( args + [ by: 0 ], rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
                .map { restructure(*it) }
        case 'subset':
            // Return meta keys subset
            return lhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() }
                .combine( args + [ by: 0 ], rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
        default:
            assert args.meta in ['lhs','rhs','merge','subset']
    }
}

def joinByMetaKeys( Map args = [:], lhs, rhs ) {
    assert args.keySet != null
    switch(args.meta) {
        case 'lhs':
            // Return meta map from lhs
            return lhs.map{ [ it[0].subMap(args.keySet) ] + it }
                .join( args + [ by: 0 ], rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
                .map { it.tail() }
        case 'rhs':
            // Return meta map from rhs
            def restructure = { key, clip, ... tail -> def list = tail.toList(); [ list.removeAt(clip) ] + list }
            return lhs.map{ [ it[0].subMap(args.keySet), it.size()-1 ] + it.tail() }
                .join( args + [ by: 0 ], rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
                .map { restructure(*it) }
        case 'merge':
            // Return merged meta map
            def restructure = { key, clip, ... tail -> def list = tail.toList(); [ list.head().deepMerge(list.removeAt(clip)) ] + list.tail() }
            return lhs.map{ [ it[0].subMap(args.keySet), it.size() ] + it }
                .join( args + [ by: 0 ], rhs.map{ [ it[0].subMap(args.keySet) ] + it } )
                .map { restructure(*it) }
        case 'subset':
            // Return meta keys subset
            return lhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() }
                .join( args + [ by: 0 ], rhs.map{ [ it[0].subMap(args.keySet) ] + it.tail() } )
        default:
            assert args.meta in ['lhs','rhs','merge','subset']
    }
}
