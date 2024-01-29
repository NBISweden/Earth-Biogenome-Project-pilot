def dropMeta( channel ){
    channel.map{ it[0] }
}

def combineByMetaKeys( Map args = [:], lhs, rhs, keys ){
    lhs.map{ [ it[0].subMap(keys) ] + it }
        .combine( args + [ by: 0 ], rhs.map{ [ it[0].subMap(keys) ] + it.tail() } )
        .map { it.tail() }
}

def joinByMetaKeys( Map args = [:], lhs, rhs, keys ) {
    lhs.map{ [ it[0].subMap(keys) ] + it }
        .join( args + [ by: 0 ], rhs.map{ [ it[0].subMap(keys) ] + it.tail() } )
        .map{ it.tail() }
}
