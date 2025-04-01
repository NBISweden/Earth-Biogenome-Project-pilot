/*
Custom method for Map class called deepMerge.

Returns a new Map with entries merged.
For each key, if the value is a map, deepMerge the subMaps
If the key is not a map, see if it exists in rhs and replace current value in lhs (same behaviour as Map + Map )
Add the missing keys from rhs to lhs map.
*/
Map.metaClass.deepMerge = { Map rhs ->
    def lhs = delegate
    lhs.collectEntries { k, v -> rhs[k] instanceof Map ? [ (k): v.deepMerge( rhs[k] ) ] : [ (k): rhs[k]?:v ] } + rhs.subMap(rhs.keySet()-lhs.keySet())
}

/* Custom method for String class called commonPrefix.

Returns a string of the longest common prefix between two strings
*/
String.metaClass.commonPrefix = { String last ->
    def first = delegate
    def min_length = Math.min( first.length(), last.length() )
    def prefix_length = 0
    while ( prefix_length < min_length && first.charAt(prefix_length) == last.charAt(prefix_length) ){
        prefix_length++
    }
    // return common prefix without trailing alphanumeric characters.
    first.substring(0, prefix_length).replaceAll(/[^a-zA-Z0-9]+$/, '')
}