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
