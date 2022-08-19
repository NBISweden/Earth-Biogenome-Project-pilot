# Code adapted from https://github.com/marbl/merqury/blob/master/eval/false_duplications.sh
$1 == 1 {
    singlecopy_freq[$2] = $3 # Defer processing single copy k-mers until end when cutoff is fixed
    if ( max < $3 ){
        max  = $3
        cutoff = int(sprintf("%.0f",$2 * 1.5))
    }
}
$1 > 1 && $2 < cutoff {
    $1 == ">4" ? idx = 5 : idx = $1
    cp_sum[idx] += $NF
}
END {
    for (i = 1; i < cutoff; i++){
        cp_sum[1] += singlecopy_freq[i]
    }
    for (i = 2; i < 6; i++){
        dup += cp_sum[i]
    }
    OFS="\t"
    print "hist", "cutoff", "1", "2", "3", "4", ">4", "dup(>1)", "all", "dups%"
    print FILENAME, cutoff, cp_sum[1], cp_sum[2], cp_sum[3], cp_sum[4], cp_sum[5], dup, dup+cp_sum[1], (100*dup)/(dup+cp_sum[1])
}
