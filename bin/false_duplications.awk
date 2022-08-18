# Code adapted from https://github.com/marbl/merqury/blob/master/eval/false_duplications.sh
$1==1 && max < $3 {
    max  = $3
    cutoff = sprintf("%.0f",$2 * 1.5)
}
$1 > 0 && $2 < cutoff {
    $1 > 4 ? idx = 5 : idx = $1
    cp_sum[idx] += $NF
}
END {
    for (i = 2; i < length(cp_sum); i++ ){
        dup += cp_sum[i]
    }
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "hist", "cutoff", "1", "2", "3", "4", ">4", "dup", "all", "%dups"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", FILENAME, cutoff, cp_sum[1], cp_sum[2], cp_sum[3], cp_sum[4], cp_sum[5], dup, dup+cp_sum[1], (100*dup)/(dup+cp_sum[1])
}
