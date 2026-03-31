#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = query_fasta
# $4 = n_top_hits (number of top hits to return per sequence)

N_TOP_HITS=$4

# replace each space in sequence headers of fasta files with the string "!?!?"
sed '/^>/ s/ /!?!?/g' $3 > query.fasta

# search
lastal \
    -P 1 \
    -m 10000 \
    -k 2 \
    -u 3 \
    -D 100 \
    -N 10000 \
    -f BlastTab \
    DB \
    query.fasta \
    > results_pre.tsv

# convert tsv back to original sequence names
sed 's/!?!?/ /g' results_pre.tsv > results.tsv
rm -f results_pre.tsv

# remove output lines starting with '#' (aka. convert to true .tsv), remove unneeded columns, then split by query into individual files
grep '^#' -v results.tsv | \
    cut -f1,2 |
    sed -e 's/!?!?/ /g' |
    awk  -v RS="\n" -v ORS="\n" -v FS="\t" -v OFS="\t" '
        BEGIN {
            prev_query = ""
            file_index = 0
        }
        { 
        if ( $1 != prev_query ) {
            file_index++
        }
        print $0 > file_index".split" 
        prev_query = $1
    }' - 

# loop through each one-query table, removing lines with targets in the same genus, and appending to combined file
touch results_filtered.tsv
for i in *.split; do 
    awk -v n_top_hits="$N_TOP_HITS" -v RS="\n" -v ORS="\n" -v FS="\t" -v OFS="\t" ' 
        BEGIN { hit_count = 1 }
        {
        # split query and target into arrays
        split($1, query_array, ";" )
        split($2, target_array, ";" )
        # get genus string of query
        q_genus = query_array[2] ";" query_array[3] ";" query_array[4] ";" query_array[5] ";" query_array[6] ";" query_array[7]
        # get genus string of target
        t_genus = target_array[2] ";" target_array[3] ";" target_array[4] ";" target_array[5] ";" target_array[6] ";" target_array[7]
        # only keep matches if:
        if (hit_count < n_top_hits && ((q_genus != t_genus && query_array[7] !~ /Unclassified/ && target_array[7] !~ /Unclassified/ ) || (query_array[6] != target_array[6]) || (query_array[5] != target_array[5]) || (query_array[4] != target_array[4]) || (query_array[3] != target_array[3]) || (query_array[2] != target_array[2]) ) ) print $1 "\t" $2;
        hit_count++
        }' \
        $i  \
    >> results_filtered.tsv
    rm -f $i
done 

rm -f results.tsv
rm -f *split
