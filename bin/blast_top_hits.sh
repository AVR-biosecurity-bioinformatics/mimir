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

# run blast
blastn \
	-query query.fasta \
	-db blast_db \
	-out results_pre.tsv \
	-strand plus \
	-task dc-megablast \
	-word_size 11 \
	-outfmt "6 qseqid sseqid" \
	-num_threads $2 \
	-max_target_seqs 1000

# convert tsv back to original sequence names
sed 's/!?!?/ /g' results_pre.tsv > results.tsv
rm -f results_pre.tsv

## filter results to only keep target hits outside of the query's genus

# split into files by query
### NOTE: This code assumes the results table is sorted by query and then by target hit significance, which is the default of mmseqs2
cat results.tsv | \
    awk '
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
    awk -v n_top_hits="$N_TOP_HITS" ' 
        BEGIN { hit_count = 1 }
        {
        # get genus lineage string of query
        split($1, query_array, ";" )
        query_string = query_array[2] ";" query_array[3] ";" query_array[4] ";" query_array[5] ";" query_array[6] ";" query_array[7]
        # get genus lineage string of target
        split($2, target_array, ";" )
        target_string = target_array[2] ";" target_array[3] ";" target_array[4] ";" target_array[5] ";" target_array[6] ";" target_array[7]
        if ( query_string != target_string && hit_count < n_top_hits ) print $0;
        hit_count++
        }' \
        $i \
    >> results_filtered.tsv
done 

echo "Finished filtering hits"

# remove unneeded files
rm -f *.split
rm -f results.tsv
rm -f query.fasta