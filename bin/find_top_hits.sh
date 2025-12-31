#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = task memory in kilobytes
# $4 = query_fasta
# $5 = target_fasta
# $6 = n_top_hits (number of top hits to return per sequence)

# parse memory limit
TASK_MEMORY_KB=$3
SPLIT_MEMORY_LIMIT=$(( TASK_MEMORY_KB * 9 / 10 )) # 90% of total memory gos to mmseqs2

N_TOP_HITS=$6

# replace each space in sequence headers of fasta files with the string "!?!?"
sed '/^>/ s/ /!?!?/g' $4 > query.fasta
sed '/^>/ s/ /!?!?/g' $5 > target.fasta

# create databases for query and target fasta files
mmseqs createdb \
	query.fasta \
	queryDB \
	--dbtype 2

mmseqs createdb \
	target.fasta \
	targetDB \
	--dbtype 2

# compute prefiltering scores
mmseqs prefilter \
	queryDB \
	targetDB \
	pfDB \
	-s 7.5 \
	--max-seqs 1000 \
	--threads $2 \
	--split-memory-limit ${SPLIT_MEMORY_LIMIT}K

# align prefiltered sequences
mmseqs align \
	queryDB \
	targetDB \
	pfDB \
	alDB \
	--threads $2 

# convert alignment output to table
mmseqs convertalis \
	queryDB \
	targetDB \
	alDB \
	results_pre.tsv \
	--format-output "query,target"

# convert tsv back to original sequence names
sed 's/!?!?/ /g' results_pre.tsv > results.tsv

# remove unneeded files from searching
rm -f alDB*
rm -f pfDB*
rm -f query*
rm -f target*
rm -f results_pre.tsv
rm -rf tmp

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
        split($1, query_array, ";", seps )
        query_string=""
        for (i = 2; i <= 7; i++ )
            query_string = query_string ";" query_array[i]
        # get genus lineage string of target
        split($2, target_array, ";", seps )
        target_string=""
        for (i = 2; i <= 7; i++ )
            target_string = target_string ";" target_array[i]
        if ( query_string != target_string ){
            if ( hit_count < n_top_hits ){
                print $0
            }
        }
        hit_count++
    	}' \
		$i \
		>> results_filtered.tsv
done 

echo "Finished filtering hits"

# remove unneeded files
rm -f *.split
rm -f results.tsv