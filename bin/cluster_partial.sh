#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = task memory in kilobytes
# $4 = fasta_files
# $5 = threshold csv

# parse memory limit
TASK_MEMORY_KB=$3
SPLIT_MEMORY_LIMIT=$(( TASK_MEMORY_KB * 9 / 10 )) # 90% of total memory should go to mmseqs2

# extract family_max threshold from the thresholds .csv file
FAMILY_MAX=$( awk -F, '(NR>1) && ($1=="family")'  $5 | cut -f4 -d, )

## loop through list of .fasta files, appending final clustering results to a single file
touch clusters_pre.tsv
for FILE in $4; do 
	# replace each space in sequence headers of fasta with the string "!?!?"
	sed '/^>/ s/ /!?!?/g' $FILE > renamed.fasta

	# create nucleotide DB
	mmseqs createdb \
		renamed.fasta \
		DB \
		--dbtype 2

	# clustering
	mmseqs cluster \
		DB \
		DB.clustered \
		tmp \
		--min-seq-id $FAMILY_MAX \
		--threads $2 \
		-s 7.5 \
		--split-memory-limit ${SPLIT_MEMORY_LIMIT}K

	# tsv output
	mmseqs createtsv \
		DB \
		DB \
		DB.clustered \
		tmp.tsv \
		--threads $2 

	cat tmp.tsv >> clusters_pre.tsv
	rm renamed.fasta
	rm DB*

done

# convert tsv back to original sequence names
sed 's/!?!?/ /g' clusters_pre.tsv > clusters.tsv

rm -f clusters_pre.tsv
rm -f renamed.fasta