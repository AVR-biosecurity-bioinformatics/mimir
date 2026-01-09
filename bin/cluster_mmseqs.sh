#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = task memory in kilobytes
# $4 = fasta_files
# $5 = threshold csv
# $6 = type of clustering process, 'partial', 'large_genus', 'component'

# parse memory limit
TASK_MEMORY_KB=$3
SPLIT_MEMORY_LIMIT=$(( TASK_MEMORY_KB * 8 / 10 )) # 80% of total memory goes to mmseqs2

# determine the min sequence identity for the clusters
if [[ $6 == "partial" ]]; then
	# extract median genus identity from the thresholds .csv file
	MIN_SEQ_ID=$( awk -F, '(NR>1) && ($1=="genus")'  $5 | cut -f2 -d, )
elif [[ $6 == "large_genus" ]]; then
	# use median species identity
	MIN_SEQ_ID=$( awk -F, '(NR>1) && ($1=="species")'  $5 | cut -f2 -d, )
elif [[ $6 == "component" ]]; then
	# also use median genus identity
	MIN_SEQ_ID=$( awk -F, '(NR>1) && ($1=="genus")'  $5 | cut -f2 -d, )
else 
	echo "${6} is an incorrect value for the type of clustering process required"
fi


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
		--min-seq-id $MIN_SEQ_ID \
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
	rm tmp.tsv

done

# convert tsv back to original sequence names
sed 's/!?!?/ /g' clusters_pre.tsv > clusters.tsv

rm -f clusters_pre.tsv
rm -f renamed.fasta
rm -rf tmp