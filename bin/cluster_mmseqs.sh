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
	echo "MIN_SEQ_ID will be determined per .fasta file"
	# get median for genus and species for calculations
	MED_GENUS=$( awk -F, '(NR>1) && ($1=="genus")'  $5 | cut -f2 -d, )
	MED_SPECIES=$( awk -F, '(NR>1) && ($1=="species")'  $5 | cut -f2 -d, )
	### HIGH_LIMIT=$MED_SPECIES
	HIGH_LIMIT="0.99"
	# minimum records to change min_seq_id 
	MIN_N=50000
	# records at which scaling stops at median species ID
	MAX_N=500000
elif [[ $6 == "large_genus" ]]; then
	# use median species identity
	MIN_SEQ_ID=$( awk -F, '(NR>1) && ($1=="species")'  $5 | cut -f2 -d, )
elif [[ $6 == "component" ]]; then
	# use median genus identity 
	MIN_SEQ_ID=$( awk -F, '(NR>1) && ($1=="genus")'  $5 | cut -f2 -d, )
else 
	echo "${6} is an incorrect value for the type of clustering process required"
fi


## loop through list of .fasta files, appending final clustering results to a single file
touch clusters_pre.tsv
for FILE in $4; do 
	# determine sequence identity based on number of sequences, for partial clustering only
	if [[ $6 == "partial" ]]; then
		# count number of records in file
		N_SEQS=$( grep -c ">" $4 )
		if [[ $N_SEQS < $MIN_N ]]; then 
			# use genus if few sequences
			MIN_SEQ_ID=$MED_GENUS
		elif [[ $N_SEQS > $MAX_N ]]; then
			# use species if many sequences
			MIN_SEQ_ID=$HIGH_LIMIT
		else 
			# scale between genus and high limit on log-scale if in-between
			MIN_SEQ_ID=$( awk \
				-v high_limit="$HIGH_LIMIT" \
				-v genus_id="$MED_GENUS" \
				-v n_seqs="$N_SEQS" \
				-v min_n="$MIN_N" \
				-v max_n="$MAX_N" \
				' BEGIN { 
					scaling_factor=(log(n_seqs) - log(min_n)) / (log(max_n) - log(min_n))
					print genus_id + ((high_limit - genus_id) * scaling_factor) 
				} ' )
			echo "Clustering threshold for file ${4} (${N_SEQS} sequences) set to $MIN_SEQ_ID"
		fi 
	fi
	
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