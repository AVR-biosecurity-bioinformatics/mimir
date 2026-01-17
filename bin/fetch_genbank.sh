#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = acc_list (txt file)
# $4 = entrez_key (either the key value or "no_key")

## process variables

if [[ "${4}" != "no_key" ]]; then
	export NCBI_API_KEY="${4}"
fi 

ACC_LIST=${3}

## fetch sequences as GenBank flat file

# pause between 0-10 sec
RAND_N=$(shuf -i 1-1000 -n1)
SLEEP_SEC=$(awk -v rand_n="$RAND_N" 'BEGIN { printf("%.2f", rand_n / 100) }')
sleep $SLEEP_SEC

efetch \
	-input "$ACC_LIST" \
	-db nuccore \
	-format gb \
	> sequences.gb

# check number of records in .gb is the same as the number of accessions in accession file
if [[ $(grep -c "^//" sequences.gb) == $(cat accessions.txt | wc -l) ]]; then
	echo "Number of sequences fetched equals number of input accessions"
	exit 0
else 
	echo "Number of sequences fetched does not equal number of input accessions -- will retry"
	exit 140 # to retry task with timeout code
fi
