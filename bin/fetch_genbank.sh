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
epost \
	-input $ACC_LIST \
	-db nuccore \
	| \
    efetch \
	-format gb \
	> sequences.gb

# check number of records in .gb is the same as the number of accessions in accession file
if [[ $(grep -c "^//" sequences.gb) == $(cat accessions.txt | wc -l) ]]; then
	echo "Number of sequences fetched equals number of input accessions"
	exit 0
else 
	echo "Number of sequences fetched does not equal number of input accessions"
	sleep 10 # sleep 10 seconds to try to fend off too many queries
	exit 140 # to retry task
fi

# ## fetch sequences as .fasta file
# epost \
# 	-input $ACC_LIST \
# 	-db nuccore \
# 	| \
#     efetch \
# 	-format fasta \
# 	> ${NAME_NOEXT}.fasta

# ## fetch tab-sep accession + taxid 
# cat $ACC_LIST \
#     | \
#     epost \
# 	-db nuccore \
#     | \
#     esummary \
# 	-db nuccore \
#     | \
#     xtract \
# 	-pattern DocumentSummary \
# 	-element Caption,TaxId \
# 	> ${NAME_NOEXT}.taxids.txt

# if [ -s ${NAME_NOEXT}.taxids.txt ]; then 
# 	# if file is not empty

# 	# check number of .fasta sequences and number of rows in taxid file are the same
# 	N_SEQS=$( grep -c "^>" ${NAME_NOEXT}.fasta )
# 	N_TAXIDS=$( cat ${NAME_NOEXT}.taxids.txt | wc -l )

# 	if [[ $N_SEQS != $N_TAXIDS ]]; then
# 		echo "ERROR: number of .fasta sequences ($N_SEQS) and number of fetched taxids ($N_TAXIDS) is not equal"
# 		exit 140 # do this to retry process
# 	else 
# 		exit 0
# 	fi 

# else 
# 	# if file is empty 
# 	echo "ERROR: ${NAME_NOEXT}.taxids.txt is empty"
# 	exit 140 # do this to retry process
# fi
