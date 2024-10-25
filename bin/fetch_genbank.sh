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

BASENAME="${ACC_LIST##*/}"

NAME_NOEXT="${ACC_LIST%.*}"

## fetch sequences as .fasta file
epost \
	-input $ACC_LIST \
	-db nuccore \
	| \
    efetch \
	-format fasta \
	> ${NAME_NOEXT}.fasta

## fetch tab-sep accession + taxid 
cat $ACC_LIST \
    | \
    epost \
	-db nuccore \
    | \
    esummary \
	-db nuccore \
    | \
    xtract \
	-pattern DocumentSummary \
	-element Caption,TaxId \
	> ${NAME_NOEXT}.taxids.txt

if [ -s ${NAME_NOEXT}.taxids.txt ]; then 
	# if file is not empty
	exit 0
else 
	# if file is empty 
	exit 140
fi
