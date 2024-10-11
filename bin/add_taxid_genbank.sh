#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = fasta
# $4 = taxids

## process variables

FASTA=${3}

TAXIDS=${4}

BASENAME="${FASTA##*/}"

NAME_NOEXT="${FASTA%.*}"

## rename sequences using taxid key-value tsv 
# starts by removing everything but the sequence accession/version (eg. "LC823819.1")
# then replaces 

seqkit seq -i $FASTA \
    | \
seqkit replace \
	--pattern '^(.+)(\.\d)$' \
	--replacement '${1}${2}|{kv}' \
	--kv-file $TAXIDS \
	--line-width 0 \
    > ${NAME_NOEXT}.taxid.fasta