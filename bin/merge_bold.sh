#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = fasta_list
# $4 = matching_taxids_list
# $5 = synchanges_list

### input Groovy lists to bash lists
FASTA_LIST=$(echo $3 | tr -d ',[]')

MATCHING_TAXIDS_LIST=$(echo $4 | tr -d ',[]')

SYNCHANGES_LIST=$(echo $5 | tr -d ',[]')

### concatenate .fasta files into a single file
cat $FASTA_LIST > bold_seqs.merged.fasta

### concatenate matching_taxids .csv files into a single file
awk 'FNR==1 && NR!=1{next;}{print}' $MATCHING_TAXIDS_LIST > matching_taxids.merged.csv

### concatenate synchanges .csv files into a single file
awk 'FNR==1 && NR!=1{next;}{print}' $SYNCHANGES_LIST > synchanges.merged.csv

