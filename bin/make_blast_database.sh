#!/bin/bash
set -e
set -u
## args are the following:
# $1 = projectDir 
# $2 = threads
# $3 = target_fasta

# replace each space in sequence headers of fasta file with the string "!?!?"
sed '/^>/ s/ /!?!?/g' $3 > target.fasta

# create database
makeblastdb \
	-in target.fasta \
	-input_type fasta \
	-out blast_db \
	-dbtype nucl

# remove unneeded files
rm -f target.fasta
